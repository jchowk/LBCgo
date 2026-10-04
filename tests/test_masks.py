"""
Tests for LBCgo.masks and the mask/weight sidecars written by the pipeline.

go_overscan  -> <base>_over.mask.fits   (SATURATED from raw ADU)
go_flatfield -> <base>_flat.mask.fits, <base>_flat.weight.fits
go_extractchips -> <base>_<chip>.mask.fits, <base>_<chip>.weight.fits
"""

import numpy as np
import pytest
from pathlib import Path
from astropy.io import fits
from ccdproc import ImageFileCollection

from LBCgo import masks
from conftest import (write_lbc_file, write_master_flat, LBC_KEYWORDS,
                      N_CHIPS, NX_SCIENCE, NY)


# ---------------------------------------------------------------------------
# Unit tests: masks module
# ---------------------------------------------------------------------------

def test_sidecar_name():
    assert masks.sidecar_name('a/b_over.fits', 'mask') == 'a/b_over.mask.fits'
    assert masks.sidecar_name('b_1.fits', 'weight') == 'b_1.weight.fits'
    with pytest.raises(ValueError):
        masks.sidecar_name('b.fit', 'mask')


def test_saturation_mask_threshold_and_grow():
    data = np.zeros((20, 20))
    data[10, 10] = 60000.0
    data[2, 2] = 58000.0          # below 0.9 * 65535 = 58981.5
    m0 = masks.saturation_mask(data, 65535.0, fraction=0.9, grow=0)
    assert m0.sum() == 1 and m0[10, 10]
    m1 = masks.saturation_mask(data, 65535.0, fraction=0.9, grow=1)
    assert m1[9, 10] and m1[11, 10] and m1[10, 9] and m1[10, 11]
    assert not m1[2, 2]


def test_flat_mask_uniform_flat_is_clean():
    assert not masks.flat_mask(np.ones((50, 80))).any()


def test_flat_mask_flags_bad_column_and_pixel():
    flat = np.ones((50, 80))
    flat[:, 30] = 0.6          # low-response column
    flat[10, 60] = 1.5         # isolated hot pixel
    m = masks.flat_mask(flat, badpix_threshold=0.2)
    assert np.all(m[:, 30] & masks.BADPIX)
    assert m[10, 60] & masks.BADPIX
    # Neighbours are not flagged
    assert not (m[:, 29] & masks.BADPIX).any()
    assert not (m[:, 31] & masks.BADPIX).any()
    assert (m & masks.BADPIX).astype(bool).sum() == 50 + 1


def test_flat_mask_ignores_smooth_gradient():
    """A smooth 30% illumination gradient must not be called bad pixels."""
    yy, xx = np.mgrid[0:50, 0:80]
    flat = 0.85 + 0.3 * xx / 79.0
    assert not masks.flat_mask(flat).any()


def test_flat_mask_vignetting_and_nonpositive():
    flat = np.ones((50, 80))
    flat[:10, :10] = 0.3       # vignetted corner (smooth block, not BADPIX)
    flat[40, 40] = 0.0
    flat[41, 41] = np.nan
    m = masks.flat_mask(flat, vignette_threshold=0.5)
    assert np.all(m[:10, :10] & masks.VIGNETTED)
    assert not (m[2:8, 2:8] & masks.BADPIX).any()
    assert m[40, 40] & masks.BADPIX
    assert m[41, 41] & masks.BADPIX


def test_badpix_regions(tmp_path):
    f = tmp_path / 'bad.txt'
    f.write_text('# chip x1 x2 y1 y2\n2 5 5 1 50\n2 10 12 3 4  # block\n')
    regions = masks.read_badpix_regions(str(f))
    assert regions == {2: [(5, 5, 1, 50), (10, 12, 3, 4)]}
    m = masks.apply_badpix_regions(np.zeros((50, 80), np.uint8), regions[2])
    assert np.all(m[:, 4] == masks.BADPIX)
    assert np.all(m[2:4, 9:12] == masks.BADPIX)
    assert m.astype(bool).sum() == 50 + 6


def test_inverse_variance_weight_formula():
    flat = np.array([[1.0, 0.5, 2.0]])
    sky, gain, rn = 1000.0, 2.0, 10.0
    w = masks.inverse_variance_weight(flat, sky, gain, rn)
    expected = flat**2 / (sky * flat / gain + (rn / gain)**2)
    np.testing.assert_allclose(w, expected, rtol=1e-6)
    assert w.dtype == np.float32


def test_inverse_variance_weight_masked_and_bad_flat():
    flat = np.array([[1.0, 1.0, 0.0, np.nan]])
    mask = np.array([[0, 1, 0, 0]], np.uint8)
    w = masks.inverse_variance_weight(flat, 100.0, 1.75, 12.0, mask)
    assert w[0, 0] > 0
    assert np.all(w[0, 1:] == 0)


def test_inverse_variance_weight_matches_simulated_noise():
    """Weights must equal 1/variance of simulated flat-fielded pixels."""
    rng = np.random.default_rng(1)
    gain, rn, sky = 1.75, 12.0, 800.0      # sky in flat-fielded ADU
    n = 200000
    for f in (1.0, 0.6):
        electrons = rng.poisson(sky * f * gain, n) + rng.normal(0, rn, n)
        flattened = (electrons / gain) / f
        w = masks.inverse_variance_weight(np.array([f]), sky, gain, rn)[0]
        assert np.isclose(1.0 / flattened.var(), w, rtol=0.02)


def test_sky_level_ignores_masked_and_outliers():
    rng = np.random.default_rng(2)
    data = rng.normal(500.0, 5.0, (100, 100))
    data[:10, :10] = 1e5                   # masked bright region
    data[50, 50] = 1e6                     # unmasked outlier, clipped
    mask = np.zeros(data.shape, np.uint8)
    mask[:10, :10] = 1
    assert abs(masks.sky_level(data, mask) - 500.0) < 0.5


# ---------------------------------------------------------------------------
# Pipeline integration
# ---------------------------------------------------------------------------

def _saturated_raw(raw_dir, filename='lbcb.20230101.000001.fits'):
    """Object file whose chip 2 has a saturated pixel at (y=20, x=40)."""
    path = write_lbc_file(raw_dir, filename)
    with fits.open(path, mode='update') as hdul:
        hdul[2].data[20, 40] = 65000.0
        hdul[2].header['SATURATE'] = 65535
        hdul[2].header['GAIN'] = 2.0
        hdul[2].header['RDNOISE'] = 10.0
    return path


def _run_overscan(raw_dir, work_dir):
    """Run go_overscan; return full paths (it returns names relative to
    image_directory)."""
    from LBCgo.lbcproc import go_overscan
    ic = ImageFileCollection(str(raw_dir), keywords=LBC_KEYWORDS)
    names = go_overscan(ic, image_directory=str(work_dir) + '/',
                        raw_directory=str(raw_dir) + '/',
                        verbose=False, return_files=True)
    return [str(work_dir / name) for name in names]


def _run_flatfield(over_files, flat_value=20000.0, edit_flat=None, **kwargs):
    from LBCgo.lbcproc import go_flatfield
    over_dir = str(Path(over_files[0]).parent) + '/'
    flat = write_master_flat(over_dir, filter_name='g-SLOAN',
                             pixel_value=flat_value)
    if edit_flat is not None:
        with fits.open(flat, mode='update') as hdul:
            edit_flat(hdul)
    ic = ImageFileCollection(over_dir, keywords=LBC_KEYWORDS,
                             filenames=[Path(f).name for f in over_files])
    return go_flatfield(ic, flat_directory=over_dir, image_directory=over_dir,
                        input_directory=over_dir, cosmiccorrect=False,
                        verbose=False, return_files=True, **kwargs)


def test_overscan_writes_saturation_mask(raw_dir, work_dir):
    _saturated_raw(raw_dir)
    over = _run_overscan(raw_dir, work_dir)
    mask_file = masks.sidecar_name(over[0], 'mask')
    assert Path(mask_file).exists()
    with fits.open(mask_file) as hdul:
        assert len(hdul) == N_CHIPS + 1
        assert 'IMAGETYP' not in hdul[0].header
        assert hdul[1].data.shape == (NY, NX_SCIENCE)
        assert not hdul[1].data.any()
        m2 = hdul[2].data
        assert m2[20, 40] == masks.SATURATED
        assert m2[19, 40] == masks.SATURATED       # grown by 1 pixel
        assert m2.astype(bool).sum() == 5
        assert hdul[2].header['EXTNAME'] == 'chip2'


def test_overscan_mask_follows_trimsec_offset(raw_dir, work_dir):
    """Real LBC chips have a prescan (TRIMSEC [51:2098,...]); the mask must
    be trimmed exactly like the data."""
    path = write_lbc_file(raw_dir, 'lbcb.20230101.000001.fits')
    with fits.open(path, mode='update') as hdul:
        for ext in range(1, N_CHIPS + 1):
            hdul[ext].header['TRIMSEC'] = '[11:80,1:{0}]'.format(NY)
        hdul[1].data[5, 40] = 65000.0          # raw x=40 -> trimmed x=30
    over = _run_overscan(raw_dir, work_dir)
    with fits.open(over[0]) as dh, \
         fits.open(masks.sidecar_name(over[0], 'mask')) as mh:
        assert mh[1].data.shape == dh[1].data.shape == (NY, NX_SCIENCE - 10)
        assert mh[1].data[5, 30] == masks.SATURATED
        assert dh[1].data[5, 30] > 60000.0


def test_overscan_make_masks_false(raw_dir, work_dir):
    from LBCgo.lbcproc import go_overscan
    write_lbc_file(raw_dir, 'lbcb.20230101.000001.fits')
    ic = ImageFileCollection(str(raw_dir), keywords=LBC_KEYWORDS)
    over = go_overscan(ic, image_directory=str(work_dir) + '/',
                       raw_directory=str(raw_dir) + '/', verbose=False,
                       return_files=True, make_masks=False)
    assert (work_dir / over[0]).is_file()
    assert not (work_dir / masks.sidecar_name(over[0], 'mask')).exists()


def test_flatfield_writes_mask_and_weight(raw_dir, work_dir):
    _saturated_raw(raw_dir)
    over = _run_overscan(raw_dir, work_dir)

    def bad_flat(hdul):
        hdul[1].data[:, 30] *= 0.5          # bad column on chip 1
        hdul[3].data[:10, :10] *= 0.3       # vignetted corner on chip 3

    # An empty detector table, so the header GAIN/RDNOISE path is tested
    # (the packaged table has LBCB rows, and this file is named lbcb.*)
    from astropy.table import Table
    from LBCgo import detector
    empty = detector.write_detector_table(
        Table(names=detector.TABLE_COLUMNS,
              dtype=['U4', 'i4', 'f8', 'f8', 'f8', 'f8', 'U10']),
        str(work_dir / 'empty_detector.ecsv'))
    flat_files = _run_flatfield(over, edit_flat=bad_flat,
                                detector_table=empty)
    out = str(work_dir / flat_files[0])

    with fits.open(masks.sidecar_name(out, 'mask')) as mh, \
         fits.open(masks.sidecar_name(out, 'weight')) as wh:
        assert len(mh) == len(wh) == N_CHIPS + 1
        assert np.all(mh[1].data[:, 30] & masks.BADPIX)
        assert mh[2].data[20, 40] & masks.SATURATED
        assert np.all(mh[3].data[:10, :10] & masks.VIGNETTED)
        assert not mh[4].data.any()

        for ext in range(1, N_CHIPS + 1):
            w, m = wh[ext].data, mh[ext].data
            assert np.all(w[m != 0] == 0)
            assert np.all(w[m == 0] > 0)
            assert wh[ext].header['WGTTYPE'] == 'INVVAR'

        # Chip 2 header values propagate into the weight
        hdr = wh[2].header
        assert hdr['GAIN'] == 2.0 and hdr['RDNOISE'] == 10.0
        # Synthetic chips: 5000 ADU minus a 50 ADU overscan level
        assert abs(hdr['SKYLEVEL'] - 4950.0) < 20.0
        f = 1.0
        expected = f**2 / (hdr['SKYLEVEL'] * f / 2.0 + (10.0 / 2.0)**2)
        np.testing.assert_allclose(np.median(wh[2].data[wh[2].data > 0]),
                                   expected, rtol=1e-5)
        # Defaults used when headers lack GAIN/RDNOISE
        assert wh[1].header['GAIN'] == masks.DEFAULT_GAIN


def test_flatfield_moves_input_sidecar_to_data(raw_dir, work_dir):
    _saturated_raw(raw_dir)
    over = _run_overscan(raw_dir, work_dir)
    _run_flatfield(over)
    name = Path(masks.sidecar_name(over[0], 'mask')).name
    assert not (work_dir / name).exists()
    assert (work_dir / 'data' / name).exists()


def test_flatfield_without_overscan_mask(raw_dir, work_dir):
    """No saturation sidecar: still writes mask/weight from the flat."""
    _saturated_raw(raw_dir)
    over = _run_overscan(raw_dir, work_dir)
    Path(masks.sidecar_name(over[0], 'mask')).unlink()
    flat_files = _run_flatfield(over)
    out = str(work_dir / flat_files[0])
    with fits.open(masks.sidecar_name(out, 'mask')) as mh:
        assert not mh[2].data.any()


def test_flatfield_make_weights_false(raw_dir, work_dir):
    _saturated_raw(raw_dir)
    over = _run_overscan(raw_dir, work_dir)
    flat_files = _run_flatfield(over, make_weights=False)
    out = str(work_dir / flat_files[0])
    assert not Path(masks.sidecar_name(out, 'mask')).exists()
    assert not Path(masks.sidecar_name(out, 'weight')).exists()


def test_flatfield_badpix_file(raw_dir, work_dir, tmp_path):
    _saturated_raw(raw_dir)
    over = _run_overscan(raw_dir, work_dir)
    bpfile = tmp_path / 'badpix.txt'
    bpfile.write_text('4 7 7 1 50\n')
    flat_files = _run_flatfield(over, badpix_file=str(bpfile))
    out = str(work_dir / flat_files[0])
    with fits.open(masks.sidecar_name(out, 'mask')) as mh:
        assert np.all(mh[4].data[:, 6] == masks.BADPIX)
        assert mh[4].data.astype(bool).sum() == NY


def test_sidecars_excluded_from_object_collection(raw_dir, work_dir):
    _saturated_raw(raw_dir)
    over = _run_overscan(raw_dir, work_dir)
    _run_flatfield(over)
    ic = ImageFileCollection(str(work_dir), keywords=LBC_KEYWORDS)
    objects = ic.files_filtered(imagetyp='object')
    assert len(objects) == 1
    assert not any('.mask.' in f or '.weight.' in f for f in objects)


def test_extractchips_splits_sidecars(raw_dir, work_dir):
    from LBCgo.lbcproc import go_extractchips
    _saturated_raw(raw_dir)
    over = _run_overscan(raw_dir, work_dir)
    flat_files = _run_flatfield(over)
    flat_path = work_dir / flat_files[0]

    chip_files = go_extractchips(str(work_dir) + '/', verbose=False,
                                 return_files=True)
    assert len(chip_files) == N_CHIPS
    for chip, cf in enumerate(sorted(chip_files), start=1):
        with fits.open(masks.sidecar_name(cf, 'mask')) as mh, \
             fits.open(masks.sidecar_name(cf, 'weight')) as wh:
            assert len(mh) == len(wh) == 2
            assert wh[1].header['WGTTYPE'] == 'INVVAR'
            assert mh[1].header['EXTNAME'] == 'chip{0}'.format(chip)
            if chip == 2:
                assert mh[1].data[20, 40] & masks.SATURATED
                assert wh[1].data[20, 40] == 0

    # MEF sidecars follow the MEF file into data/
    assert (work_dir / 'data' / Path(masks.sidecar_name(str(flat_path),
                                                        'weight')).name).exists()


def test_targetdirectories_moves_sidecars(raw_dir, work_dir):
    from LBCgo.lbcproc import make_targetdirectories
    _saturated_raw(raw_dir)
    over = _run_overscan(raw_dir, work_dir)
    flat_files = _run_flatfield(over)
    ic = ImageFileCollection(str(work_dir), keywords=LBC_KEYWORDS,
                             filenames=flat_files)
    _, fltr_dirs = make_targetdirectories(ic,
                                          image_directory=str(work_dir) + '/',
                                          verbose=False)
    fdir = Path(fltr_dirs[0])
    for kind in ('mask', 'weight'):
        name = Path(masks.sidecar_name(flat_files[0], kind)).name
        assert (fdir / name).exists()
        assert not (work_dir / name).exists()


def test_bias_copies_sidecar(raw_dir, work_dir):
    from LBCgo.lbcproc import go_bias
    _saturated_raw(raw_dir)
    over = _run_overscan(raw_dir, work_dir)
    over_dir = str(work_dir) + '/'
    # Simple synthetic master bias with the trimmed chip shape
    hdul = fits.HDUList([fits.PrimaryHDU()])
    for _ in range(N_CHIPS):
        hdul.append(fits.ImageHDU(np.zeros((NY, NX_SCIENCE), 'f4'),
                                  header=fits.Header([('BUNIT', 'adu')])))
    hdul.writeto(over_dir + 'zero.fits')
    zero_files = go_bias([Path(f).name for f in over], bias_file='zero.fits',
                         image_directory=over_dir, input_directory=over_dir,
                         bias_directory=over_dir, verbose=False,
                         return_files=True)
    with fits.open(masks.sidecar_name(over_dir + zero_files[0], 'mask')) as mh:
        assert mh[2].data[20, 40] == masks.SATURATED


def test_flatfield_flags_cosmic_rays(raw_dir, work_dir):
    """With cosmiccorrect=True, pixels replaced by L.A.Cosmic get COSMIC."""
    path = _saturated_raw(raw_dir)
    with fits.open(path, mode='update') as hdul:
        hdul[4].data[25, 50] += 20000.0        # single-pixel hit on chip 4
    over = _run_overscan(raw_dir, work_dir)
    from LBCgo.lbcproc import go_flatfield
    over_dir = str(work_dir) + '/'
    write_master_flat(over_dir, filter_name='g-SLOAN')
    ic = ImageFileCollection(over_dir, keywords=LBC_KEYWORDS,
                             filenames=[Path(f).name for f in over])
    flat_files = go_flatfield(ic, flat_directory=over_dir,
                              image_directory=over_dir,
                              input_directory=over_dir, cosmiccorrect=True,
                              verbose=False, return_files=True)
    with fits.open(masks.sidecar_name(over_dir + flat_files[0], 'mask')) as mh, \
         fits.open(masks.sidecar_name(over_dir + flat_files[0], 'weight')) as wh:
        assert mh[4].data[25, 50] & masks.COSMIC
        assert wh[4].data[25, 50] == 0
        # The cleaner must not flag most of a quiet chip
        assert (mh[4].data & masks.COSMIC).astype(bool).mean() < 0.01
