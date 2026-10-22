"""
Tests for make_flatfield() and go_flatfield() in lbcproc.py.

make_flatfield(): median-combines sky-flat frames per filter, writes flat.{filter}.fits.
go_flatfield():   divides object frames by the master flat.
"""

import pytest
import numpy as np
from pathlib import Path
from astropy.io import fits
from ccdproc import ImageFileCollection

from conftest import N_CHIPS, LBC_KEYWORDS, write_lbc_file, write_master_flat


# ---------------------------------------------------------------------------
# make_flatfield: output file creation
# ---------------------------------------------------------------------------

def test_make_flat_creates_output_file(ic_with_flats, work_dir, raw_dir):
    """make_flatfield() must write flat.g-SLOAN.fits."""
    from LBCgo.lbcproc import make_flatfield
    make_flatfield(ic_with_flats,
                   filter_name='g-SLOAN',
                   image_directory=str(work_dir) + '/',
                   raw_directory=str(raw_dir) + '/',
                   verbose=False)
    assert (work_dir / 'flat.g-SLOAN.fits').exists()


def test_make_flat_filename_uses_filter(work_dir, raw_dir, three_flat_files):
    """Output filename must embed the filter name."""
    from LBCgo.lbcproc import make_flatfield
    # Write r-SLOAN flats
    for i in range(1, 3):
        write_lbc_file(raw_dir, f'lbcr.20230101.00030{i}.fits',
                       imagetyp='flat', filter_name='r-SLOAN',
                       object_name='SkyFlat', lbcobnam='SkyFlat',
                       pixel_value=18000.0)
    ic = ImageFileCollection(str(raw_dir), keywords=LBC_KEYWORDS)
    make_flatfield(ic, filter_name='r-SLOAN',
                   image_directory=str(work_dir) + '/',
                   raw_directory=str(raw_dir) + '/',
                   verbose=False)
    assert (work_dir / 'flat.r-SLOAN.fits').exists()
    assert not (work_dir / 'flat.g-SLOAN.fits').exists()


# ---------------------------------------------------------------------------
# make_flatfield: no flat frames → return None
# ---------------------------------------------------------------------------

def test_make_flat_no_matching_flats_returns_none(ic_with_flats, work_dir,
                                                   raw_dir):
    """Requesting a filter with no flat frames must return None."""
    from LBCgo.lbcproc import make_flatfield
    result = make_flatfield(ic_with_flats,
                            filter_name='i-SLOAN',   # not in fixture
                            image_directory=str(work_dir) + '/',
                            raw_directory=str(raw_dir) + '/',
                            verbose=False)
    assert result is None


# ---------------------------------------------------------------------------
# make_flatfield: saturation rejection
# ---------------------------------------------------------------------------

def test_make_flat_skips_saturated_frames(work_dir, raw_dir):
    """Flat frames with per-chip average > 55000 ADU must be excluded.

    We write 2 unsaturated flats + 1 saturated flat.
    NCOMBINE in primary header should be < 3 (= 2 unsaturated frames used,
    one per chip so NCOMBINE reflects chip 0 count = 2).
    """
    from LBCgo.lbcproc import make_flatfield

    # 2 good flats
    for i in range(1, 3):
        write_lbc_file(raw_dir, f'lbcb.20230101.00040{i}.fits',
                       imagetyp='flat', filter_name='g-SLOAN',
                       object_name='SkyFlat', lbcobnam='SkyFlat',
                       pixel_value=20000.0)
    # 1 saturated flat
    write_lbc_file(raw_dir, 'lbcb.20230101.000403.fits',
                   imagetyp='flat', filter_name='g-SLOAN',
                   object_name='SkyFlat', lbcobnam='SkyFlat',
                   pixel_value=60000.0)

    ic = ImageFileCollection(str(raw_dir), keywords=LBC_KEYWORDS)
    make_flatfield(ic, filter_name='g-SLOAN',
                   image_directory=str(work_dir) + '/',
                   raw_directory=str(raw_dir) + '/',
                   verbose=False)

    hdul = fits.open(work_dir / 'flat.g-SLOAN.fits')
    ncombine = hdul[0].header.get('NCOMBINE', None)
    hdul.close()

    assert ncombine is not None, "NCOMBINE missing from primary header"
    assert ncombine < 3, \
        f"Expected NCOMBINE < 3 (saturated frame excluded), got {ncombine}"


# ---------------------------------------------------------------------------
# make_flatfield: normalization box follows the chip shape (Stage A4)
# ---------------------------------------------------------------------------

def test_flat_norm_section_full_chip():
    """Central half of a trimmed LBC chip (4608 rows x 2048 columns)."""
    from LBCgo.lbcproc import flat_norm_section
    rows, cols = flat_norm_section((4608, 2048))
    assert (rows.start, rows.stop) == (1152, 3456)
    assert (cols.start, cols.stop) == (512, 1536)


@pytest.mark.parametrize('shape', [(50, 80), (51, 79), (3, 2), (1, 1)])
def test_flat_norm_section_small_shapes(shape):
    """Non-empty, on the chip, and centred for any shape."""
    from LBCgo.lbcproc import flat_norm_section
    section = flat_norm_section(shape)
    box = np.zeros(shape)[section]
    assert box.size > 0
    for sl, n in zip(section, shape):
        assert 0 <= sl.start < sl.stop <= n
        # Centred to within one pixel
        assert abs(sl.start - (n - sl.stop)) <= 1


def test_make_flat_normalizes_frames_of_different_level(work_dir, raw_dir):
    """Flats at 10k, 20k and 40k ADU combine to a master normalized to 1
    over the central box, with the pixel pattern preserved.

    With the old hard-coded box the section was empty on these chips, so the
    scale factors (and the master flat) were NaN.
    """
    from LBCgo.headers import HeaderCollection
    from LBCgo.lbcproc import make_flatfield, flat_norm_section

    rng = np.random.default_rng(3)
    ny, nx = 50, 80
    # Smooth illumination pattern, asymmetric so that the central box and
    # the whole chip have different means (~1 %)
    pattern = (0.9 + 0.2 * (np.arange(nx) / (nx - 1))**2)[None, :] \
        * np.ones((ny, 1))
    pattern /= pattern[flat_norm_section(pattern.shape)].mean()
    noise = 1e-3   # fractional, per pixel
    for i, level in enumerate([10000., 20000., 40000.], start=1):
        p = write_lbc_file(raw_dir, f'lbcb.20230101.00050{i}.fits',
                           imagetyp='flat', filter_name='g-SLOAN',
                           object_name='SkyFlat', lbcobnam='SkyFlat',
                           pixel_value=level)
        with fits.open(p, mode='update') as hdul:
            for chip in range(1, N_CHIPS + 1):
                data = hdul[chip].data
                data[:, :nx] = level * pattern * \
                    (1 + noise * rng.standard_normal((ny, nx)))
                data[:, nx:] = 0.0    # flat overscan

    ic = HeaderCollection(str(raw_dir), keywords=LBC_KEYWORDS)
    make_flatfield(ic, filter_name='g-SLOAN',
                   image_directory=str(work_dir) + '/',
                   raw_directory=str(raw_dir) + '/', verbose=False)

    with fits.open(work_dir / 'flat.g-SLOAN.fits') as hdul:
        assert hdul[0].header['NCOMBINE'] == 3
        for chip in range(1, N_CHIPS + 1):
            master = hdul[chip].data.astype(float)
            assert np.all(np.isfinite(master))
            # Median of 3 frames: per-pixel error ~ 1.25 * noise; the box
            # mean of 1/x adds a bias of order noise**2.
            assert np.abs(master - pattern).max() < 6 * 1.25 * noise
            box = master[flat_norm_section(master.shape)]
            assert box.mean() == pytest.approx(1.0, abs=1e-3)


# ---------------------------------------------------------------------------
# make_flatfield: header keywords
# ---------------------------------------------------------------------------

def test_make_flat_removes_bzero(ic_with_flats, work_dir, raw_dir):
    """BZERO must be absent from chip headers in the master flat."""
    from LBCgo.lbcproc import make_flatfield
    make_flatfield(ic_with_flats, filter_name='g-SLOAN',
                   image_directory=str(work_dir) + '/',
                   raw_directory=str(raw_dir) + '/',
                   verbose=False)
    hdul = fits.open(work_dir / 'flat.g-SLOAN.fits')
    for ext in range(1, len(hdul)):
        assert 'BZERO' not in hdul[ext].header, \
            f"Ext {ext}: BZERO should be removed from master flat"
    hdul.close()


def test_make_flat_ncombine_in_primary_header(ic_with_flats, work_dir,
                                               raw_dir):
    """NCOMBINE must be set in the primary (ext 0) header."""
    from LBCgo.lbcproc import make_flatfield
    make_flatfield(ic_with_flats, filter_name='g-SLOAN',
                   image_directory=str(work_dir) + '/',
                   raw_directory=str(raw_dir) + '/',
                   verbose=False)
    hdul = fits.open(work_dir / 'flat.g-SLOAN.fits')
    assert 'NCOMBINE' in hdul[0].header, \
        "NCOMBINE missing from primary header of master flat"
    assert hdul[0].header['NCOMBINE'] == 3
    hdul.close()


# ---------------------------------------------------------------------------
# make_flatfield: data type and extensions
# ---------------------------------------------------------------------------

def test_make_flat_float32(ic_with_flats, work_dir, raw_dir):
    """Master flat chip data must be 32-bit float."""
    from LBCgo.lbcproc import make_flatfield
    make_flatfield(ic_with_flats, filter_name='g-SLOAN',
                   image_directory=str(work_dir) + '/',
                   raw_directory=str(raw_dir) + '/',
                   verbose=False)
    hdul = fits.open(work_dir / 'flat.g-SLOAN.fits')
    for ext in range(1, len(hdul)):
        dt = hdul[ext].data.dtype
        assert dt.kind == 'f' and dt.itemsize == 4, \
            f"Ext {ext}: expected 32-bit float, got {dt}"
    hdul.close()


def test_make_flat_four_chip_extensions(ic_with_flats, work_dir, raw_dir):
    """Master flat must have primary + 4 chip extensions."""
    from LBCgo.lbcproc import make_flatfield
    make_flatfield(ic_with_flats, filter_name='g-SLOAN',
                   image_directory=str(work_dir) + '/',
                   raw_directory=str(raw_dir) + '/',
                   verbose=False)
    hdul = fits.open(work_dir / 'flat.g-SLOAN.fits')
    assert len(hdul) == N_CHIPS + 1, \
        f"Expected {N_CHIPS + 1} HDUs, got {len(hdul)}"
    hdul.close()


# ---------------------------------------------------------------------------
# go_flatfield: output filename
# ---------------------------------------------------------------------------
# NOTE: These tests use write_master_flat() to build the master flat directly
# rather than calling make_flatfield(), so that the expected flat is exact.

def test_go_flat_output_filename(overscan_files):
    """go_flatfield() output must be named *_flat.fits."""
    from LBCgo.lbcproc import go_flatfield

    over_dir = str(Path(overscan_files[0]).parent) + '/'
    write_master_flat(over_dir, filter_name='g-SLOAN')

    ic_over = ImageFileCollection(over_dir, keywords=LBC_KEYWORDS,
                                  filenames=[Path(f).name
                                             for f in overscan_files])
    flat_files = go_flatfield(ic_over,
                              flat_directory=over_dir,
                              image_directory=over_dir,
                              input_directory=over_dir,
                              cosmiccorrect=False,
                              verbose=False,
                              return_files=True)

    assert flat_files is not None and len(flat_files) > 0
    for f in flat_files:
        assert '_flat.fits' in f, f"Expected *_flat.fits, got {f}"


def test_go_flat_returns_filelist(overscan_files):
    """go_flatfield() with return_files=True returns existing paths."""
    from LBCgo.lbcproc import go_flatfield

    over_dir = str(Path(overscan_files[0]).parent) + '/'
    write_master_flat(over_dir, filter_name='g-SLOAN')

    ic_over = ImageFileCollection(over_dir, keywords=LBC_KEYWORDS,
                                  filenames=[Path(f).name
                                             for f in overscan_files])
    flat_files = go_flatfield(ic_over,
                              flat_directory=over_dir,
                              image_directory=over_dir,
                              input_directory=over_dir,
                              cosmiccorrect=False,
                              verbose=False,
                              return_files=True)

    assert isinstance(flat_files, list)
    for f in flat_files:
        assert Path(over_dir + f).exists(), f"Expected {f} to exist"


def test_go_flat_moves_input_to_data_dir(overscan_files):
    """go_flatfield() must move the input _over.fits file to data/."""
    from LBCgo.lbcproc import go_flatfield

    over_dir = str(Path(overscan_files[0]).parent) + '/'
    write_master_flat(over_dir, filter_name='g-SLOAN')

    over_name = Path(overscan_files[0]).name
    ic_over = ImageFileCollection(over_dir, keywords=LBC_KEYWORDS,
                                  filenames=[over_name])
    go_flatfield(ic_over,
                 flat_directory=over_dir,
                 image_directory=over_dir,
                 input_directory=over_dir,
                 cosmiccorrect=False,
                 verbose=False)

    # Original overscan file should no longer be in over_dir
    assert not Path(over_dir + over_name).exists(), \
        f"{over_name} should have been moved to data/"
