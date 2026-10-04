"""
Tests that lbcproc.py honours its directory arguments.

Every test here runs with the cwd set to an empty directory that is distinct
from image_directory / raw_directory, then asserts that outputs land in
image_directory and that nothing is written to the cwd. image_directory is
passed without a trailing slash so that string concatenation of paths
(rather than os.path.join) would also be caught.
"""

import pytest
import numpy as np
from pathlib import Path
from astropy.io import fits
from ccdproc import ImageFileCollection

from conftest import (LBC_KEYWORDS, NX_SCIENCE, write_lbc_file,
                      write_master_flat)


# ---------------------------------------------------------------------------
# Directory fixtures: cwd, image_directory and raw_directory are all distinct
# ---------------------------------------------------------------------------

@pytest.fixture
def cwd_dir(tmp_path, monkeypatch):
    """An empty directory used as the cwd; must stay empty."""
    d = tmp_path / 'cwd'
    d.mkdir()
    monkeypatch.chdir(d)
    return d


@pytest.fixture
def out_dir(tmp_path):
    """The image_directory (output location)."""
    d = tmp_path / 'reduced'
    d.mkdir()
    return d


def assert_cwd_empty(cwd_dir):
    """Fail if anything has been written to the cwd."""
    leaked = sorted(p.name for p in cwd_dir.iterdir())
    assert leaked == [], f"Files written to cwd: {leaked}"


@pytest.fixture
def overscan_names(cwd_dir, out_dir, raw_dir, ic_with_flats):
    """Run go_overscan into out_dir; return the names it reports."""
    from LBCgo.lbcproc import go_overscan
    return go_overscan(ic_with_flats,
                       image_directory=str(out_dir),
                       raw_directory=str(raw_dir),
                       verbose=False,
                       return_files=True)


# ---------------------------------------------------------------------------
# go_overscan
# ---------------------------------------------------------------------------

def test_overscan_writes_to_image_directory(overscan_names, out_dir, cwd_dir):
    """go_overscan output must be written to image_directory, not the cwd."""
    assert len(overscan_names) == 1
    for name in overscan_names:
        assert (out_dir / name).is_file()
    assert_cwd_empty(cwd_dir)


# ---------------------------------------------------------------------------
# make_bias / go_bias
# ---------------------------------------------------------------------------

def test_make_bias_writes_to_image_directory(cwd_dir, out_dir, raw_dir,
                                             ic_bias_only):
    from LBCgo.lbcproc import make_bias
    make_bias(ic_bias_only, image_directory=str(out_dir),
              raw_directory=str(raw_dir), verbose=False)
    assert (out_dir / 'zero.fits').is_file()
    assert_cwd_empty(cwd_dir)


def test_go_bias_default_name_and_directories(overscan_names, out_dir,
                                              cwd_dir, raw_dir,
                                              three_bias_files):
    """go_bias with bias_file=None must find the 'zero.fits' that make_bias
    writes, read inputs from input_directory, and write to image_directory."""
    from LBCgo.lbcproc import make_bias, go_bias

    # Only the bias frames: one_object_file also lives in raw_dir
    ic_bias = ImageFileCollection(str(raw_dir), keywords=LBC_KEYWORDS,
                                  filenames=[p.name for p in three_bias_files])
    make_bias(ic_bias, image_directory=str(out_dir),
              raw_directory=str(raw_dir), verbose=False)

    zero_files = go_bias(overscan_names,
                         input_directory=str(out_dir),
                         bias_directory=str(out_dir),
                         image_directory=str(out_dir),
                         verbose=False, return_files=True)

    assert len(zero_files) == 1
    for name in zero_files:
        assert name.endswith('_zero.fits')
        assert (out_dir / name).is_file()
    assert_cwd_empty(cwd_dir)


# ---------------------------------------------------------------------------
# make_flatfield / go_flatfield
# ---------------------------------------------------------------------------

def test_make_flatfield_writes_to_image_directory(cwd_dir, out_dir, raw_dir,
                                                  ic_with_flats):
    from LBCgo.lbcproc import make_flatfield
    make_flatfield(ic_with_flats, filter_name='g-SLOAN',
                   image_directory=str(out_dir),
                   raw_directory=str(raw_dir), verbose=False)
    assert (out_dir / 'flat.g-SLOAN.fits').is_file()
    assert_cwd_empty(cwd_dir)


def test_go_flatfield_uses_directories(overscan_names, out_dir, cwd_dir):
    """go_flatfield must read from input_directory, write to image_directory,
    and move inputs into image_directory/data/ (not ./data/)."""
    from LBCgo.lbcproc import go_flatfield

    write_master_flat(out_dir, filter_name='g-SLOAN')
    ic_over = ImageFileCollection(str(out_dir), keywords=LBC_KEYWORDS,
                                  filenames=overscan_names)
    flat_files = go_flatfield(ic_over,
                              image_directory=str(out_dir),
                              input_directory=str(out_dir),
                              flat_directory=str(out_dir),
                              cosmiccorrect=False, verbose=False,
                              return_files=True)

    assert len(flat_files) == 1
    for name in flat_files:
        assert (out_dir / name).is_file()
    for name in overscan_names:
        assert not (out_dir / name).exists()
        assert (out_dir / 'data' / name).is_file()
    assert_cwd_empty(cwd_dir)


def test_go_flatfield_overwrites_existing_data_copy(overscan_names, out_dir,
                                                    cwd_dir):
    """A stale copy in data/ must be replaced, as 'mv' did previously."""
    from LBCgo.lbcproc import go_flatfield

    write_master_flat(out_dir, filter_name='g-SLOAN')
    (out_dir / 'data').mkdir()
    stale = out_dir / 'data' / overscan_names[0]
    stale.write_bytes(b'stale')

    ic_over = ImageFileCollection(str(out_dir), keywords=LBC_KEYWORDS,
                                  filenames=overscan_names)
    go_flatfield(ic_over, image_directory=str(out_dir),
                 input_directory=str(out_dir), flat_directory=str(out_dir),
                 cosmiccorrect=False, verbose=False)

    # The moved file is a valid FITS file, not the stale placeholder.
    with fits.open(stale) as hdul:
        assert len(hdul) == 5
    assert_cwd_empty(cwd_dir)


# ---------------------------------------------------------------------------
# make_targetdirectories / go_extractchips
# ---------------------------------------------------------------------------

def test_targetdirs_and_extractchips_use_image_directory(cwd_dir, out_dir):
    from LBCgo.lbcproc import make_targetdirectories, go_extractchips

    name = 'lbcb.20230101.000001_flat.fits'
    write_lbc_file(out_dir, name, imagetyp='object', filter_name='g-SLOAN',
                   object_name='NGC891', nx=NX_SCIENCE)
    ic = ImageFileCollection(str(out_dir), keywords=LBC_KEYWORDS)

    _, fltr_dirs = make_targetdirectories(ic, image_directory=str(out_dir),
                                          verbose=False)
    assert [Path(d) for d in fltr_dirs] == [out_dir / 'NGC891' / 'g-SLOAN']

    chip_files = go_extractchips(fltr_dirs, image_directory=str(out_dir),
                                 verbose=False, return_files=True)

    filter_dir = out_dir / 'NGC891' / 'g-SLOAN'
    assert len(chip_files) == 4
    for chip in range(1, 5):
        assert (filter_dir / f'lbcb.20230101.000001_{chip}.fits').is_file()
    # The multi-extension input goes to image_directory/data/
    assert (out_dir / 'data' / name).is_file()
    assert not (filter_dir / name).exists()
    assert_cwd_empty(cwd_dir)


def test_extractchips_accepts_dir_without_trailing_slash(cwd_dir, out_dir):
    """filter_directories without a trailing slash must still be globbed."""
    from LBCgo.lbcproc import go_extractchips

    filter_dir = out_dir / 'NGC891' / 'g-SLOAN'
    filter_dir.mkdir(parents=True)
    write_lbc_file(filter_dir, 'lbcb.20230101.000001_flat.fits',
                   imagetyp='object', filter_name='g-SLOAN',
                   object_name='NGC891', nx=NX_SCIENCE)

    chip_files = go_extractchips(str(filter_dir), lbc_chips=[1, 2],
                                 image_directory=str(out_dir),
                                 verbose=False, return_files=True)
    assert len(chip_files) == 2
    assert_cwd_empty(cwd_dir)


# ---------------------------------------------------------------------------
# lbcgo end to end (no astrometry)
# ---------------------------------------------------------------------------

def test_lbcgo_writes_only_to_image_directory(cwd_dir, out_dir, raw_dir,
                                              ic_with_flats, monkeypatch):
    """Run lbcgo (no astrometry) with an existing flat in image_directory.

    The flat must be found there and reused (make_flatfield's normalization
    box is hard-coded for full-size LBC chips, so it cannot run on the
    synthetic data), and every output must land under image_directory.
    """
    import LBCgo.lbcproc as lbcproc

    write_master_flat(out_dir, filter_name='g-SLOAN')

    def fail(*args, **kwargs):
        raise AssertionError('make_flatfield should not be called')
    monkeypatch.setattr(lbcproc, 'make_flatfield', fail)

    lbcproc.lbcgo(raw_directory=str(raw_dir), image_directory=str(out_dir),
                  do_astrometry=False, clean=False, verbose=False)

    filter_dir = out_dir / 'NGC891' / 'g-SLOAN'
    for chip in range(1, 5):
        assert (filter_dir / f'lbcb.20230101.000001_{chip}.fits').is_file()
    assert (out_dir / 'data' / 'lbcb.20230101.000001_over.fits').is_file()
    assert (out_dir / 'data' / 'lbcb.20230101.000001_flat.fits').is_file()
    assert_cwd_empty(cwd_dir)


# ---------------------------------------------------------------------------
# Mask/weight sidecars (LBCgo.masks) follow their images
# ---------------------------------------------------------------------------

def test_go_bias_copies_sidecar_to_image_directory(overscan_names, out_dir,
                                                   cwd_dir, raw_dir,
                                                   three_bias_files):
    from LBCgo.lbcproc import make_bias, go_bias
    from LBCgo.masks import sidecar_name

    assert (out_dir / sidecar_name(overscan_names[0], 'mask')).is_file()
    ic_bias = ImageFileCollection(str(raw_dir), keywords=LBC_KEYWORDS,
                                  filenames=[p.name for p in three_bias_files])
    make_bias(ic_bias, image_directory=str(out_dir),
              raw_directory=str(raw_dir), verbose=False)
    zero_files = go_bias(overscan_names, input_directory=str(out_dir),
                         bias_directory=str(out_dir),
                         image_directory=str(out_dir),
                         verbose=False, return_files=True)
    assert (out_dir / sidecar_name(zero_files[0], 'mask')).is_file()
    assert_cwd_empty(cwd_dir)


def test_lbcgo_sidecars_follow_images(cwd_dir, out_dir, raw_dir,
                                      ic_with_flats, monkeypatch):
    """Every image lbcgo leaves behind has its mask (and, after flat
    fielding, weight) sidecar beside it, under image_directory."""
    import LBCgo.lbcproc as lbcproc

    write_master_flat(out_dir, filter_name='g-SLOAN')

    def fail(*args, **kwargs):
        raise AssertionError('make_flatfield should not be called')
    monkeypatch.setattr(lbcproc, 'make_flatfield', fail)

    lbcproc.lbcgo(raw_directory=str(raw_dir), image_directory=str(out_dir),
                  do_astrometry=False, clean=False, verbose=False)

    filter_dir = out_dir / 'NGC891' / 'g-SLOAN'
    data_dir = out_dir / 'data'
    base = 'lbcb.20230101.000001'
    for chip in range(1, 5):
        assert (filter_dir / f'{base}_{chip}.mask.fits').is_file()
        assert (filter_dir / f'{base}_{chip}.weight.fits').is_file()
    assert (data_dir / f'{base}_over.mask.fits').is_file()
    assert (data_dir / f'{base}_flat.mask.fits').is_file()
    assert (data_dir / f'{base}_flat.weight.fits').is_file()
    # No sidecars stranded at the top level of image_directory
    assert sorted(p.name for p in out_dir.glob('*.mask.fits')) == []
    assert sorted(p.name for p in out_dir.glob('*.weight.fits')) == []
    assert_cwd_empty(cwd_dir)
