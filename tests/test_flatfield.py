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
# rather than calling make_flatfield(), because make_flatfield's normalization
# section (data[500:1500, 2000:2500]) is hardcoded for full LBC chip dimensions
# and returns an empty array on our small (50×80) synthetic data.

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
