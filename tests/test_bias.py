"""
Tests for make_bias() and go_bias() in lbcproc.py.

make_bias(): median-combines zero frames, writes zero.fits.
go_bias():   subtracts master bias from object frames.
"""

import pytest
import numpy as np
from pathlib import Path
from astropy.io import fits

from conftest import N_CHIPS


# ---------------------------------------------------------------------------
# make_bias: output file creation
# ---------------------------------------------------------------------------

def test_make_bias_creates_output_file(ic_bias_only, work_dir, raw_dir):
    """make_bias() must write zero.fits by default."""
    from LBCgo.lbcproc import make_bias
    make_bias(ic_bias_only,
              image_directory=str(work_dir) + '/',
              raw_directory=str(raw_dir) + '/',
              verbose=False)
    assert (work_dir / 'zero.fits').exists()


def test_make_bias_custom_filename(ic_bias_only, work_dir, raw_dir):
    """bias_filename parameter controls the output filename."""
    from LBCgo.lbcproc import make_bias
    make_bias(ic_bias_only,
              bias_filename='mybias.fits',
              image_directory=str(work_dir) + '/',
              raw_directory=str(raw_dir) + '/',
              verbose=False)
    assert (work_dir / 'mybias.fits').exists()
    assert not (work_dir / 'zero.fits').exists()


# ---------------------------------------------------------------------------
# make_bias: header content
# ---------------------------------------------------------------------------

def test_make_bias_ncombine_header(ic_bias_only, work_dir, raw_dir):
    """NCOMBINE in chip header must equal the number of combined frames."""
    from LBCgo.lbcproc import make_bias
    make_bias(ic_bias_only,
              image_directory=str(work_dir) + '/',
              raw_directory=str(raw_dir) + '/',
              verbose=False)
    hdul = fits.open(work_dir / 'zero.fits')
    # Fixture has 3 bias frames
    assert hdul[1].header['NCOMBINE'] == 3, \
        f"Expected NCOMBINE=3, got {hdul[1].header['NCOMBINE']}"
    hdul.close()


def test_make_bias_removes_biassec(ic_bias_only, work_dir, raw_dir):
    """BIASSEC must be absent from chip headers in the master bias."""
    from LBCgo.lbcproc import make_bias
    make_bias(ic_bias_only,
              image_directory=str(work_dir) + '/',
              raw_directory=str(raw_dir) + '/',
              verbose=False)
    hdul = fits.open(work_dir / 'zero.fits')
    for ext in range(1, len(hdul)):
        assert 'BIASSEC' not in hdul[ext].header, \
            f"Ext {ext}: BIASSEC should be removed from master bias"
    hdul.close()


# ---------------------------------------------------------------------------
# make_bias: data type and values
# ---------------------------------------------------------------------------

def test_make_bias_float32(ic_bias_only, work_dir, raw_dir):
    """Master bias chip data must be 32-bit float."""
    from LBCgo.lbcproc import make_bias
    make_bias(ic_bias_only,
              image_directory=str(work_dir) + '/',
              raw_directory=str(raw_dir) + '/',
              verbose=False)
    hdul = fits.open(work_dir / 'zero.fits')
    for ext in range(1, len(hdul)):
        dt = hdul[ext].data.dtype
        assert dt.kind == 'f' and dt.itemsize == 4, \
            f"Ext {ext}: expected 32-bit float, got {dt}"
    hdul.close()


def test_make_bias_values_near_input(ic_bias_only, work_dir, raw_dir):
    """Master bias values should be close to the input pixel_value (~500 ADU)."""
    from LBCgo.lbcproc import make_bias
    make_bias(ic_bias_only,
              image_directory=str(work_dir) + '/',
              raw_directory=str(raw_dir) + '/',
              verbose=False)
    hdul = fits.open(work_dir / 'zero.fits')
    median_val = float(np.median(hdul[1].data))
    assert 400.0 < median_val < 600.0, \
        f"Expected master bias ~500 ADU, got median {median_val:.1f}"
    hdul.close()


# ---------------------------------------------------------------------------
# make_bias: correct number of extensions
# ---------------------------------------------------------------------------

def test_make_bias_four_chip_extensions(ic_bias_only, work_dir, raw_dir):
    """Master bias must have primary + 4 chip extensions."""
    from LBCgo.lbcproc import make_bias
    make_bias(ic_bias_only,
              image_directory=str(work_dir) + '/',
              raw_directory=str(raw_dir) + '/',
              verbose=False)
    hdul = fits.open(work_dir / 'zero.fits')
    assert len(hdul) == N_CHIPS + 1, \
        f"Expected {N_CHIPS + 1} HDUs, got {len(hdul)}"
    hdul.close()


# ---------------------------------------------------------------------------
# go_bias: output filename
# ---------------------------------------------------------------------------

def test_go_bias_output_filename(overscan_files, raw_dir, ic_bias_only):
    """go_bias() output should be named *_zero.fits."""
    from LBCgo.lbcproc import make_bias, go_bias

    # The overscan files live in the work_dir used by that fixture;
    # derive its directory from the first returned path.
    over_dir = str(Path(overscan_files[0]).parent) + '/'

    # Build master bias into the same directory
    make_bias(ic_bias_only,
              image_directory=over_dir,
              raw_directory=str(raw_dir) + '/',
              verbose=False)

    zero_files = go_bias([Path(f).name for f in overscan_files],
                         bias_file='zero.fits',
                         bias_directory=over_dir,
                         input_directory=over_dir,
                         image_directory=over_dir,
                         verbose=False,
                         return_files=True)

    assert zero_files is not None and len(zero_files) > 0
    for f in zero_files:
        assert '_zero.fits' in f, f"Expected *_zero.fits, got {f}"


# ---------------------------------------------------------------------------
# go_bias: subtraction correctness
# ---------------------------------------------------------------------------

def test_go_bias_subtracts_correctly(overscan_files, raw_dir, ic_bias_only):
    """Output pixel values should be approximately input − bias level."""
    from LBCgo.lbcproc import make_bias, go_bias

    over_dir = str(Path(overscan_files[0]).parent) + '/'

    make_bias(ic_bias_only,
              image_directory=over_dir,
              raw_directory=str(raw_dir) + '/',
              verbose=False)

    # Record the median of the overscan'd object before bias subtraction
    hdul_before = fits.open(overscan_files[0])
    median_before = float(np.median(hdul_before[1].data))
    hdul_before.close()

    zero_files = go_bias([Path(f).name for f in overscan_files],
                         bias_file='zero.fits',
                         bias_directory=over_dir,
                         input_directory=over_dir,
                         image_directory=over_dir,
                         verbose=False,
                         return_files=True)

    hdul_after = fits.open(Path(overscan_files[0]).parent / zero_files[0])
    median_after = float(np.median(hdul_after[1].data))
    hdul_after.close()

    # Object is ~5000 ADU, bias is ~500 ADU → result should be ~4500 ADU
    expected = median_before - 500.0
    assert abs(median_after - expected) < 200.0, \
        f"After bias sub: expected ~{expected:.0f}, got {median_after:.0f}"
