"""
Tests for go_overscan() in lbcproc.py.

Synthetic LBC data: NX=100 cols, NY=50 rows per chip.
  BIASSEC = [81:100,1:50]  →  20 overscan columns
  TRIMSEC = [1:80,1:50]    →  80 science columns
After trim: output shape = (50, 80).
"""

import pytest
import numpy as np
from pathlib import Path
from astropy.io import fits

from conftest import NX_SCIENCE, NY, N_CHIPS


# ---------------------------------------------------------------------------
# Output file naming
# ---------------------------------------------------------------------------

def test_overscan_output_filename(overscan_files):
    """Output files should be named *_over.fits."""
    assert len(overscan_files) > 0
    for f in overscan_files:
        assert '_over.fits' in f


def test_overscan_creates_files(overscan_files, work_dir):
    """Output files must exist on disk."""
    for f in overscan_files:
        assert Path(f).exists(), f"Expected {f} to exist"


# ---------------------------------------------------------------------------
# Science-region shape after overscan trim
# ---------------------------------------------------------------------------

def test_overscan_output_shape(overscan_files):
    """Each chip should have NX_SCIENCE=80 columns after trim."""
    hdul = fits.open(overscan_files[0])
    for ext in range(1, len(hdul)):
        ny, nx = hdul[ext].data.shape
        assert nx == NX_SCIENCE, \
            f"Ext {ext}: expected {NX_SCIENCE} cols after trim, got {nx}"
        assert ny == NY, \
            f"Ext {ext}: expected {NY} rows, got {ny}"
    hdul.close()


# ---------------------------------------------------------------------------
# Header keyword removal
# ---------------------------------------------------------------------------

def test_overscan_removes_biassec(overscan_files):
    """BIASSEC must be absent from all chip headers after overscan."""
    hdul = fits.open(overscan_files[0])
    for ext in range(1, len(hdul)):
        assert 'BIASSEC' not in hdul[ext].header, \
            f"Ext {ext}: BIASSEC should be removed"
    hdul.close()


def test_overscan_removes_trimsec(overscan_files):
    """TRIMSEC must be absent from all chip headers after overscan."""
    hdul = fits.open(overscan_files[0])
    for ext in range(1, len(hdul)):
        assert 'TRIMSEC' not in hdul[ext].header, \
            f"Ext {ext}: TRIMSEC should be removed"
    hdul.close()


def test_overscan_removes_datasec(overscan_files):
    """DATASEC must be absent from all chip headers after overscan."""
    hdul = fits.open(overscan_files[0])
    for ext in range(1, len(hdul)):
        assert 'DATASEC' not in hdul[ext].header, \
            f"Ext {ext}: DATASEC should be removed"
    hdul.close()


def test_overscan_removes_rotangle(overscan_files):
    """ROTANGLE must be absent from all chip headers after overscan."""
    hdul = fits.open(overscan_files[0])
    for ext in range(1, len(hdul)):
        assert 'ROTANGLE' not in hdul[ext].header, \
            f"Ext {ext}: ROTANGLE should be removed"
    hdul.close()


def test_overscan_removes_cstar_keywords(overscan_files):
    """Keywords matching C*A pattern (e.g. CTYPE1A) must be removed."""
    hdul = fits.open(overscan_files[0])
    for ext in range(1, len(hdul)):
        for key in hdul[ext].header.keys():
            is_cstar_a = (len(key) > 2
                          and key.startswith('C')
                          and key.endswith('A'))
            assert not is_cstar_a, \
                f"Ext {ext}: C*A keyword '{key}' should have been removed"
    hdul.close()


# ---------------------------------------------------------------------------
# Data type
# ---------------------------------------------------------------------------

def test_overscan_float32_output(overscan_files):
    """Output chip data must be 32-bit float (FITS stores big-endian >f4)."""
    hdul = fits.open(overscan_files[0])
    for ext in range(1, len(hdul)):
        dt = hdul[ext].data.dtype
        assert dt.kind == 'f' and dt.itemsize == 4, \
            f"Ext {ext}: expected 32-bit float, got {dt}"
    hdul.close()


# ---------------------------------------------------------------------------
# Return value
# ---------------------------------------------------------------------------

def test_overscan_returns_filelist(overscan_files):
    """With return_files=True, a non-empty list of paths is returned."""
    assert isinstance(overscan_files, list)
    assert len(overscan_files) > 0


# ---------------------------------------------------------------------------
# Chip count
# ---------------------------------------------------------------------------

def test_overscan_all_chips_default(overscan_files):
    """Default lbc_chips=True → 4 science extensions (primary + 4)."""
    hdul = fits.open(overscan_files[0])
    assert len(hdul) == N_CHIPS + 1, \
        f"Expected {N_CHIPS + 1} HDUs (primary + {N_CHIPS} chips), got {len(hdul)}"
    hdul.close()


def test_overscan_subset_chips(ic_with_flats, work_dir, raw_dir):
    """lbc_chips=[1,3] → only 2 science extensions in output."""
    from LBCgo.lbcproc import go_overscan
    files = go_overscan(ic_with_flats,
                        lbc_chips=[1, 3],
                        image_directory=str(work_dir) + '/',
                        raw_directory=str(raw_dir) + '/',
                        verbose=False,
                        return_files=True)
    assert len(files) > 0
    hdul = fits.open(files[0])
    assert len(hdul) == 3, \
        f"Expected 3 HDUs (primary + chips 1,3), got {len(hdul)}"
    hdul.close()


# ---------------------------------------------------------------------------
# objects_only filtering
# ---------------------------------------------------------------------------

def test_overscan_objects_only_skips_flats(ic_with_flats, work_dir, raw_dir):
    """objects_only=True (default) must skip flat frames.

    Fixture has 1 object + 3 flats.  Only the object should be overscan'd.
    """
    from LBCgo.lbcproc import go_overscan
    files = go_overscan(ic_with_flats,
                        objects_only=True,
                        image_directory=str(work_dir) + '/',
                        raw_directory=str(raw_dir) + '/',
                        verbose=False,
                        return_files=True)
    assert len(files) == 1, \
        f"Expected 1 overscan'd object file, got {len(files)}"
