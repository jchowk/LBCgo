"""
Tests for go_extractchips() in lbcproc.py.

go_extractchips() reads *_flat.fits MEF files and writes one single-chip
file per chip: *_1.fits, *_2.fits, *_3.fits, *_4.fits.
The input flat files are moved to a data/ subdirectory.
"""

import pytest
import numpy as np
from pathlib import Path
from astropy.io import fits

from conftest import N_CHIPS, NX_SCIENCE, NY, write_lbc_file


# ---------------------------------------------------------------------------
# Fixture: a filter subdirectory containing one flat-fielded MEF file
# ---------------------------------------------------------------------------

@pytest.fixture
def filter_dir(tmp_path):
    """Isolated filter directory — avoids sharing state with work_dir."""
    d = tmp_path / 'NGC891' / 'g-SLOAN'
    d.mkdir(parents=True)
    return d


@pytest.fixture(autouse=True)
def chdir_to_tmp(tmp_path, monkeypatch):
    """go_extractchips creates data/ relative to cwd; set cwd to tmp_path
    so each test gets its own isolated data/ directory."""
    monkeypatch.chdir(tmp_path)


@pytest.fixture
def flat_file(filter_dir):
    """One synthetic flat-fielded MEF file in filter_dir.

    Uses NX_SCIENCE columns (no overscan) to simulate post-flatfield data.
    """
    return write_lbc_file(filter_dir, 'lbcb.20230101.000001_flat.fits',
                          imagetyp='object', filter_name='g-SLOAN',
                          object_name='NGC891', nx=NX_SCIENCE)


# ---------------------------------------------------------------------------
# Output file count and naming
# ---------------------------------------------------------------------------

def test_extractchips_creates_four_files(filter_dir, flat_file):
    """go_extractchips() must create *_1.fits through *_4.fits."""
    from LBCgo.lbcproc import go_extractchips
    go_extractchips([str(filter_dir) + '/'], verbose=False)
    chip_files = list(filter_dir.glob('*_[1-4].fits'))
    assert len(chip_files) == N_CHIPS, \
        f"Expected {N_CHIPS} chip files, found {len(chip_files)}"


def test_extractchips_output_filename_convention(filter_dir, flat_file):
    """Output filenames follow the *_{chip}.fits pattern."""
    from LBCgo.lbcproc import go_extractchips
    go_extractchips([str(filter_dir) + '/'], verbose=False)
    for chip in range(1, N_CHIPS + 1):
        expected = filter_dir / f'lbcb.20230101.000001_{chip}.fits'
        assert expected.exists(), f"Expected {expected.name} to exist"


# ---------------------------------------------------------------------------
# Output file structure
# ---------------------------------------------------------------------------

def test_extractchips_single_extension_output(filter_dir, flat_file):
    """Each chip file must have exactly 2 HDUs: primary + 1 chip."""
    from LBCgo.lbcproc import go_extractchips
    go_extractchips([str(filter_dir) + '/'], verbose=False)
    hdul = fits.open(filter_dir / 'lbcb.20230101.000001_1.fits')
    assert len(hdul) == 2, \
        f"Expected 2 HDUs (primary + chip), got {len(hdul)}"
    hdul.close()


def test_extractchips_nextend_equals_one(filter_dir, flat_file):
    """Primary header NEXTEND must be 1 in each chip file."""
    from LBCgo.lbcproc import go_extractchips
    go_extractchips([str(filter_dir) + '/'], verbose=False)
    for chip in range(1, N_CHIPS + 1):
        hdul = fits.open(filter_dir / f'lbcb.20230101.000001_{chip}.fits')
        assert hdul[0].header['NEXTEND'] == 1, \
            f"Chip {chip}: expected NEXTEND=1, got {hdul[0].header['NEXTEND']}"
        hdul.close()


def test_extractchips_chip_data_matches_source(filter_dir, flat_file):
    """Chip 2 data in *_2.fits must match extension 2 in the source MEF."""
    from LBCgo.lbcproc import go_extractchips

    # Read chip 2 from source before extraction
    src = fits.open(flat_file)
    chip2_data = src[2].data.copy()
    src.close()

    go_extractchips([str(filter_dir) + '/'], verbose=False)

    out = fits.open(filter_dir / 'lbcb.20230101.000001_2.fits')
    assert np.allclose(out[1].data, chip2_data), \
        "Chip 2 data in output does not match source extension 2"
    out.close()


# ---------------------------------------------------------------------------
# Input file moved to data/
# ---------------------------------------------------------------------------

def test_extractchips_moves_flat_to_data_dir(filter_dir, flat_file):
    """The input *_flat.fits must be moved to a data/ directory."""
    from LBCgo.lbcproc import go_extractchips
    go_extractchips([str(filter_dir) + '/'], verbose=False)

    # File should no longer be in filter_dir
    assert not flat_file.exists(), \
        f"{flat_file.name} should have been moved out of filter_dir"

    # It should exist somewhere under a data/ directory
    found = list(Path('.').rglob('data/lbcb.20230101.000001_flat.fits'))
    assert len(found) > 0, \
        "Moved flat file not found under any data/ directory"


# ---------------------------------------------------------------------------
# Subset chips
# ---------------------------------------------------------------------------

def test_extractchips_subset_chips(filter_dir, flat_file):
    """lbc_chips=[1,3] must produce only *_1.fits and *_3.fits."""
    from LBCgo.lbcproc import go_extractchips
    go_extractchips([str(filter_dir) + '/'], lbc_chips=[1, 3], verbose=False)

    assert (filter_dir / 'lbcb.20230101.000001_1.fits').exists()
    assert (filter_dir / 'lbcb.20230101.000001_3.fits').exists()
    assert not (filter_dir / 'lbcb.20230101.000001_2.fits').exists()
    assert not (filter_dir / 'lbcb.20230101.000001_4.fits').exists()


# ---------------------------------------------------------------------------
# Return value and single-string input
# ---------------------------------------------------------------------------

def test_extractchips_returns_filelist(filter_dir, flat_file):
    """return_files=True must return a list of existing chip file paths."""
    from LBCgo.lbcproc import go_extractchips
    files = go_extractchips([str(filter_dir) + '/'],
                             verbose=False, return_files=True)
    assert isinstance(files, list)
    assert len(files) == N_CHIPS
    for f in files:
        assert Path(f).exists(), f"Expected {f} to exist"


def test_extractchips_single_dir_as_string(filter_dir, flat_file):
    """Passing a single string (not a list) must work without error."""
    from LBCgo.lbcproc import go_extractchips
    # Pass string, not list
    go_extractchips(str(filter_dir) + '/', verbose=False)
    assert (filter_dir / 'lbcb.20230101.000001_1.fits').exists()
