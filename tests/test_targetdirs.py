"""
Tests for make_targetdirectories() in lbcproc.py.

make_targetdirectories() reads an ImageFileCollection of flat-fielded files,
creates per-object and per-filter subdirectories, and moves files into them.
"""

import pytest
from pathlib import Path
from ccdproc import ImageFileCollection

from conftest import LBC_KEYWORDS, NX_SCIENCE, write_lbc_file


# ---------------------------------------------------------------------------
# Fixture: flat-fielded files written directly into work_dir
# ---------------------------------------------------------------------------

@pytest.fixture
def flat_files(work_dir):
    """Two synthetic flat-fielded object files in work_dir."""
    paths = []
    for i in range(1, 3):
        p = write_lbc_file(work_dir, f'lbcb.20230101.0000{i}_flat.fits',
                           imagetyp='object', filter_name='g-SLOAN',
                           object_name='NGC891', nx=NX_SCIENCE)
        paths.append(p)
    return paths


@pytest.fixture
def ic_flat_obj(work_dir, flat_files):
    """ImageFileCollection over work_dir containing the flat-fielded files."""
    return ImageFileCollection(str(work_dir), keywords=LBC_KEYWORDS)


# ---------------------------------------------------------------------------
# Directory creation
# ---------------------------------------------------------------------------

def test_targetdirs_creates_object_directory(ic_flat_obj, work_dir):
    """A directory named after the object must be created."""
    from LBCgo.lbcproc import make_targetdirectories
    make_targetdirectories(ic_flat_obj,
                           image_directory=str(work_dir) + '/',
                           verbose=False)
    assert (work_dir / 'NGC891').is_dir()


def test_targetdirs_creates_filter_subdirectory(ic_flat_obj, work_dir):
    """A filter subdirectory inside the object directory must be created."""
    from LBCgo.lbcproc import make_targetdirectories
    make_targetdirectories(ic_flat_obj,
                           image_directory=str(work_dir) + '/',
                           verbose=False)
    assert (work_dir / 'NGC891' / 'g-SLOAN').is_dir()


# ---------------------------------------------------------------------------
# File moves
# ---------------------------------------------------------------------------

def test_targetdirs_moves_files(ic_flat_obj, work_dir, flat_files):
    """Flat-fielded files must end up in object/filter subdirectory."""
    from LBCgo.lbcproc import make_targetdirectories
    make_targetdirectories(ic_flat_obj,
                           image_directory=str(work_dir) + '/',
                           verbose=False)
    dest = work_dir / 'NGC891' / 'g-SLOAN'
    moved = list(dest.glob('*_flat.fits'))
    assert len(moved) == len(flat_files), \
        f"Expected {len(flat_files)} files moved, found {len(moved)}"


def test_targetdirs_originals_removed_from_work_dir(ic_flat_obj, work_dir,
                                                     flat_files):
    """Original files must no longer exist in work_dir after moving."""
    from LBCgo.lbcproc import make_targetdirectories
    make_targetdirectories(ic_flat_obj,
                           image_directory=str(work_dir) + '/',
                           verbose=False)
    for f in flat_files:
        assert not f.exists(), \
            f"{f.name} should have been moved out of work_dir"


# ---------------------------------------------------------------------------
# Return values
# ---------------------------------------------------------------------------

def test_targetdirs_returns_correct_lists(ic_flat_obj, work_dir):
    """Must return (object_dirs, filter_dirs) with correct paths."""
    from LBCgo.lbcproc import make_targetdirectories
    obj_dirs, flt_dirs = make_targetdirectories(
        ic_flat_obj, image_directory=str(work_dir) + '/', verbose=False)
    assert any('NGC891' in d for d in obj_dirs), \
        f"NGC891 not found in object_dirs: {obj_dirs}"
    assert any('g-SLOAN' in d for d in flt_dirs), \
        f"g-SLOAN not found in filter_dirs: {flt_dirs}"


# ---------------------------------------------------------------------------
# Idempotency
# ---------------------------------------------------------------------------

def test_targetdirs_existing_directory_no_error(ic_flat_obj, work_dir):
    """Re-running on pre-existing directories must not raise."""
    from LBCgo.lbcproc import make_targetdirectories
    (work_dir / 'NGC891').mkdir(exist_ok=True)
    (work_dir / 'NGC891' / 'g-SLOAN').mkdir(exist_ok=True)
    # Should not raise even if directories already exist
    make_targetdirectories(ic_flat_obj,
                           image_directory=str(work_dir) + '/',
                           verbose=False)


# ---------------------------------------------------------------------------
# Space stripping from object name
# ---------------------------------------------------------------------------

def test_targetdirs_spaces_stripped_from_object_name(work_dir):
    """Object name 'NGC 891' (with space) must produce directory 'NGC891/'."""
    p = write_lbc_file(work_dir, 'lbcb.20230101.000099_flat.fits',
                       imagetyp='object', filter_name='g-SLOAN',
                       object_name='NGC 891', nx=NX_SCIENCE)
    ic = ImageFileCollection(str(work_dir), keywords=LBC_KEYWORDS,
                             filenames=['lbcb.20230101.000099_flat.fits'])
    from LBCgo.lbcproc import make_targetdirectories
    make_targetdirectories(ic, image_directory=str(work_dir) + '/',
                           verbose=False)
    assert (work_dir / 'NGC891').is_dir(), \
        "Expected directory 'NGC891', not found"
    assert not (work_dir / 'NGC 891').is_dir(), \
        "Directory 'NGC 891' (with space) should not exist"
