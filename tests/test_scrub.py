"""
Tests for go_scrub() in lbcproc.py (astromatic plan Stage A13).

go_scrub used to run ``rm data/`` (no ``-r``, so nothing was removed)
relative to the cwd. It now removes ``<image_directory>/data/`` and
everything in it, and nothing else.
"""

import os

import pytest

from conftest import write_master_flat


def _fill_data_dir(root):
    data = root / 'data'
    (data / 'sub').mkdir(parents=True)
    for name in ['a_over.fits', 'a_flat.fits', 'a_flat.weight.fits']:
        (data / name).write_text('x')
    (data / 'sub' / 'nested.fits').write_text('x')
    return data


def test_scrub_removes_data_dir_of_image_directory(tmp_path, monkeypatch):
    from LBCgo.lbcproc import go_scrub
    out = tmp_path / 'reduced'
    out.mkdir()
    data = _fill_data_dir(out)
    keep = out / 'flat.g-SLOAN.fits'
    keep.write_text('x')
    # A data/ directory in the cwd is not touched
    cwd = tmp_path / 'cwd'
    cwd.mkdir()
    cwd_data = _fill_data_dir(cwd)
    monkeypatch.chdir(cwd)

    go_scrub(image_directory=str(out), verbose=False)

    assert not data.exists()
    assert keep.is_file()
    assert (cwd_data / 'a_over.fits').is_file()


def test_scrub_without_data_dir_is_harmless(tmp_path):
    from LBCgo.lbcproc import go_scrub
    go_scrub(image_directory=str(tmp_path), verbose=False)
    assert list(tmp_path.iterdir()) == []


def test_scrub_leaves_symlinked_data_alone(tmp_path):
    from LBCgo.lbcproc import go_scrub
    target = tmp_path / 'elsewhere'
    target.mkdir()
    (target / 'precious.fits').write_text('x')
    out = tmp_path / 'reduced'
    out.mkdir()
    os.symlink(target, out / 'data')

    go_scrub(image_directory=str(out), verbose=False)

    assert (out / 'data').is_symlink()
    assert (target / 'precious.fits').is_file()


def test_lbcgo_clean_removes_data_dir(tmp_path, raw_dir, ic_with_flats,
                                      monkeypatch):
    """lbcgo(clean=True) removes data/ under image_directory and keeps the
    per-chip images."""
    import LBCgo.lbcproc as lbcproc
    out = tmp_path / 'reduced'
    out.mkdir()
    write_master_flat(out, filter_name='g-SLOAN')

    lbcproc.lbcgo(raw_directory=str(raw_dir), image_directory=str(out),
                  do_astrometry=False, clean=True, verbose=False)

    assert not (out / 'data').exists()
    filter_dir = out / 'NGC891' / 'g-SLOAN'
    for chip in range(1, 5):
        assert (filter_dir / f'lbcb.20230101.000001_{chip}.fits').is_file()
