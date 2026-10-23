"""
Tests for go_scrub() in lbcproc.py (astromatic plan Stage A13).

go_scrub used to run ``rm data/`` (no ``-r``, so nothing was removed)
relative to the cwd. It now removes the intermediate files (``_over``,
``_zero``, ``_flat`` and their sidecars) from ``<image_directory>/data/``,
and ``data/`` itself only if that leaves it empty. ``data/`` may also be
the raw directory or a target directory, whose files must survive.
"""

import os

from conftest import write_lbc_file, write_master_flat

INTERMEDIATES = ['a_over.fits', 'a_over.mask.fits', 'a_zero.fits',
                 'a_flat.fits', 'a_flat.mask.fits', 'a_flat.weight.fits']


def _fill_data_dir(root):
    data = root / 'data'
    data.mkdir(parents=True)
    for name in INTERMEDIATES:
        (data / name).write_text('x')
    return data


def test_scrub_removes_intermediates_of_image_directory(tmp_path,
                                                        monkeypatch):
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


def test_scrub_keeps_other_files_and_subdirectories(tmp_path):
    """Only intermediates go; anything else keeps data/ in place."""
    from LBCgo.lbcproc import go_scrub
    data = _fill_data_dir(tmp_path)
    others = [data / 'lbcb.20230101.000001.fits',      # raw frame
              data / 'g-SLOAN' / 'lbcb.20230101.000001_1.fits',  # chip
              data / 'notes.txt']
    others[1].parent.mkdir()
    for p in others:
        p.write_text('x')

    go_scrub(image_directory=str(tmp_path), verbose=False)

    for name in INTERMEDIATES:
        assert not (data / name).exists()
    for p in others:
        assert p.is_file()


def test_scrub_without_data_dir_is_harmless(tmp_path):
    from LBCgo.lbcproc import go_scrub
    go_scrub(image_directory=str(tmp_path), verbose=False)
    assert list(tmp_path.iterdir()) == []


def test_scrub_leaves_symlinked_data_alone(tmp_path):
    from LBCgo.lbcproc import go_scrub
    target = tmp_path / 'elsewhere'
    target.mkdir()
    (target / 'a_over.fits').write_text('x')
    out = tmp_path / 'reduced'
    out.mkdir()
    os.symlink(target, out / 'data')

    go_scrub(image_directory=str(out), verbose=False)

    assert (out / 'data').is_symlink()
    assert (target / 'a_over.fits').is_file()


def test_lbcgo_clean_removes_data_dir(tmp_path, raw_dir, ic_with_flats):
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


def test_lbcgo_clean_keeps_raw_frames_in_data(tmp_path):
    """raw_directory = <image_directory>/data: the raw frames survive."""
    import LBCgo.lbcproc as lbcproc
    raw = tmp_path / 'data'
    raw.mkdir()
    raw_files = [write_lbc_file(raw, 'lbcb.20230101.000001.fits')]
    for i in range(2, 5):
        raw_files.append(write_lbc_file(
            raw, f'lbcb.20230101.00000{i}.fits', imagetyp='flat',
            object_name='SkyFlat', lbcobnam='SkyFlat', pixel_value=20000.0))

    lbcproc.lbcgo(raw_directory=str(raw), image_directory=str(tmp_path),
                  do_astrometry=False, clean=True, verbose=False)

    for p in raw_files:
        assert p.is_file()
    assert not list(raw.glob('*_over*.fits'))
    assert not list(raw.glob('*_flat*.fits'))
    for chip in range(1, 5):
        assert (tmp_path / 'NGC891' / 'g-SLOAN' /
                f'lbcb.20230101.000001_{chip}.fits').is_file()


def test_lbcgo_clean_keeps_target_named_data(tmp_path, raw_dir):
    """OBJECT = 'data': the final chip images in data/<filter>/ survive."""
    import LBCgo.lbcproc as lbcproc
    write_lbc_file(raw_dir, 'lbcb.20230101.000001.fits', object_name='data',
                   lbcobnam='data')
    out = tmp_path / 'reduced'
    out.mkdir()
    write_master_flat(out, filter_name='g-SLOAN')

    lbcproc.lbcgo(raw_directory=str(raw_dir), image_directory=str(out),
                  do_astrometry=False, clean=True, verbose=False)

    for chip in range(1, 5):
        assert (out / 'data' / 'g-SLOAN' /
                f'lbcb.20230101.000001_{chip}.fits').is_file()
    assert not list((out / 'data').glob('*_over*.fits'))
