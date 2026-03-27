"""
Tests for go_sextractor(), go_scamp(), go_swarp(), and go_register()
in lbcregister.py.

All subprocess calls (Popen) are mocked — no external binaries required.
shutil.which is mocked to control binary-presence checks.
"""

import pytest
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock, call

from conftest import write_lbc_file, NX_SCIENCE, NY


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_mock_process():
    """A mock subprocess that completes immediately."""
    mock = MagicMock()
    mock.wait.return_value = 0
    return mock


def make_mock_votable():
    """Mock astropy votable result consumed by go_scamp after last iteration."""
    mock_vot = MagicMock()
    mock_array = MagicMock()
    # AstromSigma_Reference.data must be a 2D array
    mock_array.__getitem__ = MagicMock(
        return_value=MagicMock(data=np.array([[0.05, 0.05]])))
    mock_vot.get_first_table.return_value.array = mock_array
    return mock_vot


# ---------------------------------------------------------------------------
# go_sextractor — binary-not-found guard
# ---------------------------------------------------------------------------

def test_sextractor_raises_if_no_binary():
    """RuntimeError must be raised if 'sex' binary is absent."""
    from LBCgo.lbcregister import go_sextractor
    with patch('shutil.which', return_value=None):
        with pytest.raises(RuntimeError, match="SExtractor"):
            go_sextractor('dummy.fits')


# ---------------------------------------------------------------------------
# go_sextractor — Popen called with correct command
# ---------------------------------------------------------------------------

def test_sextractor_calls_popen(tmp_path):
    """go_sextractor must call Popen with a command starting with 'sex'."""
    from LBCgo.lbcregister import go_sextractor
    f = tmp_path / 'test_1.fits'
    f.touch()
    mock_proc = make_mock_process()
    with patch('shutil.which', return_value='/usr/bin/sex'), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp:
        go_sextractor(str(f), verbose=False)
    assert mp.called, "Popen was not called"
    cmd = mp.call_args[0][0]   # shlex.split result
    assert cmd[0] == 'sex', f"Expected command 'sex', got '{cmd[0]}'"


def test_sextractor_catalog_name_in_command(tmp_path):
    """The output catalog path (*_1.cat) must appear in the SExtractor command."""
    from LBCgo.lbcregister import go_sextractor
    f = tmp_path / 'test_1.fits'
    f.touch()
    mock_proc = make_mock_process()
    with patch('shutil.which', return_value='/usr/bin/sex'), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp:
        go_sextractor(str(f), verbose=False)
    cmd_str = ' '.join(mp.call_args[0][0])
    assert 'test_1.cat' in cmd_str, \
        f"Expected catalog 'test_1.cat' in command, got: {cmd_str}"


def test_sextractor_missing_paramfile_returns_none(tmp_path):
    """go_sextractor must return None if the parameter file does not exist."""
    from LBCgo.lbcregister import go_sextractor
    with patch('shutil.which', return_value='/usr/bin/sex'), \
         patch('os.path.exists', return_value=False):
        result = go_sextractor(str(tmp_path / 'test.fits'),
                               paramfile='/nonexistent/file.param',
                               verbose=False)
    assert result is None


# ---------------------------------------------------------------------------
# go_scamp — binary-not-found guard
# ---------------------------------------------------------------------------

def test_scamp_raises_if_no_binary():
    """RuntimeError must be raised if 'scamp' binary is absent."""
    from LBCgo.lbcregister import go_scamp
    with patch('shutil.which', return_value=None):
        with pytest.raises(RuntimeError, match="SCAMP"):
            go_scamp('dummy.fits')


# ---------------------------------------------------------------------------
# go_scamp — minimum iterations
# ---------------------------------------------------------------------------

def test_scamp_forces_min_two_iterations(tmp_path, capsys):
    """Passing num_iterations=1 must be silently upgraded to 2."""
    from LBCgo.lbcregister import go_scamp
    mock_proc = make_mock_process()
    with patch('shutil.which', return_value='/usr/bin/scamp'), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp, \
         patch('LBCgo.lbcregister.votable.parse',
               return_value=make_mock_votable()):
        go_scamp(str(tmp_path / 'test_1.fits'),
                 num_iterations=1, verbose=False)

    captured = capsys.readouterr()
    assert 'WARNING' in captured.out, \
        "Expected WARNING about minimum iterations"
    assert mp.call_count >= 2, \
        f"Expected at least 2 Popen calls, got {mp.call_count}"


def test_scamp_calls_popen_n_times(tmp_path):
    """Popen must be called exactly num_iterations times."""
    from LBCgo.lbcregister import go_scamp
    mock_proc = make_mock_process()
    with patch('shutil.which', return_value='/usr/bin/scamp'), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp, \
         patch('LBCgo.lbcregister.votable.parse',
               return_value=make_mock_votable()):
        go_scamp(str(tmp_path / 'test_1.fits'),
                 num_iterations=3, verbose=False)
    assert mp.call_count == 3, \
        f"Expected 3 Popen calls for 3 iterations, got {mp.call_count}"


def test_scamp_uses_gaia_dr3_by_default(tmp_path):
    """Default astrometric catalog must be GAIA-DR3."""
    from LBCgo.lbcregister import go_scamp
    mock_proc = make_mock_process()
    with patch('shutil.which', return_value='/usr/bin/scamp'), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp, \
         patch('LBCgo.lbcregister.votable.parse',
               return_value=make_mock_votable()):
        go_scamp(str(tmp_path / 'test_1.fits'), verbose=False)
    all_cmds = ' '.join(' '.join(c[0][0]) for c in mp.call_args_list)
    assert 'GAIA-DR3' in all_cmds, \
        f"Expected 'GAIA-DR3' in scamp commands, got: {all_cmds}"


# ---------------------------------------------------------------------------
# go_swarp — binary-not-found guard
# ---------------------------------------------------------------------------

def test_swarp_raises_if_no_binary(tmp_path):
    """RuntimeError must be raised if 'swarp' binary is absent."""
    from LBCgo.lbcregister import go_swarp
    f = write_lbc_file(tmp_path, 'img_1.fits',
                       imagetyp='object', filter_name='g-SLOAN',
                       nx=NX_SCIENCE)
    with patch('shutil.which', return_value=None):
        with pytest.raises(RuntimeError, match="SWarp"):
            go_swarp([str(f)])


# ---------------------------------------------------------------------------
# go_swarp — mixed-filter guard
# ---------------------------------------------------------------------------

def test_swarp_raises_on_mixed_filters(tmp_path):
    """ValueError must be raised when input files have different FILTER keywords."""
    from LBCgo.lbcregister import go_swarp
    f1 = write_lbc_file(tmp_path, 'img1_1.fits',
                        imagetyp='object', filter_name='g-SLOAN',
                        nx=NX_SCIENCE)
    f2 = write_lbc_file(tmp_path, 'img2_1.fits',
                        imagetyp='object', filter_name='r-SLOAN',
                        nx=NX_SCIENCE)
    with patch('shutil.which', return_value='/usr/bin/swarp'):
        with pytest.raises(ValueError, match="filter"):
            go_swarp([str(f1), str(f2)])


# ---------------------------------------------------------------------------
# go_swarp — output filename derived from header
# ---------------------------------------------------------------------------

def test_swarp_output_filename_from_header(tmp_path):
    """SWARP output filename must be derived from OBJECT and FILTER headers."""
    from LBCgo.lbcregister import go_swarp
    f1 = write_lbc_file(tmp_path, 'img1_1.fits',
                        imagetyp='object', filter_name='g-SLOAN',
                        object_name='NGC891', nx=NX_SCIENCE)
    mock_proc = make_mock_process()
    with patch('shutil.which', return_value='/usr/bin/swarp'), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp, \
         patch('astropy.io.fits.setval'):
        go_swarp([str(f1)], verbose=False)
    cmd_str = ' '.join(mp.call_args[0][0])
    assert 'NGC891' in cmd_str, \
        f"Expected object name 'NGC891' in swarp command: {cmd_str}"
    assert '.mos.fits' in cmd_str, \
        f"Expected '.mos.fits' in swarp command: {cmd_str}"


def test_swarp_weight_filename(tmp_path):
    """SWARP weight output filename must be *.mos.weight.fits."""
    from LBCgo.lbcregister import go_swarp
    f1 = write_lbc_file(tmp_path, 'img1_1.fits',
                        imagetyp='object', filter_name='g-SLOAN',
                        object_name='NGC891', nx=NX_SCIENCE)
    mock_proc = make_mock_process()
    with patch('shutil.which', return_value='/usr/bin/swarp'), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp, \
         patch('astropy.io.fits.setval'):
        go_swarp([str(f1)], verbose=False)
    cmd_str = ' '.join(mp.call_args[0][0])
    assert '.mos.weight.fits' in cmd_str, \
        f"Expected '.mos.weight.fits' in swarp command: {cmd_str}"


# ---------------------------------------------------------------------------
# go_register — chip selection
# ---------------------------------------------------------------------------

def test_register_lbc_chips_true_uses_all_four(tmp_path):
    """go_register with lbc_chips=True must run sextractor on all 4 chip files."""
    from LBCgo.lbcregister import go_register

    flt_dir = tmp_path / 'NGC891' / 'g-SLOAN'
    flt_dir.mkdir(parents=True)
    for chip in range(1, 5):
        write_lbc_file(flt_dir, f'lbcb.20230101.000001_{chip}.fits',
                       imagetyp='object', filter_name='g-SLOAN',
                       object_name='NGC891', nx=NX_SCIENCE)

    mock_proc = make_mock_process()
    with patch('shutil.which', return_value='/usr/bin/sex'), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp:
        go_register([str(flt_dir) + '/'],
                    lbc_chips=True,
                    do_sextractor=True,
                    do_scamp=False,
                    do_swarp=False,
                    verbose=False)

    # One Popen call per chip file for sextractor
    assert mp.call_count == 4, \
        f"Expected 4 sextractor calls (one per chip), got {mp.call_count}"


def test_register_subset_chips(tmp_path):
    """go_register with lbc_chips=[1,2] must run sextractor on only 2 files."""
    from LBCgo.lbcregister import go_register

    flt_dir = tmp_path / 'NGC891' / 'g-SLOAN'
    flt_dir.mkdir(parents=True)
    for chip in range(1, 5):
        write_lbc_file(flt_dir, f'lbcb.20230101.000001_{chip}.fits',
                       imagetyp='object', filter_name='g-SLOAN',
                       object_name='NGC891', nx=NX_SCIENCE)

    mock_proc = make_mock_process()
    with patch('shutil.which', return_value='/usr/bin/sex'), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp:
        go_register([str(flt_dir) + '/'],
                    lbc_chips=[1, 2],
                    do_sextractor=True,
                    do_scamp=False,
                    do_swarp=False,
                    verbose=False)

    assert mp.call_count == 2, \
        f"Expected 2 sextractor calls (chips 1 and 2), got {mp.call_count}"
