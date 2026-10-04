"""
Tests for go_sextractor(), go_scamp(), go_swarp(), and go_register()
in lbcregister.py.

All subprocess calls (Popen) are mocked — no external binaries required.
shutil.which is mocked to control binary-presence checks.
"""

import pytest
import numpy as np
from astropy.table import Table
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
    conf = tmp_path / 'swarp.conf'
    conf.write_text('')
    with patch('shutil.which', return_value='/usr/bin/swarp'), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp, \
         patch('astropy.io.fits.setval'):
        go_swarp([str(f1)], configfile=str(conf), verbose=False)
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
    conf = tmp_path / 'swarp.conf'
    conf.write_text('')
    with patch('shutil.which', return_value='/usr/bin/swarp'), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp, \
         patch('astropy.io.fits.setval'):
        go_swarp([str(f1)], configfile=str(conf), verbose=False)
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


# ---------------------------------------------------------------------------
# go_sextractor — Phase 0 §5.3: weights, flags, tunables, tool names
# ---------------------------------------------------------------------------

def _run_sextractor(tmp_path, which=lambda n: '/usr/bin/' + n, **kwargs):
    """Run go_sextractor on tmp_path/test_1.fits with Popen mocked; return cmd."""
    from LBCgo.lbcregister import go_sextractor
    f = tmp_path / 'test_1.fits'
    f.touch()
    mock_proc = make_mock_process()
    with patch('shutil.which', side_effect=which), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp:
        go_sextractor(str(f), verbose=False, **kwargs)
    return mp.call_args[0][0]


def _opt(cmd, flag):
    """Value following ``flag`` in a command list."""
    return cmd[cmd.index(flag) + 1]


def test_sextractor_alignment_defaults(tmp_path):
    """Defaults: thresh 5/5 (no stray 8), mesh 32, filter 3, mincont 1e-4."""
    cmd = _run_sextractor(tmp_path)
    assert _opt(cmd, '-DETECT_THRESH') == '5.0'
    assert _opt(cmd, '-ANALYSIS_THRESH') == '5.0'
    assert _opt(cmd, '-BACK_SIZE') == '32'
    assert _opt(cmd, '-BACK_FILTERSIZE') == '3'
    assert float(_opt(cmd, '-DEBLEND_MINCONT')) == 1e-4


def test_sextractor_tunables_are_arguments(tmp_path):
    cmd = _run_sextractor(tmp_path, detect_thresh=3, analysis_thresh=2,
                          back_size=64, back_filtersize=5,
                          deblend_mincont=0.005)
    assert _opt(cmd, '-DETECT_THRESH') == '3'
    assert _opt(cmd, '-ANALYSIS_THRESH') == '2'
    assert _opt(cmd, '-BACK_SIZE') == '64'
    assert _opt(cmd, '-BACK_FILTERSIZE') == '5'
    assert float(_opt(cmd, '-DEBLEND_MINCONT')) == 0.005


def test_sextractor_uses_weight_and_flag_sidecars(tmp_path):
    for kind in ('weight', 'mask'):
        (tmp_path / f'test_1.{kind}.fits').touch()
    cmd = _run_sextractor(tmp_path)
    assert _opt(cmd, '-WEIGHT_TYPE') == 'MAP_WEIGHT'
    assert _opt(cmd, '-WEIGHT_IMAGE').endswith('test_1.weight.fits')
    assert _opt(cmd, '-FLAG_IMAGE').endswith('test_1.mask.fits')


def test_sextractor_no_sidecars_runs_unweighted(tmp_path):
    cmd = _run_sextractor(tmp_path)
    assert '-WEIGHT_TYPE' not in cmd
    assert '-FLAG_IMAGE' not in cmd


def test_sextractor_sidecars_can_be_disabled(tmp_path):
    for kind in ('weight', 'mask'):
        (tmp_path / f'test_1.{kind}.fits').touch()
    cmd = _run_sextractor(tmp_path, use_weight=False, use_flags=False)
    assert '-WEIGHT_TYPE' not in cmd
    assert '-FLAG_IMAGE' not in cmd


def test_sextractor_flag_columns_added_only_with_flag_image(tmp_path):
    """IMAFLAGS_ISO must be requested iff a FLAG_IMAGE is passed (SExtractor
    fails otherwise)."""
    (tmp_path / 'test_1.mask.fits').touch()
    captured = {}

    def grab(cmd, **kw):
        with open(_opt(cmd, '-PARAMETERS_NAME')) as fh:
            captured['param'] = fh.read()
        return make_mock_process()

    from LBCgo.lbcregister import go_sextractor
    f = tmp_path / 'test_1.fits'
    f.touch()
    with patch('shutil.which', return_value='/usr/bin/sex'), \
         patch('LBCgo.lbcregister.Popen', side_effect=grab):
        go_sextractor(str(f), verbose=False)
    lines = [ln.split('#')[0].strip() for ln in captured['param'].splitlines()]
    assert 'IMAFLAGS_ISO(1)' in lines and 'NIMAFLAGS_ISO(1)' in lines


def test_sextractor_subtracted_image_mode(tmp_path):
    """Extended-target mode runs on the subtracted image with BACK_TYPE MANUAL
    but still names the catalog after the original chip file."""
    sub = tmp_path / 'test_1.sub.fits'
    cmd = _run_sextractor(tmp_path, subtracted_image=str(sub))
    assert Path(cmd[1]).name == sub.name  # may be a staged symlink
    assert _opt(cmd, '-BACK_TYPE') == 'MANUAL'
    assert _opt(cmd, '-BACK_VALUE') == '0'
    assert _opt(cmd, '-CATALOG_NAME').endswith('test_1.cat')


def test_sextractor_path_with_spaces_is_staged(tmp_path):
    """SExtractor cannot read paths with spaces; inputs are symlinked into a
    space-free directory and the catalog is moved back."""
    from LBCgo.lbcregister import go_sextractor
    d = tmp_path / 'with space'
    d.mkdir()
    f = d / 'test_1.fits'
    f.touch()
    seen = {}

    def fake_sex(cmd, **kw):
        assert all(' ' not in a for a in cmd), cmd
        seen['image'] = cmd[1]
        Path(_opt(cmd, '-CATALOG_NAME')).write_text('cat')
        return make_mock_process()

    with patch('shutil.which', return_value='/usr/bin/sex'), \
         patch('LBCgo.lbcregister.Popen', side_effect=fake_sex):
        go_sextractor(str(f), verbose=False)
    assert (d / 'test_1.cat').read_text() == 'cat'


def test_sextractor_accepts_source_extractor_name(tmp_path):
    """Debian/Ubuntu install the binary as 'source-extractor'."""
    only = lambda n: '/usr/bin/source-extractor' if n == 'source-extractor' else None
    cmd = _run_sextractor(tmp_path, which=only)
    assert cmd[0] == 'source-extractor'


def test_find_astromatic_tool_alternate_names():
    from LBCgo.lbcregister import find_astromatic_tool
    with patch('shutil.which',
               side_effect=lambda n: '/x/SWarp' if n == 'SWarp' else None):
        assert find_astromatic_tool('swarp') == 'SWarp'
    with patch('shutil.which', return_value=None):
        assert find_astromatic_tool('scamp') is None


def test_register_passes_sextractor_args(tmp_path):
    from LBCgo.lbcregister import go_register
    d = tmp_path / 'NGC891' / 'g'
    d.mkdir(parents=True)
    write_lbc_file(d, 'lbcb.20230101.000001_1.fits', imagetyp='object',
                   filter_name='g-SLOAN', object_name='NGC891',
                   nx=NX_SCIENCE)
    with patch('shutil.which', return_value='/usr/bin/sex'), \
         patch('LBCgo.lbcregister.Popen',
               return_value=make_mock_process()) as mp:
        go_register([str(d) + '/'], lbc_chips=[1], do_scamp=False,
                    do_swarp=False, verbose=False,
                    sextractor_args=dict(detect_thresh=3))
    assert _opt(mp.call_args[0][0], '-DETECT_THRESH') == '3'


# ---------------------------------------------------------------------------
# Joint SCAMP — Phase 0 §5.4
# ---------------------------------------------------------------------------

FIXTURES = Path(__file__).parent / 'fixtures' / 'scamp'


def _write_ldac(path, n_objects, chip_marker):
    """Minimal FITS_LDAC catalog: primary + (LDAC_IMHEAD, LDAC_OBJECTS).

    The IMHEAD is one row holding the header as 80-character cards padded
    with spaces, like SExtractor writes it (including FITSEXT/FITSNEXT).
    """
    from astropy.io import fits
    cards = [fits.Card('FITSFILE', 'chip_%d.fits' % chip_marker, 'File name').image,
             fits.Card('FITSEXT', 1, 'FITS Extension number').image,
             fits.Card('FITSNEXT', 1, 'Number of FITS image extensions').image,
             fits.Card('CHIPNO', chip_marker, 'test marker').image,
             'END'.ljust(80)]
    imhead = fits.BinTableHDU.from_columns(
        [fits.Column(name='Field Header Card', format='%dA' % (80 * len(cards)),
                     array=np.array([''.join(cards)]))], name='LDAC_IMHEAD')
    objs = fits.BinTableHDU.from_columns(
        [fits.Column(name='NUMBER', format='J', array=np.arange(n_objects))],
        name='LDAC_OBJECTS')
    fits.HDUList([fits.PrimaryHDU(), imhead, objs]).writeto(path, overwrite=True)


def _imhead_cards(hdu):
    """The 80-character cards of an LDAC_IMHEAD HDU, as raw bytes."""
    raw = hdu.data.tobytes()
    return [raw[i:i + 80] for i in range(0, len(raw), 80)]


def test_group_chips_by_exposure_orders_chips():
    from LBCgo.lbcregister import group_chips_by_exposure
    files = ['d/a_2.fits', 'd/b_1.fits', 'd/a_1.fits', 'd/b_2.fits']
    g = group_chips_by_exposure(files)
    assert g == {'d/a': ['d/a_1.fits', 'd/a_2.fits'],
                 'd/b': ['d/b_1.fits', 'd/b_2.fits']}


def test_group_chips_rejects_non_chip_name():
    from LBCgo.lbcregister import group_chips_by_exposure
    with pytest.raises(ValueError):
        group_chips_by_exposure(['d/not_a_chip.fits'])


def test_merge_ldac_concatenates_pairs_in_order(tmp_path):
    from astropy.io import fits
    from LBCgo.lbcregister import merge_ldac
    cats = []
    for chip, n in ((1, 3), (2, 5), (4, 2)):
        c = tmp_path / f'x_{chip}.cat'
        _write_ldac(c, n, chip)
        cats.append(str(c))
    out = merge_ldac(cats, str(tmp_path / 'x_exp.cat'))
    with fits.open(out) as h:
        assert [x.name for x in h] == ['PRIMARY'] + ['LDAC_IMHEAD', 'LDAC_OBJECTS'] * 3
        assert [len(h[i].data) for i in (2, 4, 6)] == [3, 5, 2]
        for ext, i in enumerate((1, 3, 5), start=1):
            cards = _imhead_cards(h[i])
            keyed = {c[:8].decode().strip(): fits.Card.fromstring(c.decode())
                     for c in cards if c[:8].strip() not in (b'END', b'')}
            assert keyed['CHIPNO'].value == (1, 2, 4)[ext - 1]
            # SCAMP needs FITSEXT/FITSNEXT to describe the merged file
            assert keyed['FITSEXT'].value == ext
            assert keyed['FITSNEXT'].value == 3
            assert keyed['FITSFILE'].value == 'x_exp.fits'
            # Every other byte is carried over untouched. (Re-writing the
            # table through astropy turns SExtractor's space padding into
            # NULs, which makes SCAMP fault.)
            with fits.open(cats[ext - 1]) as src:
                for new, old in zip(cards, _imhead_cards(src[1])):
                    if new[:8].strip() not in (b'FITSFILE', b'FITSEXT', b'FITSNEXT'):
                        assert new == old


def test_split_head_writes_one_file_per_section(tmp_path):
    from LBCgo.lbcregister import split_head
    head = tmp_path / 'e.head'
    head.write_text('CRPIX1  = 1\nEND\nCRPIX1  = 2\nEND\n')
    outs = [str(tmp_path / 'a_1.head'), str(tmp_path / 'a_2.head')]
    split_head(str(head), outs)
    assert Path(outs[0]).read_text() == 'CRPIX1  = 1\nEND\n'
    assert Path(outs[1]).read_text() == 'CRPIX1  = 2\nEND\n'


def test_split_head_section_count_mismatch(tmp_path):
    from LBCgo.lbcregister import split_head
    head = tmp_path / 'e.head'
    head.write_text('CRPIX1  = 1\nEND\n')
    with pytest.raises(ValueError, match='sections'):
        split_head(str(head), [str(tmp_path / 'a_1.head'),
                               str(tmp_path / 'a_2.head')])


def _qa_from_fixture(tmp_path, monkeypatch, **kw):
    """QA table from a real SCAMP 2.14.1 XML/head (3 exposures, 2 chips)."""
    import shutil
    from LBCgo.lbcregister import scamp_qa_table
    for n in range(3):
        shutil.copy(FIXTURES / 'lbcb.20140101.000000_exp.head',
                    tmp_path / f'lbcb.20140101.00000{n}_exp.head')
    shutil.copy(FIXTURES / 'scamp.xml', tmp_path / 'scamp.xml')
    monkeypatch.chdir(tmp_path)
    groups = {f'lbcb.20140101.00000{n}': [f'lbcb.20140101.00000{n}_1.fits',
                                          f'lbcb.20140101.00000{n}_2.fits']
              for n in range(3)}
    return scamp_qa_table('scamp.xml', groups, **kw)


def test_qa_table_values_from_real_scamp_output(tmp_path, monkeypatch):
    qa = _qa_from_fixture(tmp_path, monkeypatch)
    assert len(qa) == 6
    assert list(qa['chip'][:2]) == [1, 2]
    # ASTRRMS is in degrees in the .head; the table is in arcsec
    assert qa['ref_rms_x'][0] == pytest.approx(1.204308070295e-05 * 3600)
    assert qa['xy_contrast'][0] == pytest.approx(4.794, abs=1e-3)
    assert qa.meta['scamp_version'] == '2.14.1'
    assert qa.meta['astref_catalog'] == 'GAIA-DR3'
    assert not any(qa['bad'])


def test_qa_flags_by_threshold(tmp_path, monkeypatch):
    qa = _qa_from_fixture(tmp_path, monkeypatch, max_ref_rms=0.01,
                          min_xy_contrast=10.0)
    assert all(qa['bad'])
    assert all('ref_rms' in r and 'low_contrast' in r for r in qa['reason'])


def test_qa_flags_missing_reference_match(tmp_path, monkeypatch):
    import re
    _qa_from_fixture(tmp_path, monkeypatch)
    head = tmp_path / 'lbcb.20140101.000000_exp.head'
    head.write_text(re.sub(r'(ASTRRMS[12]=)\s*\S+', r'\1 0.0', head.read_text()))
    from LBCgo.lbcregister import scamp_qa_table
    groups = {'lbcb.20140101.000000': ['lbcb.20140101.000000_1.fits',
                                       'lbcb.20140101.000000_2.fits']}
    qa = scamp_qa_table('scamp.xml', groups)
    assert all(qa['reason'] == 'no_ref_match')


def _joint_setup(tmp_path, n_exp=2, chips=(1, 2)):
    from conftest import write_lbc_file
    files = []
    for k in range(n_exp):
        for chip in chips:
            name = f'lbcb.20230101.00000{k}_{chip}.fits'
            write_lbc_file(tmp_path, name, imagetyp='object',
                           filter_name='g-SLOAN', object_name='X',
                           nx=NX_SCIENCE)
            _write_ldac(tmp_path / name.replace('.fits', '.cat'), 4, chip)
            files.append(str(tmp_path / name))
    return files


def test_joint_scamp_single_run_per_iteration(tmp_path):
    """One SCAMP call per iteration over all exposure catalogs, with the
    focal-plane / instrument options set and the chip heads split out."""
    from LBCgo.lbcregister import go_scamp_joint
    files = _joint_setup(tmp_path, n_exp=2)

    def fake_scamp(cmd, **kw):
        for cat in [c for c in cmd if c.endswith('_exp.cat')]:
            head = cat.replace('.cat', '.head')
            (Path(kw['cwd']) / head).write_text('CRPIX1  = 1\nEND\nCRPIX1  = 2\nEND\n')
        return make_mock_process()

    with patch('shutil.which', return_value='/usr/bin/scamp'), \
         patch('LBCgo.lbcregister.Popen', side_effect=fake_scamp) as mp, \
         patch('LBCgo.lbcregister.scamp_qa_table',
               return_value=Table({'bad': [False]})):
        go_scamp_joint(files, num_iterations=3, verbose=False, qa_file=None)
    assert mp.call_count == 3
    cmds = [c[0][0] for c in mp.call_args_list]
    assert [_opt(c, '-MOSAIC_TYPE') for c in cmds] == \
        ['LOOSE', 'FIX_FOCALPLANE', 'FIX_FOCALPLANE']
    for c in cmds:
        assert sum(a.endswith('_exp.cat') for a in c) == 2
        assert _opt(c, '-STABILITY_TYPE') == 'INSTRUMENT'
        assert _opt(c, '-ASTRINSTRU_KEY') == 'FILTER'
        assert _opt(c, '-ASTREFEPOCH_TYPE') == 'FIELDS_AVERAGE'
    for f in files:
        assert Path(f.replace('.fits', '.head')).exists()


def test_joint_scamp_raises_on_scamp_failure(tmp_path):
    from LBCgo.lbcregister import go_scamp_joint
    files = _joint_setup(tmp_path, n_exp=1)
    crashed = MagicMock()
    crashed.wait.return_value = -10
    crashed.returncode = -10
    with patch('shutil.which', return_value='/usr/bin/scamp'), \
         patch('LBCgo.lbcregister.Popen', return_value=crashed):
        with pytest.raises(RuntimeError, match='status -10'):
            go_scamp_joint(files, verbose=False)


def test_joint_scamp_requires_catalogs(tmp_path):
    from LBCgo.lbcregister import go_scamp_joint
    with patch('shutil.which', return_value='/usr/bin/scamp'):
        with pytest.raises(FileNotFoundError):
            go_scamp_joint([str(tmp_path / 'a_1.fits')], verbose=False)


def test_register_uses_joint_scamp_by_default(tmp_path):
    from LBCgo.lbcregister import go_register
    files = _joint_setup(tmp_path, n_exp=2)
    with patch('LBCgo.lbcregister.go_scamp_joint') as joint, \
         patch('LBCgo.lbcregister.go_scamp') as legacy:
        go_register([str(tmp_path) + '/'], lbc_chips=[1, 2],
                    do_sextractor=False, do_swarp=False, verbose=False)
    assert joint.call_count == 1 and legacy.call_count == 0
    assert len(joint.call_args[0][0]) == 4


def test_register_legacy_per_chip_scamp(tmp_path):
    from LBCgo.lbcregister import go_register
    _joint_setup(tmp_path, n_exp=2)
    with patch('LBCgo.lbcregister.go_scamp_joint') as joint, \
         patch('LBCgo.lbcregister.go_scamp') as legacy:
        go_register([str(tmp_path) + '/'], lbc_chips=[1, 2],
                    do_sextractor=False, do_swarp=False,
                    scamp_joint=False, verbose=False)
    assert joint.call_count == 0 and legacy.call_count == 4


def test_scamp_legacy_passes_mosaic_type(tmp_path):
    """The per-iteration MOSAIC_TYPE is no longer dead code."""
    from LBCgo.lbcregister import go_scamp
    with patch('shutil.which', return_value='/usr/bin/scamp'), \
         patch('LBCgo.lbcregister.Popen',
               return_value=make_mock_process()) as mp, \
         patch('LBCgo.lbcregister.votable.parse',
               return_value=make_mock_votable()):
        go_scamp(str(tmp_path / 'test_1.fits'), num_iterations=3, verbose=False)
    mos = [_opt(c[0][0], '-MOSAIC_TYPE') for c in mp.call_args_list]
    assert mos == ['LOOSE', 'FIX_FOCALPLANE', 'FIX_FOCALPLANE']


# ---------------------------------------------------------------------------
# go_swarp — §5.5 options (weights, flux scale, combine type, background)
# ---------------------------------------------------------------------------

def _swarp_cmd(tmp_path, names=('img1_1.fits',), sidecars=(), heads=(),
               **kwargs):
    """Run go_swarp with a mocked Popen and return (cmd list, Popen kwargs)."""
    from LBCgo.lbcregister import go_swarp
    files = [str(write_lbc_file(tmp_path, n, imagetyp='object',
                                filter_name='g-SLOAN', nx=NX_SCIENCE))
             for n in names]
    for n in sidecars:
        (tmp_path / n).write_bytes(b'')
    for n, txt in heads:
        (tmp_path / n).write_text(txt)
    # The packaged config lives under a path with a space on this machine,
    # which would trigger staging; use a space-free copy unless overridden.
    if 'configfile' not in kwargs:
        conf = tmp_path / 'swarp.conf'
        conf.write_text('')
        kwargs['configfile'] = str(conf)
    with patch('shutil.which', return_value='/usr/bin/swarp'), \
         patch('LBCgo.lbcregister.Popen',
               return_value=make_mock_process()) as mp, \
         patch('astropy.io.fits.setval'):
        go_swarp(files, verbose=False, **kwargs)
    return mp.call_args[0][0], mp.call_args[1]


def _opt(cmd, flag):
    return cmd[cmd.index(flag) + 1]


def test_swarp_defaults_clipped_background_1024(tmp_path):
    cmd, _ = _swarp_cmd(tmp_path)
    assert _opt(cmd, '-COMBINE_TYPE') == 'CLIPPED'
    assert _opt(cmd, '-SUBTRACT_BACK') == 'Y'
    assert _opt(cmd, '-BACK_SIZE') == '1024'
    # No weight sidecar and no head -> unweighted, unscaled
    assert _opt(cmd, '-WEIGHT_TYPE') == 'NONE'
    assert _opt(cmd, '-FSCALE_KEYWORD') == 'NONE'


def test_swarp_uses_weight_maps_when_all_present(tmp_path):
    cmd, _ = _swarp_cmd(tmp_path, names=('a_1.fits', 'a_2.fits'),
                        sidecars=('a_1.weight.fits', 'a_2.weight.fits'))
    assert _opt(cmd, '-WEIGHT_TYPE') == 'MAP_WEIGHT'
    assert _opt(cmd, '-WEIGHT_SUFFIX') == '.weight.fits'


def test_swarp_unweighted_if_any_weight_missing(tmp_path):
    cmd, _ = _swarp_cmd(tmp_path, names=('a_1.fits', 'a_2.fits'),
                        sidecars=('a_1.weight.fits',))
    assert _opt(cmd, '-WEIGHT_TYPE') == 'NONE'


def test_swarp_fscale_from_head(tmp_path):
    head = 'FLXSCALE=              1.0231 / relative flux scale\nEND\n'
    cmd, _ = _swarp_cmd(tmp_path, heads=[('img1_1.head', head)])
    assert _opt(cmd, '-FSCALE_KEYWORD') == 'FLXSCALE'
    cmd, _ = _swarp_cmd(tmp_path, heads=[('img1_1.head', head)],
                        use_fscale=False)
    assert _opt(cmd, '-FSCALE_KEYWORD') == 'NONE'


def test_swarp_options_forwarded(tmp_path):
    cmd, _ = _swarp_cmd(tmp_path, combine_type='median',
                        subtract_back=False, clip_sigma=3.0)
    assert _opt(cmd, '-COMBINE_TYPE') == 'MEDIAN'
    assert _opt(cmd, '-SUBTRACT_BACK') == 'N'
    assert '-BACK_SIZE' not in cmd
    assert _opt(cmd, '-CLIP_SIGMA') == '3.0'


def test_swarp_rejects_bad_combine_type(tmp_path):
    from LBCgo.lbcregister import go_swarp
    with pytest.raises(ValueError, match="combine_type"):
        go_swarp(['x_1.fits'], combine_type='SUM')


def test_swarp_nonzero_exit_raises(tmp_path):
    from LBCgo.lbcregister import go_swarp
    f = write_lbc_file(tmp_path, 'img1_1.fits', imagetyp='object',
                       filter_name='g-SLOAN', nx=NX_SCIENCE)
    bad = make_mock_process()
    bad.wait.return_value = 1
    with patch('shutil.which', return_value='/usr/bin/swarp'), \
         patch('LBCgo.lbcregister.Popen', return_value=bad):
        with pytest.raises(RuntimeError, match="SWarp"):
            go_swarp([str(f)], verbose=False)


def test_swarp_stages_paths_with_spaces(tmp_path):
    spaced = tmp_path / 'with space'
    spaced.mkdir()
    conf = spaced / 'swarp.conf'
    spaced.mkdir(exist_ok=True)
    conf.write_text('')
    from LBCgo.lbcregister import go_swarp
    f = write_lbc_file(spaced, 'img1_1.fits', imagetyp='object',
                       filter_name='g-SLOAN', nx=NX_SCIENCE)
    out = {}

    def fake_popen(cmd, **kw):
        # Pretend SWarp wrote its products where it was told to
        for flag in ('-IMAGEOUT_NAME', '-WEIGHTOUT_NAME'):
            open(cmd[cmd.index(flag) + 1], 'wb').close()
        out['cmd'], out['kw'] = cmd, kw
        return make_mock_process()

    with patch('shutil.which', return_value='/usr/bin/swarp'), \
         patch('LBCgo.lbcregister.Popen', side_effect=fake_popen), \
         patch('astropy.io.fits.setval'):
        go_swarp([str(f)], configfile=str(conf), verbose=False)
    assert out['kw']['cwd'] is not None
    assert not any(' ' in c for c in out['cmd'])
    assert (tmp_path / 'NGC891.g.mos.fits').exists()
    assert (tmp_path / 'NGC891.g.mos.weight.fits').exists()


def test_register_forwards_swarp_args(tmp_path):
    from LBCgo.lbcregister import go_register
    with patch('LBCgo.lbcregister.go_swarp') as gs:
        write_lbc_file(tmp_path, 'o_1.fits', imagetyp='object',
                       filter_name='g-SLOAN', nx=NX_SCIENCE)
        go_register(str(tmp_path), lbc_chips=[1], do_sextractor=False,
                    do_scamp=False, swarp_args=dict(combine_type='MEDIAN'),
                    verbose=False)
    assert gs.call_args[1]['combine_type'] == 'MEDIAN'
