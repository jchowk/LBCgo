"""
Tests for calibration/make_inputs.py: building a calibration inputs.ecsv from
a directory of raw LBC frames with ccdproc's ImageFileCollection, and feeding
it to calibration/gain_rdnoise_example/run.py.
"""

import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

REPO = Path(__file__).resolve().parents[1]
CALIB = REPO / 'calibration'

spec = importlib.util.spec_from_file_location('make_inputs',
                                              CALIB / 'make_inputs.py')
make_inputs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(make_inputs)

SHAPE = (200, 200)
INSTRUME = {'lbcb': 'LBC_BLUE', 'lbcr': 'LBC-RED '}


def write_raw(directory, name, imagetyp, obj, filt, exptime, mjd,
              frames=None, lbcobnam='ob', level=1000.0):
    """Raw-format LBC MEF (20 px prescan, 20 px overscan per chip)."""
    channel = name[:4]
    primary = fits.Header({
        'IMAGETYP': imagetyp, 'OBJECT': obj, 'FILTER': filt,
        'EXPTIME': exptime, 'MJD_OBS': mjd, 'OBS_ID': name[:-5],
        'LBCOBNAM': lbcobnam, 'PROPID': 'TEST', 'AIRMASS': 1.1,
        'INSTRUME': INSTRUME[channel]})
    hdul = fits.HDUList([fits.PrimaryHDU(header=primary)])
    for chip in range(4):
        sci = (frames[chip] if frames is not None
               else np.full(SHAPE, level, dtype='f4'))
        ny, nx = sci.shape
        full = np.empty((ny, nx + 40), dtype='float32')
        full[:, :20] = 1000.0
        full[:, 20:20 + nx] = sci
        full[:, 20 + nx:] = 1000.0
        hdul.append(fits.ImageHDU(full, header=fits.Header({
            'EXTNAME': f'LBCCHIP{chip + 1}', 'BUNIT': 'adu',
            'SATURATE': 65536, 'TRIMSEC': f'[21:{20 + nx},1:{ny}]',
            'BIASSEC': f'[{21 + nx}:{40 + nx},1:{ny}]'})))
    hdul.writeto(Path(directory) / name)


@pytest.fixture
def raw_dir(tmp_path):
    d = tmp_path / 'raw'
    d.mkdir()
    # LBCB night 56981: 2 sky flats, 2 biases, 1 science frame, 1 test flat
    write_raw(d, 'lbcb.20141120.010000.fits', 'flat', 'SkyFlat', 'g-SLOAN',
              5.0, 56981.04, level=21000.0)
    write_raw(d, 'lbcb.20141120.010100.fits', 'flat', 'SkyFlat', 'g-SLOAN',
              5.0, 56981.05, level=22000.0)
    write_raw(d, 'lbcb.20141120.120000.fits', 'zero', 'Bias', 'g-SLOAN',
              0.0, 56981.50)
    write_raw(d, 'lbcb.20141120.120100.fits', 'zero', 'Bias', 'g-SLOAN',
              0.0, 56981.51)
    write_raw(d, 'lbcb.20141120.065509.fits', 'object', 'NGC 891',
              'SDT_Uspec', 200.0, 56981.29)
    write_raw(d, 'lbcb.20141120.010200.fits', 'flat', 'SkyFlat', 'g-SLOAN',
              5.0, 56981.06, lbcobnam='SkyFlatTest_x', level=20000.0)
    # LBCR night 57020: 2 flats, 2 biases
    for t, typ, mjd in (('132703', 'flat', 57020.56), ('132733', 'flat', 57020.57),
                        ('140000', 'zero', 57020.58), ('140030', 'zero', 57020.59)):
        write_raw(d, f'lbcr.20141229.{t}.fits', typ,
                  'SkyFlat' if typ == 'flat' else 'Bias', 'I-BESSEL',
                  20.0 if typ == 'flat' else 0.0, mjd, level=25000.0
                  if typ == 'flat' else 1000.0)
    # Not an LBC raw frame: must be ignored
    fits.PrimaryHDU(header=fits.Header({'IMAGETYP': 'flat'})).writeto(
        d / 'other.fits')
    return d


def test_build_inputs_default(raw_dir):
    t = make_inputs.build_inputs(raw_dir)
    assert len(t) == 9                     # SkyFlatTest and other.fits dropped
    assert set(t.colnames) >= {'filename', 'obs_id', 'mjd_obs', 'channel',
                               'filter', 'exptime', 'role', 'set', 'notes'}
    assert 'lbcb.20141120.010200.fits' not in t['filename']
    assert 'other.fits' not in t['filename']
    roles = dict(zip(t['filename'], t['role']))
    assert roles['lbcb.20141120.010000.fits'] == 'flat'
    assert roles['lbcb.20141120.120000.fits'] == 'bias'
    assert roles['lbcb.20141120.065509.fits'] == 'science'
    # one set per channel + UT date
    assert set(t['set'][t['channel'] == 'LBCB']) == {1}
    assert set(t['set'][t['channel'] == 'LBCR']) == {2}
    assert all(n == '' for n in t['notes'])
    row = t[t['filename'] == 'lbcb.20141120.065509.fits'][0]
    assert (row['object'], row['filter'], row['exptime']) == \
        ('NGC 891', 'SDT_Uspec', 200.0)
    assert row['obs_id'] == 'lbcb.20141120.065509'


def test_build_inputs_filters(raw_dir):
    t = make_inputs.build_inputs(raw_dir, channels=['LBCR'],
                                 imagetyp=['flat', 'zero'])
    assert len(t) == 4 and set(t['channel']) == {'LBCR'}
    t = make_inputs.build_inputs(raw_dir, imagetyp=['object'],
                                 objects=['ngc 891'])          # case-insensitive
    assert list(t['filename']) == ['lbcb.20141120.065509.fits']
    t = make_inputs.build_inputs(raw_dir, filters=['I-BESSEL'])
    assert set(t['channel']) == {'LBCR'}
    assert len(make_inputs.build_inputs(raw_dir, keep_tests=True)) == 10


def test_build_inputs_levels(raw_dir):
    t = make_inputs.build_inputs(raw_dir, channels=['LBCB'], levels=True)
    lev = dict(zip(t['filename'], t['median_adu']))
    assert abs(lev['lbcb.20141120.010000.fits'] - 20000.0) < 1    # 21000-1000
    assert abs(lev['lbcb.20141120.120000.fits']) < 1
    assert np.isnan(lev['lbcb.20141120.065509.fits'])


def test_main_writes_and_protects_output(raw_dir, tmp_path, capsys):
    out = tmp_path / 'inputs.ecsv'
    assert make_inputs.main([str(raw_dir), '-o', str(out)]) == 0
    assert len(Table.read(out, format='ascii.ecsv')) == 9
    assert 'Wrote 9 frames' in capsys.readouterr().out
    with pytest.raises(SystemExit):
        make_inputs.main([str(raw_dir), '-o', str(out)])
    assert make_inputs.main([str(raw_dir), '-o', str(out), '--overwrite',
                             '--imagetyp', 'object']) == 0
    assert len(Table.read(out, format='ascii.ecsv')) == 1
    assert make_inputs.main([str(raw_dir), '-o', str(tmp_path / 'x.ecsv'),
                             '--object', 'nothing']) == 1


def test_inputs_feed_gain_example(tmp_path):
    """make_inputs -> run.py end to end, with known gain/read noise."""
    rng = np.random.default_rng(5)
    truth = {'lbcb': ([1.96, 2.09, 2.06, 1.98], [11.4, 11.6, 11.6, 11.2]),
             'lbcr': ([2.08, 2.14, 2.13, 2.09], [9.8, 9.5, 9.9, 9.6])}
    raw = tmp_path / 'raw'
    raw.mkdir()
    for ch, mjd, date in (('lbcb', 56981.0, '20141120'),
                          ('lbcr', 57020.0, '20141229')):
        gains, rns = truth[ch]
        prnu = [1 + 0.01 * rng.standard_normal(SHAPE) for _ in gains]
        for k, (typ, signal) in enumerate((('flat', 30000.), ('flat', 31500.),
                                           ('zero', 0.), ('zero', 0.))):
            frames = [1000.0 + ((rng.poisson(signal * p) if signal else 0.0)
                                + rng.normal(0, rn, SHAPE)) / g
                      for g, rn, p in zip(gains, rns, prnu)]
            write_raw(raw, f'{ch}.{date}.00000{k}.fits', typ,
                      'SkyFlat' if typ == 'flat' else 'Bias', 'r-SLOAN',
                      5.0 if signal else 0.0, mjd + 0.01 * k, frames=frames)

    product = tmp_path / 'gain_rdnoise_test'
    product.mkdir()
    shutil.copy(CALIB / 'gain_rdnoise_example' / 'run.py', product / 'run.py')
    assert make_inputs.main([str(raw), '-o', str(product / 'inputs.ecsv'),
                             '--imagetyp', 'flat', 'zero']) == 0

    env = dict(os.environ, LBCGO_RAW=str(raw),
               PYTHONPATH=os.pathsep.join(
                   [str(REPO), os.environ.get('PYTHONPATH', '')]))
    proc = subprocess.run([sys.executable, str(product / 'run.py')],
                          env=env, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    rows = Table.read(product / 'detector_rows.ecsv', format='ascii.ecsv')
    assert len(rows) == 8
    for row in rows:
        gains, rns = truth[row['channel'].lower()]
        assert abs(row['gain'] / gains[row['chip'] - 1] - 1) < 0.03
        assert abs(row['rdnoise'] / rns[row['chip'] - 1] - 1) < 0.05


def test_sets_split_by_ut_date(tmp_path):
    """Same channel on two UT dates -> two sets; numbering follows time."""
    for name, mjd in (('lbcb.20141120.010000.fits', 56981.04),
                      ('lbcb.20141120.120000.fits', 56981.50),
                      ('lbcb.20141121.010000.fits', 56982.04),
                      ('lbcr.20141120.010000.fits', 56981.04)):
        write_raw(tmp_path, name, 'flat', 'SkyFlat', 'g-SLOAN', 5.0, mjd)
    t = make_inputs.build_inputs(tmp_path)
    sets = dict(zip(t['filename'], t['set']))
    assert sets['lbcb.20141120.010000.fits'] == sets['lbcb.20141120.120000.fits']
    assert sets['lbcb.20141121.010000.fits'] != sets['lbcb.20141120.010000.fits']
    assert sets['lbcr.20141120.010000.fits'] not in (
        sets['lbcb.20141120.010000.fits'], sets['lbcb.20141121.010000.fits'])
    assert sorted(set(t['set'])) == [1, 2, 3]
