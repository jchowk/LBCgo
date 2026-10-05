"""
Smoke test for calibration/gain_rdnoise_example/run.py: run it on synthetic
raw frames with known gain and read noise and check what it writes.

The script writes its outputs next to itself, so it is copied into tmp_path
together with a synthetic inputs.ecsv; nothing is written into the repository.
"""

import json
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
EXAMPLE = REPO / 'calibration' / 'gain_rdnoise_example'

TRUTH = {'LBCB': ([1.96, 2.09, 2.06, 1.98], [11.4, 11.6, 11.6, 11.2]),
         'LBCR': ([2.08, 2.14, 2.13, 2.09], [9.8, 9.5, 9.9, 9.6])}
INSTRUME = {'LBCB': 'LBC_BLUE', 'LBCR': 'LBC-RED '}


def _write_raw(path, frames, instrume, mjd, imagetyp):
    """Raw-format MEF: 20 px prescan + science + 20 px overscan per chip."""
    primary = fits.Header({'INSTRUME': instrume, 'MJD_OBS': mjd,
                           'IMAGETYP': imagetyp})
    hdul = fits.HDUList([fits.PrimaryHDU(header=primary)])
    for chip, sci in enumerate(frames, start=1):
        ny, nx = sci.shape
        full = np.empty((ny, nx + 40), dtype='float32')
        full[:, :20] = 1000.0
        full[:, 20:20 + nx] = sci
        full[:, 20 + nx:] = 1000.0
        hdul.append(fits.ImageHDU(full, header=fits.Header({
            'EXTNAME': f'LBCCHIP{chip}', 'BUNIT': 'adu', 'SATURATE': 65536,
            'TRIMSEC': f'[21:{20 + nx},1:{ny}]',
            'BIASSEC': f'[{21 + nx}:{40 + nx},1:{ny}]'})))
    hdul.writeto(path)


@pytest.fixture
def example_run(tmp_path):
    """Copy run.py, make synthetic raw frames + inputs.ecsv, run it."""
    rng = np.random.default_rng(11)
    shape = (300, 300)
    workdir = tmp_path / 'gain_rdnoise_test'
    rawdir = tmp_path / 'raw'
    workdir.mkdir()
    rawdir.mkdir()
    shutil.copy(EXAMPLE / 'run.py', workdir / 'run.py')

    rows = []
    for set_id, (channel, mjd) in enumerate([('LBCB', 56981.3),
                                             ('LBCR', 57020.6)], start=1):
        gains, rns = TRUTH[channel]
        prnu = [1 + 0.01 * rng.standard_normal(shape) for _ in gains]
        for role, signals in (('flat', (30000.0, 31500.0)),
                              ('bias', (0.0, 0.0))):
            for k, signal in enumerate(signals):
                frames = []
                for g, rn, p in zip(gains, rns, prnu):
                    e = rng.poisson(signal * p) if signal > 0 else 0.0
                    frames.append(1000.0 + (e + rng.normal(0, rn, shape)) / g)
                name = f'{channel.lower()}.set{set_id}.{role}{k}.fits'
                _write_raw(rawdir / name, frames, INSTRUME[channel], mjd,
                           'flat' if role == 'flat' else 'zero')
                rows.append((name, f'{channel}{set_id}{role}{k}', mjd,
                             channel, 'r-SLOAN', 5.0, role, set_id, ''))
    Table(rows=rows, names=['filename', 'obs_id', 'mjd_obs', 'channel',
                            'filter', 'exptime', 'role', 'set', 'notes']
          ).write(workdir / 'inputs.ecsv', format='ascii.ecsv')

    env = dict(os.environ, LBCGO_RAW=str(rawdir),
               PYTHONPATH=os.pathsep.join(
                   [str(REPO), os.environ.get('PYTHONPATH', '')]))
    proc = subprocess.run([sys.executable, str(workdir / 'run.py')],
                          env=env, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    return workdir


def test_example_recovers_gain_and_read_noise(example_run):
    rows = Table.read(example_run / 'detector_rows.ecsv', format='ascii.ecsv')
    assert len(rows) == 8
    for row in rows:
        gains, rns = TRUTH[row['channel']]
        assert abs(row['gain'] / gains[row['chip'] - 1] - 1) < 0.03
        assert abs(row['rdnoise'] / rns[row['chip'] - 1] - 1) < 0.05
        assert row['source'].startswith('calibration/gain_rdnoise_test')


def test_example_writes_provenance(example_run):
    log = json.loads((example_run / 'run_log.json').read_text())
    assert set(log) >= {'date_utc', 'lbcgo_version', 'lbcgo_commit',
                        'params', 'inputs_sha256', 'outputs_sha256'}
    assert len(log['inputs_sha256']) == 8
    per_set = Table.read(example_run / 'results_per_set.ecsv',
                         format='ascii.ecsv')
    assert len(per_set) == 8 and set(per_set['set']) == {1, 2}


def test_example_requires_data_path(tmp_path):
    shutil.copy(EXAMPLE / 'run.py', tmp_path / 'run.py')
    env = {k: v for k, v in os.environ.items() if k != 'LBCGO_RAW'}
    proc = subprocess.run([sys.executable, str(tmp_path / 'run.py')],
                          env=env, capture_output=True, text=True)
    assert proc.returncode != 0
    assert 'LBCGO_RAW' in proc.stderr


def test_committed_example_inputs_have_required_columns():
    t = Table.read(EXAMPLE / 'inputs.ecsv', format='ascii.ecsv')
    assert set(t.colnames) >= {'filename', 'obs_id', 'mjd_obs', 'channel',
                               'filter', 'exptime', 'role', 'set', 'notes'}
    for set_id in set(t['set']):
        roles = list(t['role'][t['set'] == set_id])
        assert sorted(roles) == ['bias', 'bias', 'flat', 'flat']
