"""Measure per-chip gain and read noise (photon-transfer) for one product.

Product: candidate rows for LBCgo/conf/lbc_detector.ecsv.
Inputs:  inputs.ecsv in this directory. One row per raw frame, with a `set`
         column grouping two flats + two biases of one channel/epoch.
Usage:   LBCGO_RAW=/path/to/raw python run.py [--raw DIR] [--box 1000]
         (--box 0: whole trimmed chip)

This script only orchestrates: the measurement is
LBCgo.detector.measure_gain_rdnoise_files and the combination of sets is
LBCgo.detector.summarize_gain_rdnoise (both tested in tests/test_detector.py).
The apparent photon-transfer gain rises with flat level (brighter-fatter
effect), so with PARAMS['gain_model'] = 'linear' the product gain is the
zero-level intercept of a straight-line fit of gain against level per chip,
and the read noise is that gain times the median read noise in ADU. That
per-pixel gain sets the per-pixel variance (weight maps). The flux gain
gain_flux (electrons per ADU of a flux summed over pixels) is the median
gain_sum of the sets whose two flats are at most PARAMS['flux_max_dt']
seconds apart; it differs from the per-pixel gain where the readout
correlates neighbouring pixels.
It writes, next to itself:
  results_per_set.ecsv   every measurement (one row per set and chip, with
                         level, read noise in ADU and correlation diagnostics)
  gain_fit.ecsv          the fit per channel/chip: intercept, slope, scatter,
                         and the brighter-fatter diagnostics. The slope of
                         gain_sum (covariances summed to max_lag added back)
                         separates brighter-fatter (slope ~0) from
                         non-linearity (slope ~ that of gain)
  detector_rows.ecsv     one row per channel/chip (gain, rdnoise, gain_flux),
                         ready to be reviewed and merged into
                         conf/lbc_detector.ecsv
  run_log.json           LBCgo version/commit, parameters, input checksums
Installing the rows into the package table is a separate, reviewed step.
"""

import argparse
import datetime
import hashlib
import json
import os
import subprocess
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack

import LBCgo
from LBCgo import detector

HERE = Path(__file__).resolve().parent

# Parameters that define this product (change -> new product directory)
PARAMS = {
    'box': 1000,          # central region of each trimmed chip [px];
                          # None = whole chip (less noise in gain_sum)
    'sigma': 4.0,         # clipping threshold
    'cell': 50,           # block size for the variances [px]
    'max_lag': 3,         # lags summed for gain_sum (0 = skip)
    'flux_max_dt': 60.0,  # gain_flux: only pairs with flats <= this many
                          # seconds apart (None = all pairs)
    'gain_model': 'linear',  # 'linear': intercept of gain vs level;
                             # 'median': median over sets
    'min_sets': 3,        # fewer sets per chip -> median instead of the fit
    'mjd_start': 55273.0,  # validity range written to the product rows
                           # (2010-03-18 UT, first night of the data)
    'mjd_end': np.nan,
}


def sha256(path, blocksize=2**20):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(blocksize), b''):
            h.update(block)
    return h.hexdigest()


def git_commit():
    """Commit of the LBCgo code actually imported, '+dirty' if modified."""
    pkg = Path(LBCgo.__file__).resolve().parent
    try:
        run = lambda *a: subprocess.check_output(
            ['git', *a], cwd=pkg, text=True, stderr=subprocess.DEVNULL).strip()
        return run('rev-parse', 'HEAD') + ('+dirty' if run('status',
                                           '--porcelain', '.') else '')
    except (OSError, subprocess.CalledProcessError):
        return 'unknown (LBCgo not run from a git checkout)'


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--raw', default=os.environ.get('LBCGO_RAW'),
                    help='directory holding the raw frames (or $LBCGO_RAW)')
    ap.add_argument('--box', type=int, default=PARAMS['box'],
                    help='central region [px]; 0 = whole chip')
    args = ap.parse_args()
    if not args.raw:
        ap.error('give --raw or set LBCGO_RAW')
    raw = Path(args.raw)
    params = dict(PARAMS, box=args.box if args.box and args.box > 0
                  else None)

    inputs = Table.read(HERE / 'inputs.ecsv', format='ascii.ecsv')
    results = []
    for set_id in sorted(set(inputs['set'])):
        rows = inputs[inputs['set'] == set_id]
        flats = [raw / f for f in rows[rows['role'] == 'flat']['filename']]
        biases = [raw / f for f in rows[rows['role'] == 'bias']['filename']]
        if len(flats) != 2 or len(biases) != 2:
            raise ValueError(f'set {set_id}: need 2 flats and 2 biases')

        # Sanity check: inputs.ecsv must agree with the headers
        channels = {detector.lbc_channel(fits.getheader(p), str(p))
                    for p in flats + biases}
        if channels != set(rows['channel']):
            raise ValueError(f'set {set_id}: channel mismatch {channels}')

        t = detector.measure_gain_rdnoise_files(
            *map(str, flats + biases), box=params['box'],
            sigma=params['sigma'], cell=params['cell'],
            max_lag=params['max_lag'])
        t['set'] = set_id
        t['mjd'] = float(np.mean(rows['mjd_obs']))
        flat_mjd = rows['mjd_obs'][rows['role'] == 'flat']
        t['flat_dt'] = float(abs(flat_mjd[1] - flat_mjd[0]) * 86400.0)
        t['flat_dt'].unit = 's'
        if 'lbcobnam' in rows.colnames:
            obs = set(rows['lbcobnam'][rows['role'] == 'flat'])
            if len(obs) > 1:
                print(f'WARNING set {set_id}: flats from different OBs '
                      f'{sorted(obs)} (rotator angle?)')
        results.append(t)
        print(f'set {set_id}:', ', '.join(
            f"chip {r['chip']} level={r['level']:.0f} g={r['gain']:.3f} "
            f"RN={r['rdnoise']:.2f}"
            for r in t))

    results = vstack(results)
    results.write(HERE / 'results_per_set.ecsv', format='ascii.ecsv',
                  overwrite=True)

    # One product row per channel/chip (see the module docstring)
    commit = git_commit()
    source = f'calibration/{HERE.name} (LBCgo {commit[:10]})'
    product, fit = detector.summarize_gain_rdnoise(
        results, model=params['gain_model'], min_sets=params['min_sets'],
        source=source, mjd_start=params['mjd_start'],
        mjd_end=params['mjd_end'], flux_max_dt=params['flux_max_dt'])
    fit.write(HERE / 'gain_fit.ecsv', format='ascii.ecsv', overwrite=True)
    for r in fit:
        print(f"{r['channel']} chip {r['chip']}: {r['model']} n={r['n']} "
              f"g0={r['gain0']:.3f}+-{r['gain0_err']:.3f} "
              f"slope={r['slope_pct_per_10k']:+.2f}%/10k ADU "
              f"(gain_nn {r['gain_nn_slope_pct_per_10k']:+.2f}, gain_sum "
              f"{r['gain_sum_slope_pct_per_10k']:+.2f}"
              f"+-{r['gain_sum_slope_pct_err']:.2f} %/10k) "
              f"RN={r['rdnoise']:.2f} e- gain_flux={r['gain_flux']:.3f}"
              f"+-{r['gain_flux_err']:.3f} (n={r['n_flux']})")
    detector.write_detector_table(product, str(HERE / 'detector_rows.ecsv'),
                                  overwrite=True)

    log = {
        'date_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'lbcgo_version': LBCgo.__version__,
        'lbcgo_commit': commit,
        'params': {k: (None if isinstance(v, float) and np.isnan(v) else v)
                   for k, v in params.items()},
        'inputs_sha256': {str(f): sha256(raw / f)
                          for f in inputs['filename']},
        'outputs_sha256': {name: sha256(HERE / name) for name in
                           ('results_per_set.ecsv', 'gain_fit.ecsv',
                            'detector_rows.ecsv')},
    }
    (HERE / 'run_log.json').write_text(json.dumps(log, indent=2) + '\n')
    print(f'wrote {len(product)} rows to detector_rows.ecsv')


if __name__ == '__main__':
    main()
