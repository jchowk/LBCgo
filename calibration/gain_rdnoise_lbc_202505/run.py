"""Measure per-chip gain and read noise (photon-transfer) for one product.

Product: candidate rows for LBCgo/conf/lbc_detector.ecsv.
Inputs:  inputs.ecsv in this directory. One row per raw frame, with a `set`
         column grouping two flats + two biases of one channel/epoch.
Usage:   LBCGO_RAW=/path/to/raw python run.py [--raw DIR] [--box 1000]

Writes next to itself: results_per_set.ecsv (every measurement),
detector_rows.ecsv (median per channel/chip, to review and merge into
conf/lbc_detector.ecsv) and run_log.json (code version, parameters,
checksums). Installing the rows into the package table is a separate step.
"""
import argparse, datetime, hashlib, json, os, subprocess
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
import LBCgo
from LBCgo import detector

HERE = Path(__file__).resolve().parent

# Parameters that define this product (change -> new product directory)
PARAMS = {'box': 1000, 'sigma': 4.0, 'mjd_start': np.nan, 'mjd_end': np.nan}

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
    ap.add_argument('--raw', default=os.environ.get('LBCGO_RAW'))
    ap.add_argument('--box', type=int, default=PARAMS['box'])
    args = ap.parse_args()
    if not args.raw:
        ap.error('give --raw or set LBCGO_RAW')
    raw, params = Path(args.raw), dict(PARAMS, box=args.box)

    inputs = Table.read(HERE / 'inputs.ecsv', format='ascii.ecsv')
    results = []
    for set_id in sorted(set(inputs['set'])):
        rows = inputs[inputs['set'] == set_id]
        flats = [raw / f for f in rows[rows['role'] == 'flat']['filename']]
        biases = [raw / f for f in rows[rows['role'] == 'bias']['filename']]
        if len(flats) != 2 or len(biases) != 2:
            raise ValueError(f'set {set_id}: need 2 flats and 2 biases')
        # inputs.ecsv must agree with the headers
        channels = {detector.lbc_channel(fits.getheader(p), str(p))
                    for p in flats + biases}
        if channels != set(rows['channel']):
            raise ValueError(f'set {set_id}: channel mismatch {channels}')

        t = detector.measure_gain_rdnoise_files(
            *map(str, flats + biases), box=params['box'], sigma=params['sigma'])
        t['set'] = set_id
        t['mjd'] = float(np.mean(rows['mjd_obs']))
        results.append(t)

    results = vstack(results)
    results.write(HERE / 'results_per_set.ecsv', overwrite=True)

    # One product row per channel/chip: median over sets
    commit = git_commit()
    product = Table(names=detector.TABLE_COLUMNS,
                    dtype=['U4', 'i4', 'f8', 'f8', 'f8', 'f8', 'U200'])
    source = f'calibration/{HERE.name} (LBCgo {commit[:10]})'
    for channel in sorted(set(results['channel'])):
        for chip in sorted(set(results['chip'])):
            sel = (results['channel'] == channel) & (results['chip'] == chip)
            product.add_row((channel, chip,
                             float(np.median(results['gain'][sel])),
                             float(np.median(results['rdnoise'][sel])),
                             params['mjd_start'], params['mjd_end'], source))
    detector.write_detector_table(product, str(HERE / 'detector_rows.ecsv'),
                                  overwrite=True)

    log = {'date_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
           'lbcgo_version': LBCgo.__version__, 'lbcgo_commit': commit,
           'params': {k: None if isinstance(v, float) and np.isnan(v) else v
                      for k, v in params.items()},
           'inputs_sha256': {str(f): sha256(raw / f) for f in inputs['filename']},
           'outputs_sha256': {n: sha256(HERE / n) for n in
                              ('results_per_set.ecsv', 'detector_rows.ecsv')}}
    (HERE / 'run_log.json').write_text(json.dumps(log, indent=2) + '\n')

if __name__ == '__main__':
    main()