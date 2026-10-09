"""Build a calibration inputs.ecsv from a directory of raw LBC frames.

Uses ccdproc's ImageFileCollection, the same way lbcgo() organises raw data
(same file globs, header keywords, and SkyFlatTest exclusion), and writes one
row per frame with the columns calibration/NEW_PRODUCT_TEMPLATE.md asks for:

    filename obs_id mjd_obs channel filter exptime role set notes

plus object, propid, lbcobnam, imagetyp, airmass for reviewing the selection.

`role` comes from IMAGETYP (flat -> flat, zero/bias -> bias, object ->
science). `set` numbers groups of frames from the same channel and UT date;
for a gain/read-noise product, edit the table so each set holds exactly two
flats and two biases (--levels helps choose flats at similar levels).

Examples:
    # all raw LBC frames in ./raw
    python calibration/make_inputs.py raw -o my_product/inputs.ecsv

    # LBCR sky flats and biases, with flat levels to pick pairs
    python calibration/make_inputs.py raw --channel LBCR \\
        --imagetyp flat zero --levels -o gain_rdnoise_lbcr/inputs.ecsv

    # science frames of one target in two filters
    python calibration/make_inputs.py raw --imagetyp object \\
        --object NGC891 --filter g-SLOAN r-SLOAN -o v4/inputs.ecsv
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from ccdproc.utils.slices import slice_from_string

from LBCgo.detector import lbc_channel

# Header keywords read from the primary header (lower case, as ccdproc
# reports them).
KEYWORDS = ['object', 'filter', 'exptime', 'imagetyp', 'propid', 'lbcobnam',
            'airmass', 'mjd_obs', 'obs_id', 'instrume', 'detector']

GLOBS = {'LBCB': 'lbcb.*.*.fits*', 'LBCR': 'lbcr.*.*.fits*'}

ROLES = {'flat': 'flat', 'zero': 'bias', 'bias': 'bias', 'object': 'science'}

COLUMNS = ['filename', 'obs_id', 'mjd_obs', 'channel', 'filter', 'exptime',
           'role', 'set', 'notes', 'object', 'propid', 'lbcobnam',
           'imagetyp', 'airmass']


def _scan_headers(directory, pattern, keywords=KEYWORDS):
    """Primary-header summary of the frames in `directory` matching `pattern`.

    Replaces ccdproc's ImageFileCollection, which parses *every* card in each
    header and so fails on frames with an invalid card (e.g. 2010-03-18 LBCB
    frames have an unquoted ``PA_PNT = nan``). Here only the requested
    keywords are parsed. Returns a list of dicts keyed by lower-case keyword,
    plus 'file'; keywords absent from (or unparsable in) a header are omitted.
    """
    rows = []
    for path in sorted(directory.glob(pattern)):
        header = fits.getheader(path)
        row = {'file': path.name}
        for key in keywords:
            try:
                value = header[key]
            except (KeyError, fits.verify.VerifyError):
                continue
            if isinstance(value, (str, int, float, np.number)):
                row[key] = value
        rows.append(row)
    return rows


def _value(row, key, default=''):
    """Header value from a summary row, with missing -> default."""
    if key not in row:
        return default
    value = row[key]
    return value.strip() if isinstance(value, str) else value


def _matches(value, wanted):
    """Case-insensitive, whitespace-insensitive match against a list."""
    if not wanted:
        return True
    return str(value).strip().lower() in {w.strip().lower() for w in wanted}


def chip_level(path, chip=2, box=500):
    """Overscan-subtracted median of a central box of one chip [ADU]."""
    with fits.open(path) as hdul:
        hdr = hdul[chip].header
        data = hdul[chip].data.astype(float)
    data = (data[slice_from_string(hdr['TRIMSEC'], fits_convention=True)]
            - np.median(data[slice_from_string(hdr['BIASSEC'],
                                               fits_convention=True)]))
    ny, nx = data.shape
    hy, hx = min(box, ny) // 2, min(box, nx) // 2
    return float(np.median(data[ny // 2 - hy:ny // 2 + hy,
                                nx // 2 - hx:nx // 2 + hx]))


def build_inputs(directory, channels=('LBCB', 'LBCR'), imagetyp=None,
                 objects=None, filters=None, propids=None, levels=False,
                 keep_tests=False):
    """Return the inputs table for the raw frames in `directory`."""
    directory = Path(directory)
    rows = []
    for channel in channels:
        summary = _scan_headers(directory, GLOBS[channel])
        for row in summary:
            lbcobnam = _value(row, 'lbcobnam')
            if not keep_tests and str(lbcobnam).startswith('SkyFlatTest'):
                continue
            typ = str(_value(row, 'imagetyp')).lower()
            if not (_matches(typ, imagetyp)
                    and _matches(_value(row, 'object'), objects)
                    and _matches(_value(row, 'filter'), filters)
                    and _matches(_value(row, 'propid'), propids)):
                continue
            header = {'INSTRUME': _value(row, 'instrume'),
                      'DETECTOR': _value(row, 'detector')}
            found = lbc_channel(header, row['file'])
            notes = ''
            if found is not None and found != channel:
                notes = f'header says {found}; filename says {channel}'
            rows.append({
                'filename': row['file'],
                'obs_id': _value(row, 'obs_id'),
                'mjd_obs': float(_value(row, 'mjd_obs', np.nan)),
                'channel': channel,
                'filter': _value(row, 'filter'),
                'exptime': float(_value(row, 'exptime', np.nan)),
                'role': ROLES.get(typ, typ or 'unknown'),
                'set': 0,
                'notes': notes,
                'object': _value(row, 'object'),
                'propid': _value(row, 'propid'),
                'lbcobnam': lbcobnam,
                'imagetyp': typ,
                'airmass': float(_value(row, 'airmass', np.nan)),
            })

    if not rows:
        return Table(names=COLUMNS)
    table = Table(rows=[[r[c] for c in COLUMNS] for r in rows],
                  names=COLUMNS)
    table.sort(['channel', 'mjd_obs', 'filename'])

    # set = running number per (channel, UT date)
    keys = [(c, int(np.floor(m)) if np.isfinite(m) else -1)
            for c, m in zip(table['channel'], table['mjd_obs'])]
    numbering = {k: i + 1 for i, k in enumerate(dict.fromkeys(keys))}
    table['set'] = [numbering[k] for k in keys]

    if levels:
        table['median_adu'] = [
            chip_level(directory / f) if r in ('flat', 'bias') else np.nan
            for f, r in zip(table['filename'], table['role'])]
    return table


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='\n'.join(__doc__.splitlines()[1:]))
    ap.add_argument('directory', nargs='?', default='.',
                    help='directory with raw lbcb.*/lbcr.* frames')
    ap.add_argument('-o', '--output', default='inputs.ecsv')
    ap.add_argument('--channel', nargs='+', choices=['LBCB', 'LBCR'],
                    default=['LBCB', 'LBCR'])
    ap.add_argument('--imagetyp', nargs='+',
                    help='e.g. flat zero object (default: all)')
    ap.add_argument('--object', nargs='+', help='OBJECT values to keep')
    ap.add_argument('--filter', nargs='+', help='FILTER values to keep')
    ap.add_argument('--propid', nargs='+', help='PROPID values to keep')
    ap.add_argument('--levels', action='store_true',
                    help='add median_adu: overscan-subtracted median of the '
                         'centre of chip 2 (flats and biases only)')
    ap.add_argument('--keep-tests', action='store_true',
                    help="keep LBCOBNAM 'SkyFlatTest*' frames (lbcgo drops "
                         "them)")
    ap.add_argument('--overwrite', action='store_true')
    args = ap.parse_args(argv)

    out = Path(args.output)
    if out.exists() and not args.overwrite:
        ap.error(f'{out} exists; use --overwrite')

    table = build_inputs(args.directory, channels=args.channel,
                         imagetyp=args.imagetyp, objects=args.object,
                         filters=args.filter, propids=args.propid,
                         levels=args.levels, keep_tests=args.keep_tests)
    if len(table) == 0:
        print('No frames matched.', file=sys.stderr)
        return 1
    table.meta['comments'] = [
        'Calibration inputs built by calibration/make_inputs.py from '
        f'{Path(args.directory).resolve().name}/.',
        'role from IMAGETYP; set = channel + UT date. Review and edit before '
        'use (e.g. 2 flats + 2 biases per set for gain/read noise).']
    table.write(out, format='ascii.ecsv', overwrite=args.overwrite)
    counts = ', '.join(f'{r}: {np.sum(table["role"] == r)}'
                       for r in sorted(set(table['role'])))
    print(f'Wrote {len(table)} frames to {out} ({counts}; '
          f'{len(set(table["set"]))} sets)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
