"""LBC detector parameters: channel identification, gain/read-noise lookup
and photon-transfer measurement.

The ``GAIN`` and ``RDNOISE`` keywords in LBC raw headers are the same nominal
values (1.75 e-/ADU, 12 e-) on every chip of both LBC-Blue and LBC-Red, so
they are not per-chip measurements. Measured values can be supplied in a
per-chip table (ECSV) with these columns:

    channel    'LBCB' or 'LBCR'
    chip       1-4
    gain       e-/ADU
    rdnoise    e-
    mjd_start  first MJD the row applies to (inclusive)
    mjd_end    last MJD the row applies to (exclusive)
    source     free text (e.g. which frames were used)

:func:`gain_rdnoise` looks the values up in that table and falls back to the
image headers when the table is absent or has no matching row.
:func:`measure_gain_rdnoise_files` measures them from two raw flats and two
raw biases, and :func:`write_detector_table` writes the table.

The packaged table ``conf/lbc_detector.ecsv`` is used by default; it ships
empty, so header values are used until measured values are added.
"""

import os

import numpy as np
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
import astropy.units as u
from ccdproc.utils.slices import slice_from_string

from . import masks as lbcmasks

DEFAULT_TABLE = os.path.join(os.path.dirname(__file__), 'conf',
                             'lbc_detector.ecsv')

TABLE_COLUMNS = ['channel', 'chip', 'gain', 'rdnoise', 'mjd_start',
                 'mjd_end', 'source']


def lbc_channel(header, filename=None):
    """Return 'LBCB', 'LBCR' or None for a header (and optional filename).

    Raw headers spell INSTRUME inconsistently ('LBC_BLUE' vs 'LBC-RED '),
    so this normalizes INSTRUME, then tries DETECTOR ('EEV-BLUE'/'EEV-RED'),
    then the 'lbcb.'/'lbcr.' filename prefix.
    """
    for key in ('INSTRUME', 'DETECTOR'):
        value = str(header.get(key, '')).upper()
        if 'BLUE' in value:
            return 'LBCB'
        if 'RED' in value:
            return 'LBCR'
    if filename is not None:
        base = os.path.basename(filename).lower()
        if base.startswith('lbcb'):
            return 'LBCB'
        if base.startswith('lbcr'):
            return 'LBCR'
    return None


def read_detector_table(filename=None):
    """Read a detector table.

    With ``filename=None`` the packaged default is read (None if it is
    missing). An explicitly given file that does not exist raises
    FileNotFoundError rather than silently falling back to header values.
    """
    if filename is None:
        if not os.path.exists(DEFAULT_TABLE):
            return None
        filename = DEFAULT_TABLE
    elif not os.path.exists(filename):
        raise FileNotFoundError('Detector table {0} not found.'.format(filename))
    return Table.read(filename, format='ascii.ecsv')


def lookup_detector_params(table, channel, chip, mjd=None):
    """Return (gain, rdnoise) from ``table`` for a channel/chip/date, or None.

    If several rows match, the one with the latest ``mjd_start`` wins. A row
    matches an unknown date (``mjd=None``) only if it covers all dates.
    """
    if table is None or channel is None or len(table) == 0:
        return None
    rows = table[(np.char.upper(np.asarray(table['channel'], dtype=str))
                  == channel) & (table['chip'] == chip)]
    start = np.asarray(rows['mjd_start'], dtype=float)
    end = np.asarray(rows['mjd_end'], dtype=float)
    if mjd is None:
        rows = rows[~np.isfinite(start) & ~np.isfinite(end)]
    else:
        start = np.where(np.isfinite(start), start, -np.inf)
        end = np.where(np.isfinite(end), end, np.inf)
        rows = rows[(start <= mjd) & (mjd < end)]
    if len(rows) == 0:
        return None
    start = np.asarray(rows['mjd_start'], dtype=float)
    row = rows[np.argmax(np.where(np.isfinite(start), start, -np.inf))]
    return float(row['gain']), float(row['rdnoise'])


def gain_rdnoise(chip, headers, table=None, filename=None, verbose=False):
    """Gain (e-/ADU) and read noise (e-) for one chip, with their source.

    Order of precedence: matching row of ``table`` (a Table, or the packaged
    default if None), then the GAIN/RDNOISE header keywords, then the
    LBC-Blue nominal defaults in :mod:`LBCgo.masks`.

    Parameters
    ----------
    chip : int
        Chip number (1-4).
    headers : list of fits.Header
        Headers to search, most specific first (chip header, then primary).
    table : astropy.table.Table or None, optional
        Detector table. If None, the packaged default is read.
    filename : str, optional
        Image filename, used only to identify the channel.

    Returns
    -------
    gain, rdnoise : float
    source : {'table', 'header', 'default'}
    """
    if table is None:
        table = read_detector_table()

    channel = None
    mjd = None
    for hdr in headers:
        if hdr is None:
            continue
        channel = channel or lbc_channel(hdr)
        if mjd is None and 'MJD_OBS' in hdr:
            mjd = float(hdr['MJD_OBS'])
    if channel is None and filename is not None:
        channel = lbc_channel({}, filename)

    found = lookup_detector_params(table, channel, chip, mjd)
    if found is not None:
        return found[0], found[1], 'table'

    have_header = any(h is not None and 'GAIN' in h and 'RDNOISE' in h
                      for h in headers)
    gain = lbcmasks.header_value('GAIN', headers, lbcmasks.DEFAULT_GAIN,
                                 verbose=verbose)
    rdnoise = lbcmasks.header_value('RDNOISE', headers,
                                    lbcmasks.DEFAULT_RDNOISE, verbose=verbose)
    return gain, rdnoise, ('header' if have_header else 'default')


def measure_gain_rdnoise(flat1, flat2, bias1, bias2, sigma=4.0):
    """Photon-transfer gain and read noise from two flats and two biases.

    All inputs are arrays in ADU over the same pixels (e.g. overscan-
    subtracted, trimmed central regions). With mu_i the bias-subtracted mean
    of flat i, k = mu_1 / mu_2 and sigma_r^2 = var(B1 - B2) / 2 the read
    noise variance in ADU^2, the difference F1 - k*F2 cancels the fixed
    pattern (flat-field structure) and has variance
    mu_1 (1 + k) / g + (1 + k^2) sigma_r^2, so

        gain    = mu_1 (1 + k) / (var(F1 - k F2) - (1 + k^2) sigma_r^2)
        rdnoise = gain * sigma_r

    For equal flats (k = 1) this is the usual
    (mu_1 + mu_2) / (var(F1 - F2) - var(B1 - B2)). Variances are
    sigma-clipped to reject cosmic rays and defects.

    Returns
    -------
    gain : float
        e-/ADU
    rdnoise : float
        e-
    """
    f1 = np.asarray(flat1, dtype=float)
    f2 = np.asarray(flat2, dtype=float)
    b1 = np.asarray(bias1, dtype=float)
    b2 = np.asarray(bias2, dtype=float)

    _, m_b1, _ = sigma_clipped_stats(b1, sigma=sigma)
    _, m_b2, _ = sigma_clipped_stats(b2, sigma=sigma)
    bias_level = 0.5 * (m_b1 + m_b2)
    _, mu1, _ = sigma_clipped_stats(f1 - bias_level, sigma=sigma)
    _, mu2, _ = sigma_clipped_stats(f2 - bias_level, sigma=sigma)
    if mu1 <= 0 or mu2 <= 0:
        raise ValueError('Flats must be brighter than the biases.')
    k = mu1 / mu2

    _, _, s_ff = sigma_clipped_stats((f1 - bias_level) - k * (f2 - bias_level),
                                     sigma=sigma)
    _, _, s_bb = sigma_clipped_stats(b1 - b2, sigma=sigma)
    read_var = 0.5 * s_bb ** 2

    noise = s_ff ** 2 - (1.0 + k ** 2) * read_var
    if noise <= 0:
        raise ValueError('Flat difference is not noisier than the read '
                         'noise; flats too faint?')
    gain = mu1 * (1.0 + k) / noise
    rdnoise = gain * np.sqrt(read_var)
    return float(gain), float(rdnoise)


def _overscan_corrected_chip(filename, chip, box):
    """Overscan-subtracted, trimmed data of one chip (central ``box``)."""
    ccd = CCDData.read(filename, chip, unit=u.adu)
    data = np.asarray(ccd.data, dtype=float)
    bias = slice_from_string(ccd.header['BIASSEC'], fits_convention=True)
    trim = slice_from_string(ccd.header['TRIMSEC'], fits_convention=True)
    data = data[trim] - np.median(data[bias])
    if box is not None:
        ny, nx = data.shape
        hy, hx = min(box, ny) // 2, min(box, nx) // 2
        data = data[ny // 2 - hy:ny // 2 + hy, nx // 2 - hx:nx // 2 + hx]
    return data, ccd.header


def measure_gain_rdnoise_files(flat1, flat2, bias1, bias2, lbc_chips=True,
                               box=1000, sigma=4.0):
    """Measure gain and read noise per chip from raw LBC files.

    Parameters
    ----------
    flat1, flat2 : str
        Two raw flats of the same channel at similar, unsaturated levels
        (ideally 10,000-30,000 ADU above bias).
    bias1, bias2 : str
        Two raw bias frames of the same channel and readout mode.
    lbc_chips : bool or list of int, optional
        Chips to measure. Default: all four.
    box : int or None, optional
        Use a central ``box`` x ``box`` region of each trimmed chip, avoiding
        vignetted edges. None uses the whole chip. Default: 1000
    sigma : float, optional
        Clipping threshold. Default: 4

    Returns
    -------
    astropy.table.Table
        One row per chip with the :data:`TABLE_COLUMNS` columns. mjd_start
        and mjd_end are left open (NaN); edit them before use if the values
        should apply to a limited date range.
    """
    if lbc_chips is True:
        lbc_chips = [1, 2, 3, 4]
    primary = fits.getheader(flat1)
    channel = lbc_channel(primary, flat1)

    rows = []
    for chip in lbc_chips:
        f1, chdr = _overscan_corrected_chip(flat1, chip, box)
        f2, _ = _overscan_corrected_chip(flat2, chip, box)
        b1, _ = _overscan_corrected_chip(bias1, chip, box)
        b2, _ = _overscan_corrected_chip(bias2, chip, box)
        saturate = lbcmasks.header_value('SATURATE', [chdr, primary],
                                         lbcmasks.DEFAULT_SATURATE)
        if max(np.median(f1), np.median(f2)) > 0.7 * saturate:
            raise ValueError('Chip {0}: flats are too close to saturation '
                             'for a photon-transfer measurement.'.format(chip))
        gain, rdnoise = measure_gain_rdnoise(f1, f2, b1, b2, sigma=sigma)
        source = 'PTC: {0}, {1}; {2}, {3}'.format(
            *(os.path.basename(f) for f in (flat1, flat2, bias1, bias2)))
        rows.append((channel or '', chip, gain, rdnoise, np.nan, np.nan,
                     source))
    table = Table(rows=rows, names=TABLE_COLUMNS)
    table['gain'].unit = u.electron / u.adu
    table['rdnoise'].unit = u.electron
    return table


def write_detector_table(table, filename=None, overwrite=False):
    """Write a detector table (default: the packaged ``lbc_detector.ecsv``)."""
    if filename is None:
        filename = DEFAULT_TABLE
    table.write(filename, format='ascii.ecsv', overwrite=overwrite)
    return filename
