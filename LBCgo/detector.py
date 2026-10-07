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

The packaged table ``conf/lbc_detector.ecsv`` is used by default. It is
seeded with the published LBC-Blue values (Giallongo et al. 2008, A&A 482,
349, Table 1; 2006 commissioning, open-ended validity) until measured values
replace them. It has no LBC-Red rows, so LBC-Red uses header values.
"""

import os
import warnings

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

PTC_COLUMNS = ['level', 'rdnoise_adu', 'k', 'rho_x', 'rho_y', 'bias_rho_x',
               'bias_rho_y', 'gain_nn', 'rho_sum', 'rho_sum_err',
               'bias_rho_sum', 'gain_sum', 'gain_sum_err']

TABLE_COLUMNS = ['channel', 'chip', 'gain', 'rdnoise', 'mjd_start',
                 'mjd_end', 'source']

# Optional detector-table columns (NaN when absent or unknown).
# gain_flux: electrons per ADU of a flux summed over several pixels, i.e.
# the gain with the pixel-to-pixel covariances included (median gain_sum of
# measure_ptc over suitable flat pairs). It differs from the per-pixel
# ``gain`` where the readout correlates pixels even at zero signal
# (charge-transfer inefficiency, video-chain undershoot). Use it for
# Poisson errors of source fluxes; ``gain`` sets the per-pixel variance
# (weight maps).
OPTIONAL_TABLE_COLUMNS = ['gain_flux']


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
    table = Table.read(filename, format='ascii.ecsv')
    for col in OPTIONAL_TABLE_COLUMNS:
        if col not in table.colnames:
            table[col] = np.full(len(table), np.nan)
            table[col].unit = u.electron / u.adu
    return table


def _lookup_row(table, channel, chip, mjd=None):
    """The matching row of ``table`` (see :func:`lookup_detector_params`)."""
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
    start = np.where(np.isfinite(start), start, -np.inf)
    best = start == start.max()
    if best.sum() > 1:
        warnings.warn('{0} chip {1}: {2} detector-table rows match with the '
                      'same mjd_start ({3}); using the first. Give the newer '
                      'row a finite mjd_start or end the older one.'.format(
                          channel, chip, int(best.sum()),
                          'open' if np.isinf(start.max()) else start.max()),
                      UserWarning, stacklevel=3)
    return rows[np.argmax(best)]


def lookup_detector_params(table, channel, chip, mjd=None):
    """Return (gain, rdnoise) from ``table`` for a channel/chip/date, or None.

    If several rows match, the one with the latest ``mjd_start`` wins, so a
    newer row with a finite ``mjd_start`` supersedes an open-ended older one.
    If the latest ``mjd_start`` is shared by more than one matching row
    (e.g. two open-ended rows), the choice is ambiguous: the first such row
    is used and a ``UserWarning`` is issued (see
    :func:`detector_table_conflicts`). A row matches an unknown date
    (``mjd=None``) only if it covers all dates.
    """
    row = _lookup_row(table, channel, chip, mjd)
    if row is None:
        return None
    return float(row['gain']), float(row['rdnoise'])


def lookup_gain_flux(table, channel, chip, mjd=None):
    """Flux gain (e-/ADU) for a channel/chip/date, with its source.

    Same row selection as :func:`lookup_detector_params`. Returns
    ``(gain_flux, 'gain_flux')`` if the row has a finite ``gain_flux``,
    ``(gain, 'gain')`` (the per-pixel gain) if it does not, and None if no
    row matches. See :data:`OPTIONAL_TABLE_COLUMNS` for the difference.
    """
    row = _lookup_row(table, channel, chip, mjd)
    if row is None:
        return None
    if 'gain_flux' in row.colnames and np.isfinite(float(row['gain_flux'])):
        return float(row['gain_flux']), 'gain_flux'
    return float(row['gain']), 'gain'


def detector_table_conflicts(table):
    """Rows that make :func:`lookup_detector_params` ambiguous.

    Returns a list of ``(channel, chip, mjd_start)`` for every group of two
    or more rows sharing channel, chip and ``mjd_start`` (NaN counts as one
    value): such rows always overlap in time and neither supersedes the
    other. An empty list means the table is unambiguous.
    """
    if table is None or len(table) == 0:
        return []
    counts = {}
    for row in table:
        start = float(row['mjd_start'])
        key = (str(row['channel']).strip().upper(), int(row['chip']),
               None if not np.isfinite(start) else start)
        counts[key] = counts.get(key, 0) + 1
    return sorted((k for k, n in counts.items() if n > 1),
                  key=lambda k: (k[0], k[1], -np.inf if k[2] is None else k[2]))


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


def _blocks(data, cell):
    """Split ``data`` into ``cell`` x ``cell`` blocks, shape (n, cell, cell).

    ``cell=None``, or a cell larger than the array, gives the whole array as
    one block. Edge rows/columns that do not fill a block are dropped.
    """
    data = np.asarray(data, dtype=float)
    ny, nx = data.shape
    if cell is None or cell >= min(ny, nx):
        return data[np.newaxis]
    ncy, ncx = ny // cell, nx // cell
    return (data[:ncy * cell, :ncx * cell]
            .reshape(ncy, cell, ncx, cell)
            .transpose(0, 2, 1, 3)
            .reshape(ncy * ncx, cell, cell))


def _cell_variance(data, cell, sigma):
    """Median over ``cell`` x ``cell`` blocks of the clipped variance.

    Each block's own mean is removed, so structure on scales larger than a
    block (illumination gradients, slowly varying bias) does not add
    variance. ``cell=None``, or a cell larger than the array, uses the whole
    array as one block.
    """
    data = np.asarray(data, dtype=float)
    ny, nx = data.shape
    if cell is None or cell >= min(ny, nx):
        _, _, std = sigma_clipped_stats(data, sigma=sigma)
        return float(std ** 2)
    blocks = _blocks(data, cell).reshape(-1, cell * cell)
    _, _, std = sigma_clipped_stats(blocks, sigma=sigma, axis=1)
    return float(np.nanmedian(np.asarray(std) ** 2))


def _cell_correlations(data, cell, sigma):
    """Nearest-neighbour correlation coefficients (rho_x, rho_y).

    In each ``cell`` x ``cell`` block (whole array if ``cell`` is None or too
    large) the clipped mean is removed and pixels beyond ``sigma`` are
    ignored; rho is the lag-1 covariance along x (y) divided by the
    variance. Returns the medians over blocks. Independent pixels give 0;
    the brighter-fatter effect gives positive values that grow with signal.
    """
    blocks = _blocks(data, cell)
    mean, _, std = sigma_clipped_stats(blocks, sigma=sigma, axis=(1, 2))
    mean = np.asarray(mean)[:, None, None]
    std = np.asarray(std)[:, None, None]
    d = blocks - mean
    d = np.where(np.abs(d) > sigma * std, np.nan, d)
    var = np.nanmean(d ** 2, axis=(1, 2))
    rho_x = np.nanmean(d[:, :, :-1] * d[:, :, 1:], axis=(1, 2)) / var
    rho_y = np.nanmean(d[:, :-1, :] * d[:, 1:, :], axis=(1, 2)) / var
    return float(np.nanmedian(rho_x)), float(np.nanmedian(rho_y))


def _cell_covariance_sum(data, cell, sigma, max_lag):
    """Sum of the correlation coefficients over all lags within ``max_lag``.

    S = sum of rho(dx, dy) over 0 < max(|dx|, |dy|) <= max_lag, i.e.
    (2 max_lag + 1)^2 - 1 lags, so that var * (1 + S) is the variance
    summed over that neighbourhood. The brighter-fatter effect moves charge
    between pixels without changing its total, so var * (1 + S) is free of
    it once ``max_lag`` covers the range of the effect; non-linearity
    creates no covariances and survives in it.

    In each block a plane (not just the mean) is removed, because a
    residual gradient correlates all lags and the sum over many lags would
    amplify it. Removing p = 3 fitted parameters biases each lag by about
    -p (1 + S) / n (n valid pixels in the block); that is added back.
    Pixels beyond ``sigma`` (clipped statistics of the block) are ignored.

    Returns
    -------
    S, S_err, rho : float, float, ndarray
        Clipped mean over blocks of the per-block sum (the median is 25 %
        noisier, and this sum is noise-limited); its uncertainty
        (std / sqrt(n_blocks), NaN for a single block); and the median
        correlation map,
        shape (2 max_lag + 1, 2 max_lag + 1), indexed [dy + max_lag,
        dx + max_lag], with rho[max_lag, max_lag] = 1.
    """
    blocks = _blocks(data, cell)
    nb, cy, cx = blocks.shape
    mean, _, std = sigma_clipped_stats(blocks, sigma=sigma, axis=(1, 2))
    mean = np.asarray(mean)[:, None, None]
    std = np.asarray(std)[:, None, None]
    good = np.abs(blocks - mean) <= sigma * std

    # Least-squares plane a + b x + c y per block, on the good pixels
    yy, xx = np.mgrid[:cy, :cx]
    design = np.stack([np.ones((cy, cx)), xx - (cx - 1) / 2.0,
                       yy - (cy - 1) / 2.0]).reshape(3, -1)
    w = good.reshape(nb, -1).astype(float)
    z = np.where(good, blocks, 0.0).reshape(nb, -1)
    normal = np.einsum('pk,qk,bk->bpq', design, design, w)
    rhs = np.einsum('pk,bk->bp', design, w * z)
    coef = np.linalg.solve(normal, rhs[..., None])[..., 0]
    d = blocks - (coef @ design).reshape(nb, cy, cx)
    d = np.where(good, d, np.nan)

    n = good.sum(axis=(1, 2))
    var = np.nanmean(d ** 2, axis=(1, 2))
    size = 2 * max_lag + 1
    rho = np.full((nb, size, size), np.nan)
    rho[:, max_lag, max_lag] = 1.0
    for dy in range(0, max_lag + 1):
        for dx in range(-max_lag, max_lag + 1):
            if dy == 0 and dx <= 0:
                continue
            a = d[:, dy:, max(dx, 0):cx + min(dx, 0)]
            b = d[:, :cy - dy, max(-dx, 0):cx - max(dx, 0)]
            r = np.nanmean(a * b, axis=(1, 2)) / var
            rho[:, max_lag + dy, max_lag + dx] = r
            rho[:, max_lag - dy, max_lag - dx] = r      # rho(-l) = rho(l)
    # Bias of each lag from the plane fit: -c (1 + S_total), c = 3 / n.
    # Taking S_total ~ S (covariances beyond max_lag neglected) and solving
    # S = S_raw + n_lags c (1 + S) for S:
    c = 3.0 / n
    n_lags = size * size - 1
    raw = np.nansum(rho, axis=(1, 2)) - 1.0
    per_block = (raw + n_lags * c) / (1.0 - n_lags * c)
    lagged = np.ones((size, size), dtype=bool)
    lagged[max_lag, max_lag] = False
    rho[:, lagged] += (c * (1.0 + per_block))[:, None]
    if nb == 1:
        return float(per_block[0]), np.nan, rho[0]
    S, _, S_std = sigma_clipped_stats(per_block, sigma=sigma)
    n_used = np.count_nonzero(np.abs(per_block - S) <= sigma * S_std)
    return (float(S), float(S_std / np.sqrt(n_used)),
            np.nanmedian(rho, axis=0))


def measure_ptc(flat1, flat2, bias1, bias2, sigma=4.0, cell=50,
                max_lag=3):
    """Photon-transfer measurement with diagnostics.

    Same method as :func:`measure_gain_rdnoise` (see there), returning a
    dict with:

    gain, rdnoise
        e-/ADU and e-.
    level
        Mean bias-subtracted signal of the two flats [ADU]. The apparent
        gain can depend on it (brighter-fatter effect, non-linearity), so
        measurements at several levels are fitted with
        :func:`fit_gain_vs_level`.
    rdnoise_adu
        Read noise [ADU]; independent of the gain.
    k
        mu_1 / mu_2.
    rho_x, rho_y
        Nearest-neighbour correlation coefficients of the flat difference
        (:func:`_cell_correlations`). Positive values growing with level
        indicate the brighter-fatter effect: charge pushed into
        neighbouring pixels lowers the per-pixel variance and raises the
        apparent gain.
    bias_rho_x, bias_rho_y
        The same for the bias difference (electronic correlations).
    gain_nn
        Diagnostic gain with the nearest-neighbour covariances added back to
        the variances, var * (1 + 2 rho_x + 2 rho_y). If the level
        dependence comes from the brighter-fatter effect this is much less
        level-dependent than ``gain``. It ignores longer-range covariances,
        so it is a test of the explanation, not the adopted value.
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

    diff_ff = (f1 - bias_level) - k * (f2 - bias_level)
    diff_bb = b1 - b2
    var_ff = _cell_variance(diff_ff, cell, sigma)
    read_var = 0.5 * _cell_variance(diff_bb, cell, sigma)

    noise = var_ff - (1.0 + k ** 2) * read_var
    if noise <= 0:
        raise ValueError('Flat difference is not noisier than the read '
                         'noise; flats too faint?')
    gain = mu1 * (1.0 + k) / noise

    rho_x, rho_y = _cell_correlations(diff_ff, cell, sigma)
    brho_x, brho_y = _cell_correlations(diff_bb, cell, sigma)
    noise_nn = (var_ff * (1 + 2 * (rho_x + rho_y))
                - (1.0 + k ** 2) * read_var * (1 + 2 * (brho_x + brho_y)))
    gain_nn = mu1 * (1.0 + k) / noise_nn if noise_nn > 0 else np.nan

    if max_lag and max_lag > 0:
        S, S_err, _ = _cell_covariance_sum(diff_ff, cell, sigma, max_lag)
        bS, _, _ = _cell_covariance_sum(diff_bb, cell, sigma, max_lag)
        noise_sum = var_ff * (1 + S) - (1.0 + k ** 2) * read_var * (1 + bS)
        if noise_sum > 0:
            gain_sum = mu1 * (1.0 + k) / noise_sum
            gain_sum_err = gain_sum * var_ff * S_err / noise_sum
        else:
            gain_sum = gain_sum_err = np.nan
    else:
        S = S_err = bS = gain_sum = gain_sum_err = np.nan

    return {'gain': float(gain), 'rdnoise': float(gain * np.sqrt(read_var)),
            'level': float(0.5 * (mu1 + mu2)),
            'rdnoise_adu': float(np.sqrt(read_var)), 'k': float(k),
            'rho_x': rho_x, 'rho_y': rho_y,
            'bias_rho_x': brho_x, 'bias_rho_y': brho_y,
            'gain_nn': float(gain_nn), 'rho_sum': float(S),
            'rho_sum_err': float(S_err), 'bias_rho_sum': float(bS),
            'gain_sum': float(gain_sum), 'gain_sum_err': float(gain_sum_err)}


def measure_gain_rdnoise(flat1, flat2, bias1, bias2, sigma=4.0, cell=50):
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

    The variances are measured in ``cell`` x ``cell`` blocks (each with its
    own mean removed) and the median over blocks is used. Real twilight
    flats differ slightly in illumination (sky gradient, exposure time,
    shutter, rotator angle); any difference that k cannot remove adds
    variance to F1 - k F2 and biases the gain *low*. On synthetic 1000-px
    regions a 0.5 % (1 %) peak-to-peak mismatch biased a whole-region
    variance by -6 % (-22 %) in gain; 50-px blocks recover the gain to
    < 0.5 % for mismatches up to 2 %. Structure on scales comparable to a
    block (e.g. z-band fringing) is not removed: avoid such flats.
    ``cell=None`` reproduces the whole-region estimate.

    Parameters
    ----------
    flat1, flat2, bias1, bias2 : ndarray
        Images in ADU.
    sigma : float, optional
        Clipping threshold. Default: 4
    cell : int or None, optional
        Block size in pixels for the variances. Default: 50

    Returns
    -------
    gain : float
        e-/ADU
    rdnoise : float
        e-
    """
    m = measure_ptc(flat1, flat2, bias1, bias2, sigma=sigma, cell=cell,
                    max_lag=0)
    return m['gain'], m['rdnoise']


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
                               box=1000, sigma=4.0, cell=50, max_lag=3):
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
    cell : int or None, optional
        Block size for the variances (see :func:`measure_gain_rdnoise`).
        Default: 50
    max_lag : int, optional
        Lag range of the covariance sum behind ``gain_sum`` (see
        :func:`measure_ptc`); 0 skips it. Default: 3

    Use pairs of consecutive flats from one sequence (same filter, rotator
    angle and similar exposure time). Avoid z/Y-band twilight flats, whose
    fringing changes through twilight on block-sized scales.

    Returns
    -------
    astropy.table.Table
        One row per chip with the :data:`TABLE_COLUMNS` columns plus the
        diagnostics of :func:`measure_ptc` (:data:`PTC_COLUMNS`: level,
        rdnoise_adu, k, rho_x, rho_y, bias_rho_x, bias_rho_y, gain_nn,
        rho_sum, rho_sum_err, bias_rho_sum, gain_sum, gain_sum_err).
        mjd_start and mjd_end are left open (NaN). Combine measurements at
        several levels with :func:`summarize_gain_rdnoise`.
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
        m = measure_ptc(f1, f2, b1, b2, sigma=sigma, cell=cell,
                        max_lag=max_lag)
        source = 'PTC: {0}, {1}; {2}, {3}'.format(
            *(os.path.basename(f) for f in (flat1, flat2, bias1, bias2)))
        rows.append((channel or '', chip, m['gain'], m['rdnoise'], np.nan,
                     np.nan, source) + tuple(m[c] for c in PTC_COLUMNS))
    table = Table(rows=rows, names=TABLE_COLUMNS + PTC_COLUMNS)
    table['gain'].unit = u.electron / u.adu
    table['rdnoise'].unit = u.electron
    table['level'].unit = u.adu
    table['rdnoise_adu'].unit = u.adu
    for col in ('gain_nn', 'gain_sum', 'gain_sum_err'):
        table[col].unit = u.electron / u.adu
    return table


def fit_gain_vs_level(levels, gains, min_sets=3):
    """Straight-line fit gain = gain0 + slope * level.

    The apparent photon-transfer gain rises with signal level when the
    brighter-fatter effect (or non-linearity) suppresses the variance; the
    zero-level intercept is the conversion gain. With fewer than
    ``min_sets`` points (or no spread in level) the median is returned
    instead (``model='median'``).

    Returns
    -------
    dict
        model, n, gain0, gain0_err, slope [e-/ADU per ADU], slope_err,
        rms (residual scatter), level_min, level_max. Uncertainties come
        from the residual scatter (NaN with fewer than three points).
    """
    levels = np.asarray(levels, dtype=float)
    gains = np.asarray(gains, dtype=float)
    ok = np.isfinite(levels) & np.isfinite(gains)
    levels, gains = levels[ok], gains[ok]
    n = len(gains)
    out = {'n': n, 'level_min': float(levels.min()) if n else np.nan,
           'level_max': float(levels.max()) if n else np.nan}
    if n < max(min_sets, 2) or np.ptp(levels) == 0:
        out.update(model='median', gain0=float(np.median(gains)) if n
                   else np.nan, gain0_err=np.nan, slope=np.nan,
                   slope_err=np.nan,
                   rms=float(np.std(gains, ddof=1)) if n > 1 else np.nan)
        return out
    design = np.vstack([np.ones(n), levels]).T
    coef, *_ = np.linalg.lstsq(design, gains, rcond=None)
    resid = gains - design @ coef
    if n > 2:
        s2 = np.sum(resid ** 2) / (n - 2)
        cov = s2 * np.linalg.inv(design.T @ design)
        errs = np.sqrt(np.diag(cov))
        rms = float(np.sqrt(s2))
    else:
        errs, rms = (np.nan, np.nan), np.nan
    out.update(model='linear', gain0=float(coef[0]), slope=float(coef[1]),
               gain0_err=float(errs[0]), slope_err=float(errs[1]), rms=rms)
    return out


FIT_COLUMNS = ['channel', 'chip', 'model', 'n', 'gain0', 'gain0_err',
               'slope_pct_per_10k', 'slope_pct_err', 'rms', 'level_min',
               'level_max', 'gain_median', 'rdnoise_adu', 'rdnoise_adu_std',
               'rdnoise', 'rho_slope_per_10k', 'gain_nn_slope_pct_per_10k',
               'rho_sum_slope_per_10k', 'gain_sum_median',
               'gain_sum_slope_pct_per_10k', 'gain_sum_slope_pct_err',
               'gain_flux', 'gain_flux_err', 'n_flux']


def summarize_gain_rdnoise(results, model='linear', min_sets=3,
                           source='', mjd_start=np.nan, mjd_end=np.nan,
                           flux_max_dt=60.0):
    """Combine per-set measurements into detector-table rows.

    Parameters
    ----------
    results : astropy.table.Table
        Output of :func:`measure_gain_rdnoise_files` for several sets
        (vstacked), including the ``level`` and ``rdnoise_adu`` columns.
    model : {'linear', 'median'}, optional
        'linear': gain = zero-level intercept of :func:`fit_gain_vs_level`
        (falls back to the median with fewer than ``min_sets`` sets).
        'median': median over sets (ignores any level dependence).
    source, mjd_start, mjd_end
        Written into the product rows.
    flux_max_dt : float or None, optional
        ``gain_flux`` uses only sets whose two flats were taken at most
        this many seconds apart (``flat_dt`` column of ``results``; all
        sets if that column is absent or ``flux_max_dt`` is None). In the
        2025-05 data, pairs 140-210 s apart gave gain_sum 0.5-0.8 % low,
        presumably because the twilight changes between the exposures and
        leaves small-scale structure that the covariance sum counts.
        Default: 60

    Returns
    -------
    rows : astropy.table.Table
        :data:`TABLE_COLUMNS` plus ``gain_flux`` (after ``rdnoise``), one
        row per channel/chip, with rdnoise = gain x median(rdnoise_adu) and
        gain_flux = median gain_sum of the selected sets (NaN if none).
    fit : astropy.table.Table
        :data:`FIT_COLUMNS`: the fit and diagnostics per channel/chip
        (slopes in percent of gain0 per 10,000 ADU; rho_slope is the change
        of rho_x + rho_y per 10,000 ADU; gain_nn_slope is the slope of the
        nearest-neighbour-corrected gain; rho_sum_slope and gain_sum_slope
        are the same for the covariances summed out to ``max_lag`` (see
        :func:`measure_ptc`); gain_sum_slope_pct_err is the larger of the
        fit-residual error and the error propagated from gain_sum_err.
        gain_sum_slope near zero means the trend is
        the brighter-fatter effect; a gain_sum_slope close to the gain
        slope means non-linearity). gain_flux_err is the scatter of the
        selected gain_sum values / sqrt(n_flux), or their propagated
        error if larger.
    """
    for col in ('level', 'rdnoise_adu'):
        if col not in results.colnames:
            raise ValueError("results need a '{0}' column: re-run "
                             "measure_gain_rdnoise_files".format(col))
    if model not in ('linear', 'median'):
        raise ValueError("model must be 'linear' or 'median'")

    rows, fit_rows = [], []
    channels = np.char.upper(np.asarray(results['channel'], dtype=str))
    for channel in sorted(set(channels)):
        for chip in sorted(set(results['chip'][channels == channel])):
            r = results[(channels == channel) & (results['chip'] == chip)]
            level = np.asarray(r['level'], dtype=float)
            gain = np.asarray(r['gain'], dtype=float)
            f = fit_gain_vs_level(level, gain,
                                  min_sets=min_sets if model == 'linear'
                                  else np.inf)
            rn_adu = np.asarray(r['rdnoise_adu'], dtype=float)
            g0 = f['gain0']
            rdnoise = g0 * float(np.median(rn_adu))

            def column(name):
                return (np.asarray(r[name], dtype=float)
                        if name in r.colnames else np.full(len(r), np.nan))

            def pct_slope(values):
                """Slope and error in % of the intercept per 10,000 ADU."""
                fv = fit_gain_vs_level(level, values, min_sets=3)
                if fv['model'] != 'linear' or not np.isfinite(fv['gain0']) \
                        or fv['gain0'] == 0:
                    return np.nan, np.nan
                return (1e6 * fv['slope'] / fv['gain0'],
                        1e6 * fv['slope_err'] / fv['gain0'])

            def slope_per_10k(values):
                fv = fit_gain_vs_level(level, values, min_sets=3)
                return 1e4 * fv['slope'] if fv['model'] == 'linear' \
                    else np.nan

            gsum = column('gain_sum')
            gsum_slope, gsum_slope_err = pct_slope(gsum)
            # With few sets the residual scatter can understate the error
            # of this noisy quantity; use at least the propagated error.
            gsum_err = column('gain_sum_err')
            ok = np.isfinite(gsum) & np.isfinite(gsum_err) & (gsum_err > 0)
            if np.isfinite(gsum_slope) and ok.sum() >= 2:
                w = 1.0 / gsum_err[ok] ** 2
                lw = np.sum(w * level[ok]) / np.sum(w)
                prop = 1.0 / np.sqrt(np.sum(w * (level[ok] - lw) ** 2))
                gsum0 = np.sum(w * gsum[ok]) / np.sum(w)
                gsum_slope_err = max(gsum_slope_err,
                                     1e6 * prop / gsum0)

            # Flux gain: gain_sum is level-independent, so take the median
            # over the sets whose flats are close in time
            sel = np.isfinite(gsum)
            if flux_max_dt is not None and 'flat_dt' in r.colnames:
                sel &= column('flat_dt') <= flux_max_dt
            n_flux = int(sel.sum())
            gain_flux = float(np.median(gsum[sel])) if n_flux else np.nan
            if n_flux > 1:
                gain_flux_err = max(
                    float(np.std(gsum[sel], ddof=1) / np.sqrt(n_flux)),
                    float(np.sqrt(np.nansum(gsum_err[sel] ** 2)) / n_flux))
            elif n_flux == 1:
                gain_flux_err = float(gsum_err[sel][0])
            else:
                gain_flux_err = np.nan

            rows.append((channel, int(chip), g0, rdnoise, gain_flux,
                         mjd_start, mjd_end, source))
            fit_rows.append((
                channel, int(chip), f['model'], f['n'], g0, f['gain0_err'],
                1e4 * f['slope'] / g0 * 100, 1e4 * f['slope_err'] / g0 * 100,
                f['rms'], f['level_min'], f['level_max'],
                float(np.median(gain)), float(np.median(rn_adu)),
                float(np.std(rn_adu, ddof=1)) if len(rn_adu) > 1 else np.nan,
                rdnoise,
                slope_per_10k(column('rho_x') + column('rho_y')),
                pct_slope(column('gain_nn'))[0],
                slope_per_10k(column('rho_sum')),
                float(np.nanmedian(gsum)) if np.isfinite(gsum).any()
                else np.nan,
                gsum_slope, gsum_slope_err, gain_flux, gain_flux_err,
                n_flux))

    names = TABLE_COLUMNS[:4] + ['gain_flux'] + TABLE_COLUMNS[4:]
    product = Table(rows=rows, names=names,
                    dtype=['U4', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'U200'])
    product['gain'].unit = u.electron / u.adu
    product['gain_flux'].unit = u.electron / u.adu
    product['rdnoise'].unit = u.electron
    fit = Table(rows=fit_rows, names=FIT_COLUMNS)
    for col in ('gain0', 'gain0_err', 'gain_median', 'gain_sum_median',
                'gain_flux', 'gain_flux_err'):
        fit[col].unit = u.electron / u.adu
    for col in ('level_min', 'level_max', 'rdnoise_adu', 'rdnoise_adu_std'):
        fit[col].unit = u.adu
    fit['rdnoise'].unit = u.electron
    return product, fit


def write_detector_table(table, filename=None, overwrite=False):
    """Write a detector table (default: the packaged ``lbc_detector.ecsv``)."""
    if filename is None:
        filename = DEFAULT_TABLE
    table.write(filename, format='ascii.ecsv', overwrite=overwrite)
    return filename
