"""
Tests for LBCgo.detector: channel identification, gain/read-noise table
lookup with header fallback, and photon-transfer measurement.
"""

import numpy as np
import pytest
from pathlib import Path
from astropy.io import fits
from astropy.table import Table
from ccdproc import ImageFileCollection

from LBCgo import detector, masks
from conftest import (write_lbc_file, write_master_flat, LBC_KEYWORDS,
                      N_CHIPS)


# ---------------------------------------------------------------------------
# Channel identification
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('header, expected', [
    ({'INSTRUME': 'LBC_BLUE'}, 'LBCB'),       # real LBC-Blue spelling
    ({'INSTRUME': 'LBC-RED '}, 'LBCR'),       # real LBC-Red spelling
    ({'INSTRUME': 'LBC-BLUE'}, 'LBCB'),
    ({'DETECTOR': 'EEV-RED '}, 'LBCR'),
    ({}, None),
])
def test_lbc_channel_from_header(header, expected):
    assert detector.lbc_channel(header) == expected


def test_lbc_channel_from_filename():
    assert detector.lbc_channel({}, 'raw/lbcr.20141229.132703.fits') == 'LBCR'
    assert detector.lbc_channel({}, 'lbcb.20141120.065509_over.fits') == 'LBCB'
    assert detector.lbc_channel({}, 'other.fits') is None


# ---------------------------------------------------------------------------
# Table lookup and fallback
# ---------------------------------------------------------------------------

def _table():
    nan = np.nan
    return Table(rows=[
        ('LBCB', 2, 2.09, 10.0, nan, nan, 'open-ended'),
        ('LBCB', 2, 2.20, 11.0, 58000.0, nan, 'newer'),
        ('LBCR', 2, 2.14, 10.7, 56000.0, 57500.0, 'window'),
    ], names=detector.TABLE_COLUMNS)


def test_lookup_open_ended_row():
    assert detector.lookup_detector_params(_table(), 'LBCB', 2, 56981.3) == (2.09, 10.0)


def test_lookup_latest_start_wins():
    assert detector.lookup_detector_params(_table(), 'LBCB', 2, 58500.0) == (2.20, 11.0)


def test_lookup_date_window():
    t = _table()
    assert detector.lookup_detector_params(t, 'LBCR', 2, 57020.6) == (2.14, 10.7)
    assert detector.lookup_detector_params(t, 'LBCR', 2, 57500.0) is None
    assert detector.lookup_detector_params(t, 'LBCR', 2, 55000.0) is None


def test_lookup_no_match_or_missing_info():
    t = _table()
    assert detector.lookup_detector_params(t, 'LBCB', 1, 56981.3) is None
    assert detector.lookup_detector_params(t, None, 2, 56981.3) is None
    assert detector.lookup_detector_params(None, 'LBCB', 2, 56981.3) is None
    # Unknown date only matches fully open-ended rows
    assert detector.lookup_detector_params(t, 'LBCB', 2, None) == (2.09, 10.0)
    assert detector.lookup_detector_params(t, 'LBCR', 2, None) is None


def test_gain_rdnoise_prefers_table():
    primary = fits.Header({'INSTRUME': 'LBC_BLUE', 'MJD_OBS': 56981.3,
                           'GAIN': 1.75, 'RDNOISE': 12.0})
    gain, rn, source = detector.gain_rdnoise(2, [fits.Header(), primary],
                                             table=_table())
    assert (gain, rn, source) == (2.09, 10.0, 'table')


def test_gain_rdnoise_falls_back_to_header():
    chip_hdr = fits.Header({'GAIN': 1.75, 'RDNOISE': 12.0})
    primary = fits.Header({'INSTRUME': 'LBC-RED ', 'MJD_OBS': 56981.3})
    gain, rn, source = detector.gain_rdnoise(3, [chip_hdr, primary],
                                             table=_table())
    assert (gain, rn, source) == (1.75, 12.0, 'header')


def test_gain_rdnoise_default_when_nothing_available():
    gain, rn, source = detector.gain_rdnoise(1, [fits.Header()],
                                             table=_table())
    assert (gain, rn, source) == (masks.DEFAULT_GAIN, masks.DEFAULT_RDNOISE,
                                  'default')


# Giallongo et al. 2008, A&A 482, 349, Table 1 (LBC-Blue, 2006 commissioning)
GIALLONGO_LBCB = {1: (1.96, 11.4), 2: (2.09, 11.6), 3: (2.06, 11.6),
                  4: (1.98, 11.2)}


def test_packaged_table_seeded_with_lbcb_values():
    table = detector.read_detector_table()
    assert list(table.colnames) == detector.TABLE_COLUMNS
    assert len(table) == 4
    assert set(table['channel']) == {'LBCB'}
    for chip, (gain, rn) in GIALLONGO_LBCB.items():
        # Open-ended rows match any date, and an unknown date
        assert detector.lookup_detector_params(table, 'LBCB', chip,
                                               56981.3) == (gain, rn)
        assert detector.lookup_detector_params(table, 'LBCB', chip,
                                               None) == (gain, rn)


def test_packaged_table_lbcr_falls_back_to_header():
    hdr = fits.Header({'INSTRUME': 'LBC-RED ', 'MJD_OBS': 57020.56,
                       'GAIN': 1.75, 'RDNOISE': 12.0})
    assert detector.gain_rdnoise(2, [hdr]) == (1.75, 12.0, 'header')
    hdr_b = fits.Header({'INSTRUME': 'LBC_BLUE', 'MJD_OBS': 56981.3,
                         'GAIN': 1.75, 'RDNOISE': 12.0})
    assert detector.gain_rdnoise(2, [hdr_b]) == (2.09, 11.6, 'table')


def test_missing_explicit_table_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        detector.read_detector_table(str(tmp_path / 'nope.ecsv'))


def test_table_roundtrip(tmp_path):
    out = detector.write_detector_table(_table(), str(tmp_path / 't.ecsv'))
    back = detector.read_detector_table(out)
    assert detector.lookup_detector_params(back, 'LBCB', 2, 58500.0) == (2.20, 11.0)


# ---------------------------------------------------------------------------
# Photon-transfer measurement
# ---------------------------------------------------------------------------

def _simulate(rng, gain, rn_e, signal_e, prnu, shape, bias_adu=1000.0):
    """One frame in ADU: Poisson signal with fixed PRNU + read noise + bias."""
    electrons = rng.poisson(signal_e * prnu) if signal_e > 0 else 0.0
    electrons = electrons + rng.normal(0.0, rn_e, shape)
    return bias_adu + electrons / gain


@pytest.mark.parametrize('signal1, signal2', [
    (40000.0, 40000.0),     # equal flats (standard PTC pair)
    (40000.0, 44000.0),     # 10 % illumination difference
    (1500.0, 1650.0),       # faint: read noise is ~7 % of the variance
])
def test_measure_gain_rdnoise_arrays(signal1, signal2):
    rng = np.random.default_rng(3)
    shape = (500, 500)
    gain, rn = 2.05, 10.5
    prnu = 1.0 + 0.02 * rng.standard_normal(shape)   # fixed pattern, 2 %
    f1 = _simulate(rng, gain, rn, signal1, prnu, shape)
    f2 = _simulate(rng, gain, rn, signal2, prnu, shape)
    b1 = _simulate(rng, gain, rn, 0.0, prnu, shape)
    b2 = _simulate(rng, gain, rn, 0.0, prnu, shape)
    g, r = detector.measure_gain_rdnoise(f1, f2, b1, b2)
    assert abs(g / gain - 1) < 0.015
    assert abs(r / rn - 1) < 0.02


def _write_raw(path, frames, imagetyp, instrume='LBC-RED '):
    """Raw-format MEF: 20 px prescan + science + overscan, like LBC chips."""
    primary = fits.Header({'IMAGETYP': imagetyp, 'INSTRUME': instrume,
                           'FILTER': 'I-BESSEL', 'MJD_OBS': 57020.56,
                           'GAIN': 1.75, 'RDNOISE': 12.0})
    hdul = fits.HDUList([fits.PrimaryHDU(header=primary)])
    for chip, sci in enumerate(frames, start=1):
        ny, nx = sci.shape
        full = np.empty((ny, nx + 40), dtype='float32')
        full[:, :20] = sci[:, :20]                 # prescan (unused)
        full[:, 20:20 + nx] = sci
        full[:, 20 + nx:] = 1000.0 + np.random.default_rng(chip).normal(
            0, 5.0, (ny, 20))                     # overscan at bias level
        hdr = fits.Header({'EXTNAME': 'LBCCHIP{0}'.format(chip),
                           'BUNIT': 'adu', 'SATURATE': 65536,
                           'TRIMSEC': '[21:{0},1:{1}]'.format(20 + nx, ny),
                           'BIASSEC': '[{0}:{1},1:{2}]'.format(21 + nx, 40 + nx, ny)})
        hdul.append(fits.ImageHDU(full, header=hdr))
    hdul.writeto(path)
    return str(path)


def test_measure_gain_rdnoise_files(tmp_path):
    rng = np.random.default_rng(4)
    shape = (300, 300)
    gains = [2.08, 2.14, 2.13, 2.09]
    rns = [10.4, 10.7, 11.3, 10.0]
    prnus = [1.0 + 0.01 * rng.standard_normal(shape) for _ in gains]

    def frames(signal):
        return [_simulate(rng, g, r, signal, p, shape)
                for g, r, p in zip(gains, rns, prnus)]

    f1 = _write_raw(tmp_path / 'lbcr.f1.fits', frames(30000.0), 'flat')
    f2 = _write_raw(tmp_path / 'lbcr.f2.fits', frames(30500.0), 'flat')
    b1 = _write_raw(tmp_path / 'lbcr.b1.fits', frames(0.0), 'zero')
    b2 = _write_raw(tmp_path / 'lbcr.b2.fits', frames(0.0), 'zero')

    table = detector.measure_gain_rdnoise_files(f1, f2, b1, b2, box=None)
    assert list(table['chip']) == [1, 2, 3, 4]
    assert set(table['channel']) == {'LBCR'}
    np.testing.assert_allclose(table['gain'], gains, rtol=0.03)
    np.testing.assert_allclose(table['rdnoise'], rns, rtol=0.05)
    assert np.all(np.isnan(table['mjd_start']))

    # A measured table feeds straight back into the lookup
    out = detector.write_detector_table(table, str(tmp_path / 'red.ecsv'))
    found = detector.lookup_detector_params(detector.read_detector_table(out),
                                            'LBCR', 3, 57020.56)
    assert np.isclose(found[0], table['gain'][2])


def test_measure_rejects_saturated_flats(tmp_path):
    rng = np.random.default_rng(5)
    shape = (100, 100)
    prnu = np.ones(shape)
    hot = [_simulate(rng, 1.0, 5.0, 60000.0, prnu, shape) for _ in range(4)]
    cold = [_simulate(rng, 1.0, 5.0, 0.0, prnu, shape) for _ in range(4)]
    f1 = _write_raw(tmp_path / 'a.fits', hot, 'flat')
    f2 = _write_raw(tmp_path / 'b.fits', hot, 'flat')
    b1 = _write_raw(tmp_path / 'c.fits', cold, 'zero')
    b2 = _write_raw(tmp_path / 'd.fits', cold, 'zero')
    with pytest.raises(ValueError, match='saturation'):
        detector.measure_gain_rdnoise_files(f1, f2, b1, b2, box=None)


# ---------------------------------------------------------------------------
# go_flatfield integration
# ---------------------------------------------------------------------------

def _flatfield_with_table(raw_dir, work_dir, table_path):
    from LBCgo.lbcproc import go_overscan, go_flatfield
    path = write_lbc_file(raw_dir, 'lbcb.20230101.000001.fits')
    with fits.open(path, mode='update') as hdul:
        hdul[0].header['INSTRUME'] = 'LBC_BLUE'
        hdul[0].header['MJD_OBS'] = 56981.3
        hdul[0].header['GAIN'] = 1.75
        hdul[0].header['RDNOISE'] = 12.0
    ic = ImageFileCollection(str(raw_dir), keywords=LBC_KEYWORDS)
    over = go_overscan(ic, image_directory=str(work_dir) + '/',
                       raw_directory=str(raw_dir) + '/', verbose=False)
    over_dir = str(work_dir) + '/'
    write_master_flat(over_dir, filter_name='g-SLOAN')
    ic2 = ImageFileCollection(over_dir, keywords=LBC_KEYWORDS,
                              filenames=[Path(f).name for f in over])
    flat_files = go_flatfield(ic2, flat_directory=over_dir,
                              image_directory=over_dir,
                              input_directory=over_dir, cosmiccorrect=False,
                              verbose=False, return_files=True,
                              detector_table=table_path)
    return fits.open(masks.sidecar_name(over_dir + flat_files[0], 'weight'))


def test_flatfield_uses_table_then_header(raw_dir, work_dir, tmp_path):
    table_path = detector.write_detector_table(_table(),
                                               str(tmp_path / 'det.ecsv'))
    with _flatfield_with_table(raw_dir, work_dir, table_path) as wh:
        # Chip 2 has a LBCB row valid at MJD 56981.3
        assert wh[2].header['GAINSRC'] == 'table'
        assert wh[2].header['GAIN'] == 2.09
        assert wh[2].header['RDNOISE'] == 10.0
        # Other chips fall back to the primary-header values
        for ext in (1, 3, 4):
            assert wh[ext].header['GAINSRC'] == 'header'
            assert wh[ext].header['GAIN'] == 1.75
            assert wh[ext].header['RDNOISE'] == 12.0


def test_flatfield_default_table_uses_seeded_lbcb(raw_dir, work_dir):
    """LBCB data with the packaged table: every chip gets Table 1 values."""
    with _flatfield_with_table(raw_dir, work_dir, None) as wh:
        for ext in range(1, N_CHIPS + 1):
            gain, rn = GIALLONGO_LBCB[ext]
            assert wh[ext].header['GAINSRC'] == 'table'
            assert wh[ext].header['GAIN'] == gain
            assert wh[ext].header['RDNOISE'] == rn


def test_flatfield_empty_table_uses_headers(raw_dir, work_dir, tmp_path):
    empty = detector.write_detector_table(_table()[:0],
                                          str(tmp_path / 'empty.ecsv'))
    with _flatfield_with_table(raw_dir, work_dir, empty) as wh:
        assert all(wh[ext].header['GAINSRC'] == 'header'
                   for ext in range(1, N_CHIPS + 1))


# ---------------------------------------------------------------------------
# Robustness to illumination mismatch between the two flats (cell variance)
# ---------------------------------------------------------------------------

def _mismatched_pair(rng, mismatch, shape=(600, 600), gain=2.05, rn=11.5,
                     level_adu=34000.0):
    """Two flats whose illumination differs by a linear gradient of
    `mismatch` (fractional, peak to peak across the region), plus biases."""
    ny, nx = shape
    xx = np.mgrid[0:ny, 0:nx][1] / nx
    prnu = 1.0 + 0.01 * rng.standard_normal(shape)
    signal = level_adu * gain
    f1 = (rng.poisson(signal * prnu) + rng.normal(0, rn, shape)) / gain
    f2 = (rng.poisson(0.96 * signal * prnu * (1 + mismatch * (xx - 0.5)))
          + rng.normal(0, rn, shape)) / gain
    b1 = rng.normal(0, rn, shape) / gain
    b2 = rng.normal(0, rn, shape) / gain
    return f1, f2, b1, b2


def test_gain_robust_to_illumination_mismatch():
    """A 1 % peak-to-peak gradient between flats (as in real twilight pairs)
    biases a whole-region variance strongly low; 50-px cells do not."""
    rng = np.random.default_rng(21)
    frames = _mismatched_pair(rng, 0.01)
    g_whole, rn_whole = detector.measure_gain_rdnoise(*frames, cell=None)
    g_cell, rn_cell = detector.measure_gain_rdnoise(*frames)       # cell=50
    assert g_whole / 2.05 - 1 < -0.05          # the bias seen in real data
    assert abs(g_cell / 2.05 - 1) < 0.015
    assert abs(rn_cell / 11.5 - 1) < 0.02


def test_gain_cells_match_whole_region_without_mismatch():
    rng = np.random.default_rng(22)
    frames = _mismatched_pair(rng, 0.0)
    g_whole, _ = detector.measure_gain_rdnoise(*frames, cell=None)
    g_cell, _ = detector.measure_gain_rdnoise(*frames)
    assert abs(g_cell / g_whole - 1) < 0.01
    assert abs(g_cell / 2.05 - 1) < 0.015


def test_cell_variance_falls_back_to_whole_region():
    rng = np.random.default_rng(23)
    d = rng.normal(0, 3.0, (40, 60))
    assert detector._cell_variance(d, 50, 4.0) == \
        detector._cell_variance(d, None, 4.0)
    assert abs(detector._cell_variance(rng.normal(0, 3.0, (500, 500)),
                                       50, 4.0) / 9.0 - 1) < 0.02


# ---------------------------------------------------------------------------
# Ambiguous (tied) detector-table rows
# ---------------------------------------------------------------------------

def _two_rows(start_new):
    nan = np.nan
    return Table(rows=[('LBCB', 1, 1.96, 11.4, nan, nan, 'old'),
                       ('LBCB', 1, 1.84, 8.9, start_new, nan, 'new')],
                 names=detector.TABLE_COLUMNS)


def test_lookup_warns_on_tied_open_rows():
    t = _two_rows(np.nan)
    with pytest.warns(UserWarning, match='same mjd_start'):
        assert detector.lookup_detector_params(t, 'LBCB', 1, 60822.3) == \
            (1.96, 11.4)
    assert detector.detector_table_conflicts(t) == [('LBCB', 1, None)]


def test_lookup_newer_finite_start_supersedes_without_warning(recwarn):
    t = _two_rows(60822.0)
    assert detector.lookup_detector_params(t, 'LBCB', 1, 60822.3) == (1.84, 8.9)
    assert detector.lookup_detector_params(t, 'LBCB', 1, 56981.3) == (1.96, 11.4)
    assert not [w for w in recwarn if issubclass(w.category, UserWarning)
                and 'mjd_start' in str(w.message)]
    assert detector.detector_table_conflicts(t) == []


def test_packaged_table_has_no_conflicts():
    assert detector.detector_table_conflicts(
        detector.read_detector_table()) == []


# ---------------------------------------------------------------------------
# Level dependence of the gain: diagnostics, fit and product rows
# ---------------------------------------------------------------------------

def test_measure_ptc_diagnostics():
    rng = np.random.default_rng(31)
    shape = (500, 500)
    gain, rn = 2.05, 10.5
    prnu = 1.0 + 0.02 * rng.standard_normal(shape)
    f1 = _simulate(rng, gain, rn, 40000.0, prnu, shape)
    f2 = _simulate(rng, gain, rn, 44000.0, prnu, shape)
    b1 = _simulate(rng, gain, rn, 0.0, prnu, shape)
    b2 = _simulate(rng, gain, rn, 0.0, prnu, shape)
    m = detector.measure_ptc(f1, f2, b1, b2)
    assert (m['gain'], m['rdnoise']) == detector.measure_gain_rdnoise(
        f1, f2, b1, b2)
    assert abs(m['level'] / (42000.0 / gain) - 1) < 0.005
    assert abs(m['k'] / (40000.0 / 44000.0) - 1) < 0.005
    assert abs(m['rdnoise_adu'] / (rn / gain) - 1) < 0.02
    assert np.isclose(m['rdnoise'], m['gain'] * m['rdnoise_adu'])
    # Independent pixels: no correlations, gain_nn == gain within noise
    for key in ('rho_x', 'rho_y', 'bias_rho_x', 'bias_rho_y'):
        assert abs(m[key]) < 0.01
    assert abs(m['gain_nn'] / m['gain'] - 1) < 0.03


def test_cell_correlations_known_values():
    rng = np.random.default_rng(32)
    x = rng.normal(0, 5.0, (600, 600))
    rx, ry = detector._cell_correlations(x, 50, 4.0)
    assert abs(rx) < 0.01 and abs(ry) < 0.01
    # d = x + x shifted by one column: rho_x = 1/2, rho_y = 0
    d = x + np.roll(x, 1, axis=1)
    rx, ry = detector._cell_correlations(d, 50, 4.0)
    assert abs(rx - 0.5) < 0.02 and abs(ry) < 0.01
    rx, ry = detector._cell_correlations(d.T, None, 4.0)
    assert abs(ry - 0.5) < 0.02 and abs(rx) < 0.01


def test_fit_gain_vs_level_linear():
    levels = np.array([5000., 10000., 20000., 30000.])
    gains = 1.75 * (1 + 0.03 * levels / 1e4)
    f = detector.fit_gain_vs_level(levels, gains)
    assert f['model'] == 'linear' and f['n'] == 4
    assert np.isclose(f['gain0'], 1.75)
    assert np.isclose(f['slope'], 1.75 * 0.03 / 1e4)
    assert f['rms'] < 1e-10 and f['gain0_err'] < 1e-10
    assert (f['level_min'], f['level_max']) == (5000., 30000.)
    # Noisy points: the uncertainties come from the residuals
    noisy = gains + np.array([0.002, -0.002, -0.002, 0.002])
    f = detector.fit_gain_vs_level(levels, noisy)
    assert f['rms'] > 0 and f['gain0_err'] > 0 and f['slope_err'] > 0
    # NaNs are dropped
    f = detector.fit_gain_vs_level(np.append(levels, np.nan),
                                   np.append(gains, 2.0))
    assert f['n'] == 4 and np.isclose(f['gain0'], 1.75)


@pytest.mark.parametrize('levels, gains', [
    ([10000., 20000.], [1.80, 1.85]),               # too few sets
    ([20000., 20000., 20000.], [1.80, 1.85, 1.81]),  # no spread in level
    ([20000.], [1.80]),
])
def test_fit_gain_vs_level_falls_back_to_median(levels, gains):
    f = detector.fit_gain_vs_level(levels, gains)
    assert f['model'] == 'median'
    assert np.isclose(f['gain0'], np.median(gains))
    assert np.isnan(f['slope'])


def _results(levels, gains, rn_adu, channel='LBCR', chip=1):
    n = len(levels)
    return Table({'channel': [channel] * n, 'chip': [chip] * n,
                  'gain': gains, 'rdnoise': np.multiply(gains, rn_adu),
                  'level': levels, 'rdnoise_adu': rn_adu,
                  'rho_x': np.zeros(n), 'rho_y': np.zeros(n),
                  'gain_nn': gains})


def test_summarize_gain_rdnoise_linear_and_median():
    from astropy.table import vstack
    levels = np.array([5000., 10000., 20000., 30000.])
    r = vstack([_results(levels, 1.75 * (1 + 0.03 * levels / 1e4),
                         [5.0, 5.2, 5.1, 5.3]),
                _results(levels[:2], [2.0, 2.1], [6.0, 6.0], chip=2)])
    rows, fit = detector.summarize_gain_rdnoise(r, source='test',
                                                mjd_start=60000.)
    assert list(rows.colnames) == detector.TABLE_COLUMNS
    assert list(fit.colnames) == detector.FIT_COLUMNS
    assert str(rows['gain'].unit) == 'electron / adu'
    assert str(rows['rdnoise'].unit) == 'electron'
    assert rows['source'][0] == 'test' and rows['mjd_start'][0] == 60000.
    # chip 1: intercept, not the median over sets; RN = g0 x median RN_ADU
    assert np.isclose(rows['gain'][0], 1.75)
    assert np.isclose(rows['rdnoise'][0], 1.75 * 5.15)
    assert fit['model'][0] == 'linear'
    assert np.isclose(fit['slope_pct_per_10k'][0], 3.0)
    assert np.isclose(fit['gain_median'][0], np.median(r['gain'][:4]))
    assert np.isclose(fit['gain_nn_slope_pct_per_10k'][0], 3.0)
    # chip 2: two sets -> median fallback
    assert fit['model'][1] == 'median'
    assert np.isclose(rows['gain'][1], 2.05)
    assert np.isclose(rows['rdnoise'][1], 2.05 * 6.0)
    # model='median' ignores the level dependence
    rows_m, fit_m = detector.summarize_gain_rdnoise(r, model='median')
    assert np.isclose(rows_m['gain'][0], np.median(r['gain'][:4]))
    assert set(fit_m['model']) == {'median'}


def test_summarize_gain_rdnoise_rejects_old_results():
    r = _results([1e4, 2e4, 3e4], [1.8, 1.8, 1.8], [5., 5., 5.])
    with pytest.raises(ValueError, match='model'):
        detector.summarize_gain_rdnoise(r, model='quadratic')
    r.remove_column('level')
    with pytest.raises(ValueError, match='level'):
        detector.summarize_gain_rdnoise(r)


def _bf_flat(rng, gain, signal_e, prnu, rn_e, a_per_e, bias_adu=1000.0):
    """Flat with a toy brighter-fatter effect: each pixel keeps 1 - 4a of
    its charge and gives a to each of its four neighbours (charge is
    conserved), with a proportional to the signal."""
    e = rng.poisson(signal_e * prnu).astype(float)
    a = a_per_e * signal_e
    e = (1 - 4 * a) * e + a * sum(np.roll(e, s, axis=ax)
                                  for s in (1, -1) for ax in (0, 1))
    return bias_adu + (e + rng.normal(0.0, rn_e, e.shape)) / gain


def test_brighter_fatter_level_dependence_is_recovered():
    """Smoothing by the brighter-fatter effect lowers the per-pixel
    variance, so the apparent gain rises with level (~+4 % per 10k ADU
    here, as in LBCR); the zero-level intercept recovers the true gain,
    the neighbour correlations grow with level, and gain_nn is flat."""
    rng = np.random.default_rng(33)
    shape = (500, 500)
    gain, rn = 1.75, 10.0
    a_per_e = 0.005 / (1e4 * gain)       # a = 0.005 at 10,000 ADU
    prnu = 1.0 + 0.01 * rng.standard_normal(shape)
    rows = []
    for adu in (5000., 10000., 20000., 30000.):
        s = adu * gain
        f1 = _bf_flat(rng, gain, s, prnu, rn, a_per_e)
        f2 = _bf_flat(rng, gain, 1.03 * s, prnu, rn, a_per_e)
        b1 = _simulate(rng, gain, rn, 0.0, prnu, shape)
        b2 = _simulate(rng, gain, rn, 0.0, prnu, shape)
        rows.append(detector.measure_ptc(f1, f2, b1, b2))
    t = Table(rows=[[r[c] for c in ('gain', 'level', 'rho_x', 'rho_y',
                                    'gain_nn', 'rdnoise_adu')]
                    for r in rows],
              names=['gain', 'level', 'rho_x', 'rho_y', 'gain_nn',
                     'rdnoise_adu'])
    # Apparent gain rises with level; the median over sets is biased high
    assert np.all(np.diff(t['gain']) > 0)
    assert np.median(t['gain']) / gain - 1 > 0.04
    # Correlations positive and growing (rho ~ 2a)
    assert np.all(t['rho_x'] > 0) and np.all(np.diff(t['rho_x']) > 0)
    assert abs(t['rho_x'][-1] / (2 * 0.015) - 1) < 0.15
    # Intercept recovers the true gain; gain_nn has (almost) no trend
    t['channel'], t['chip'] = 'LBCR', 1
    product, fit = detector.summarize_gain_rdnoise(t)
    assert abs(product['gain'][0] / gain - 1) < 0.015
    assert abs(product['rdnoise'][0] / rn - 1) < 0.03
    assert 3.0 < fit['slope_pct_per_10k'][0] < 5.5
    assert fit['rho_slope_per_10k'][0] > 0
    assert abs(fit['gain_nn_slope_pct_per_10k'][0]) < 0.5


# ---------------------------------------------------------------------------
# Covariances summed over lags: brighter-fatter vs non-linearity
# ---------------------------------------------------------------------------

def test_covariance_sum_white_noise_and_gradient():
    rng = np.random.default_rng(41)
    x = rng.normal(0, 5.0, (1000, 1000))
    S, S_err, rho = detector._cell_covariance_sum(x, 50, 4.0, 3)
    assert rho.shape == (7, 7) and rho[3, 3] == 1.0
    assert 0.005 < S_err < 0.02
    # Without the -p(1+S)/n correction S would be about -0.06
    assert abs(S) < 3.5 * S_err
    # A gradient (here 1.25 ADU across a block) is removed with the plane;
    # removing only the mean would add ~0.25 to S
    yy, xx = np.mgrid[:1000, :1000]
    S_ramp, _, _ = detector._cell_covariance_sum(x + 0.025 * xx + 0.01 * yy,
                                                 50, 4.0, 3)
    assert abs(S_ramp - S) < 0.002


def test_covariance_sum_known_correlation():
    rng = np.random.default_rng(42)
    x = rng.normal(0, 5.0, (1000, 1000))
    # d = x + x shifted by two columns: rho(dx=+-2) = 1/2, S = 1
    S, S_err, rho = detector._cell_covariance_sum(
        x + np.roll(x, 2, axis=1), 50, 4.0, 3)
    assert abs(S - 1.0) < 3.5 * S_err
    assert abs(rho[3, 5] - 0.5) < 0.01 and abs(rho[3, 1] - 0.5) < 0.01
    assert abs(rho[3, 4]) < 0.01 and abs(rho[5, 3]) < 0.01


def test_measure_ptc_max_lag_zero_skips_sum():
    rng = np.random.default_rng(43)
    shape = (300, 300)
    prnu = np.ones(shape)
    frames = [_simulate(rng, 2.0, 10.0, s, prnu, shape)
              for s in (20000.0, 20000.0, 0.0, 0.0)]
    m0 = detector.measure_ptc(*frames, max_lag=0)
    m3 = detector.measure_ptc(*frames)
    assert np.isnan(m0['gain_sum']) and np.isnan(m0['rho_sum'])
    assert m0['gain'] == m3['gain'] and np.isfinite(m3['gain_sum'])
    assert set(detector.PTC_COLUMNS) <= set(m3)


def _flat_bf_nl(rng, gain, signal_e, prnu, rn_e, a_per_e, b_per_e, beta):
    """Flat with a brighter-fatter effect reaching two pixels (each pixel
    gives a to its four nearest neighbours and b to the four at distance 2,
    both proportional to the signal; charge is conserved) and a sublinear
    response N (1 - beta N)."""
    e = rng.poisson(signal_e * prnu).astype(float)
    a, b = a_per_e * signal_e, b_per_e * signal_e

    def ring(k):
        return sum(np.roll(e, s, axis=ax) for s in (k, -k) for ax in (0, 1))

    e = (1 - 4 * a - 4 * b) * e + a * ring(1) + b * ring(2)
    e = e * (1 - beta * e)
    return 1000.0 + (e + rng.normal(0.0, rn_e, e.shape)) / gain


def _ptc_series(seed, a_per_e=0.0, b_per_e=0.0, beta=0.0):
    rng = np.random.default_rng(seed)
    shape = (1000, 1000)
    gain, rn = 1.75, 10.0
    prnu = 1.0 + 0.01 * rng.standard_normal(shape)
    rows = []
    for adu in (5000., 10000., 20000., 30000.):
        s = adu * gain
        m = detector.measure_ptc(
            _flat_bf_nl(rng, gain, s, prnu, rn, a_per_e, b_per_e, beta),
            _flat_bf_nl(rng, gain, 1.03 * s, prnu, rn, a_per_e, b_per_e,
                        beta),
            _simulate(rng, gain, rn, 0.0, prnu, shape),
            _simulate(rng, gain, rn, 0.0, prnu, shape))
        m.update(channel='LBCR', chip=1)
        rows.append(m)
    return detector.summarize_gain_rdnoise(Table(rows=rows))[1][0]


def test_gain_sum_removes_long_range_brighter_fatter():
    """Half of the charge sharing goes two pixels away: gain_nn keeps part
    of the slope, gain_sum (lags <= 3) removes it."""
    a = 0.003 / (1e4 * 1.75)             # a = b = 0.003 at 10,000 ADU
    f = _ptc_series(44, a_per_e=a, b_per_e=a)
    assert f['slope_pct_per_10k'] > 4.0
    assert f['gain_nn_slope_pct_per_10k'] > 1.5
    assert abs(f['gain_sum_slope_pct_per_10k']) < \
        3 * f['gain_sum_slope_pct_err']
    assert f['rho_sum_slope_per_10k'] > 0.03
    assert abs(f['gain0'] / 1.75 - 1) < 0.015
    assert abs(f['gain_sum_median'] / 1.75 - 1) < 0.015


def test_gain_sum_keeps_nonlinearity():
    """A sublinear response raises the apparent gain without creating
    covariances: gain_sum has the same slope as gain."""
    f = _ptc_series(45, beta=7.6e-7)     # ~ +4.8 % per 10,000 ADU
    assert f['slope_pct_per_10k'] > 4.0
    assert abs(f['gain_sum_slope_pct_per_10k'] - f['slope_pct_per_10k']) \
        < 3 * f['gain_sum_slope_pct_err']
    assert abs(f['rho_sum_slope_per_10k']) < 0.01
    assert abs(f['gain0'] / 1.75 - 1) < 0.015


def test_summarize_without_gain_sum_columns():
    r = _results([1e4, 2e4, 3e4], [1.8, 1.82, 1.84], [5., 5., 5.])
    _, fit = detector.summarize_gain_rdnoise(r)
    assert np.isnan(fit['gain_sum_slope_pct_per_10k'][0])
    assert np.isnan(fit['gain_sum_median'][0])
    assert list(fit.colnames) == detector.FIT_COLUMNS
