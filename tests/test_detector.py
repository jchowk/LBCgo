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
