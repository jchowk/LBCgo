"""
Tests for LBCgo/headers.py and for lbcproc on raw frames with an invalid
header card (astromatic plan Stage A1).

Some raw LBC frames carry an unquoted ``PA_PNT = nan``, which astropy
cannot parse. ccdproc's ImageFileCollection parses every card, so
one such frame stopped lbcgo before any processing. HeaderCollection reads
only the requested keywords; clean_header drops unparsable cards before a
raw header is copied into an output file.
"""

import numpy as np
import pytest
from astropy.io import fits
from astropy.io.fits.verify import VerifyError
from ccdproc import ImageFileCollection

from conftest import (LBC_KEYWORDS, NX_SCIENCE, add_unparsable_card,
                      write_lbc_file)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def mixed_dir(raw_dir):
    """Raw frames of two filters and two types, one with no PROPID."""
    write_lbc_file(raw_dir, 'lbcb.20230101.000001.fits', imagetyp='object',
                   filter_name='g-SLOAN', object_name='NGC891')
    write_lbc_file(raw_dir, 'lbcb.20230101.000002.fits', imagetyp='object',
                   filter_name='r-SLOAN', object_name='J1419+4207')
    write_lbc_file(raw_dir, 'lbcb.20230101.000003.fits', imagetyp='flat',
                   filter_name='g-SLOAN', object_name='SkyFlat',
                   lbcobnam='SkyFlat', pixel_value=20000.0)
    p = write_lbc_file(raw_dir, 'lbcb.20230101.000004.fits',
                       imagetyp='flat', filter_name='r-SLOAN',
                       object_name='SkyFlat', lbcobnam='SkyFlatTest',
                       pixel_value=20000.0)
    with fits.open(p, mode='update') as hdul:
        del hdul[0].header['PROPID']
    return raw_dir


# ---------------------------------------------------------------------------
# The failure being fixed
# ---------------------------------------------------------------------------

def test_unparsable_card_breaks_imagefilecollection(raw_dir):
    """Guard: the synthetic bad card does reproduce the real failure."""
    p = write_lbc_file(raw_dir, 'lbcb.20100318.223618.fits')
    add_unparsable_card(p)
    with pytest.raises(VerifyError, match='PA_PNT'):
        ImageFileCollection(str(raw_dir), keywords=LBC_KEYWORDS)
    with pytest.raises(VerifyError):
        fits.PrimaryHDU(header=fits.getheader(p)).writeto(
            raw_dir / 'copy.fits')


# ---------------------------------------------------------------------------
# HeaderCollection
# ---------------------------------------------------------------------------

def test_collection_scans_frame_with_unparsable_card(raw_dir):
    from LBCgo.headers import HeaderCollection
    p = write_lbc_file(raw_dir, 'lbcb.20100318.223618.fits',
                       filter_name='V-BESSEL')
    add_unparsable_card(p)
    ic = HeaderCollection(str(raw_dir), keywords=LBC_KEYWORDS + ['PA_PNT'])
    assert ic.files == ['lbcb.20100318.223618.fits']
    assert ic.summary['filter'][0] == 'V-BESSEL'
    # The bad keyword itself reads as missing
    assert ic.values('pa_pnt') == [None]
    assert list(ic.files_filtered(pa_pnt='*')) == []


def test_collection_matches_imagefilecollection(mixed_dir):
    """Same answers as ccdproc on well-formed headers, for the calls
    lbcproc makes."""
    from LBCgo.headers import HeaderCollection
    ours = HeaderCollection(str(mixed_dir), keywords=LBC_KEYWORDS)
    ref = ImageFileCollection(str(mixed_dir), keywords=LBC_KEYWORDS)

    assert ours.files == sorted(ref.files)
    queries = [dict(imagetyp='object'), dict(imagetyp='OBJECT'),
               dict(filter='g-SLOAN'), dict(object='J1419+4207'),
               dict(imagetyp='flat', object='SkyFlat', filter='g-SLOAN'),
               dict(propid='*'), dict(filter='i-SLOAN')]
    for q in queries:
        assert sorted(ours.files_filtered(**q)) == \
            sorted(ref.files_filtered(**q)), q
    for key in ['filter', 'object', 'exptime']:
        assert ours.values(key, unique=True) == \
            sorted(set(ref.values(key, unique=True))), key
    # Column-wise summary agrees, file by file
    ref_rows = {r['file']: r for r in ref.summary}
    for row in ours.summary:
        for key in ['object', 'filter', 'imagetyp', 'lbcobnam', 'exptime']:
            assert row[key] == ref_rows[row['file']][key]


def test_collection_missing_keyword_is_masked(mixed_dir):
    from LBCgo.headers import HeaderCollection
    ic = HeaderCollection(str(mixed_dir), keywords=LBC_KEYWORDS)
    missing = list(np.ma.getmaskarray(ic.summary['propid']))
    assert missing == [False, False, False, True]
    assert ic.values('propid')[3] is None
    assert ic.values('propid', unique=True) == ['TEST']
    assert len(ic.files_filtered(propid='*')) == 3


def test_collection_filenames_and_glob(mixed_dir):
    from LBCgo.headers import HeaderCollection
    ic = HeaderCollection(str(mixed_dir), keywords=['filter'],
                          filenames=['lbcb.20230101.000002.fits'])
    assert ic.files == ['lbcb.20230101.000002.fits']
    assert ic.values('filter') == ['r-SLOAN']
    write_lbc_file(mixed_dir, 'lbcr.20230101.000009.fits')
    ic = HeaderCollection(str(mixed_dir), keywords=['filter'],
                          glob_include='lbcr.*.*.fits*')
    assert ic.files == ['lbcr.20230101.000009.fits']


# ---------------------------------------------------------------------------
# clean_header
# ---------------------------------------------------------------------------

def test_clean_header_drops_only_bad_card(raw_dir):
    from LBCgo.headers import clean_header
    p = add_unparsable_card(write_lbc_file(raw_dir, 'a.fits'))
    raw = fits.getheader(p)
    clean = clean_header(raw)

    assert 'PA_PNT' not in clean
    kept = [k for k in raw.keys() if k != 'PA_PNT']
    assert [k for k in clean.keys() if k != 'HISTORY'] == kept
    assert any('PA_PNT' in h for h in clean['HISTORY'])
    assert clean['OBJECT'] == 'NGC891'
    # The cleaned header can be written
    fits.PrimaryHDU(header=clean).writeto(raw_dir / 'copy.fits')


def test_clean_header_good_header_unchanged(raw_dir):
    from LBCgo.headers import clean_header
    p = write_lbc_file(raw_dir, 'a.fits')
    raw = fits.getheader(p)
    clean = clean_header(raw)
    assert list(clean.items()) == list(raw.items())
    assert 'HISTORY' not in clean


# ---------------------------------------------------------------------------
# lbcproc on frames with the bad card
# ---------------------------------------------------------------------------

def _bad_raw_dir(raw_dir):
    """One object and three flats, every one with PA_PNT = nan in the
    primary header and in the chip 1 header."""
    paths = [write_lbc_file(raw_dir, 'lbcb.20100318.000001.fits',
                            imagetyp='object', object_name='J1419+4207',
                            lbcobnam='J1419+4207')]
    for i in range(2, 5):
        paths.append(write_lbc_file(
            raw_dir, f'lbcb.20100318.00000{i}.fits', imagetyp='flat',
            object_name='SkyFlat', lbcobnam='SkyFlat',
            pixel_value=20000.0 * (1 + 0.1 * i)))
    for p in paths:
        add_unparsable_card(p, extensions=(0, 1))
    return raw_dir


def test_lbcgo_runs_on_frames_with_unparsable_card(raw_dir, work_dir):
    """End to end (no astrometry), including make_flatfield, which now runs
    on the small synthetic chips (A4)."""
    from LBCgo import lbcproc
    _bad_raw_dir(raw_dir)
    out = work_dir / 'reduced'
    out.mkdir()
    lbcproc.lbcgo(raw_directory=str(raw_dir), image_directory=str(out),
                  do_astrometry=False, clean=False, verbose=False)

    assert (out / 'flat.g-SLOAN.fits').is_file()
    filter_dir = out / 'J1419+4207' / 'g-SLOAN'
    for chip in range(1, 5):
        assert (filter_dir / f'lbcb.20100318.000001_{chip}.fits').is_file()

    hdr = fits.getheader(out / 'data' / 'lbcb.20100318.000001_over.fits')
    assert 'PA_PNT' not in hdr
    assert any('PA_PNT' in h for h in hdr['HISTORY'])
    flat = fits.getheader(out / 'flat.g-SLOAN.fits')
    assert 'PA_PNT' not in flat


def test_make_bias_on_frames_with_unparsable_card(raw_dir, work_dir):
    from LBCgo.headers import HeaderCollection
    from LBCgo.lbcproc import make_bias
    for i in range(1, 4):
        p = write_lbc_file(raw_dir, f'lbcb.20100318.00010{i}.fits',
                           imagetyp='zero', object_name='Bias',
                           lbcobnam='Bias', pixel_value=500.0)
        add_unparsable_card(p)
    ic = HeaderCollection(str(raw_dir), keywords=LBC_KEYWORDS)
    make_bias(ic, image_directory=str(work_dir), raw_directory=str(raw_dir),
              verbose=False)
    assert 'PA_PNT' not in fits.getheader(work_dir / 'zero.fits')


def test_targetdirs_object_name_with_plus(work_dir):
    """Object names with '+' or '-' reach their filter directory.

    The old code regex-escaped the name before files_filtered, which with
    ccdproc >= 2.5 (which escapes it again) matched nothing, so the
    per-chip files were never sorted into the filter directory.
    """
    from LBCgo.headers import HeaderCollection
    from LBCgo.lbcproc import make_targetdirectories
    name = 'lbcb.20230101.000001_flat.fits'
    write_lbc_file(work_dir, name, object_name='J1419+4207-a',
                   nx=NX_SCIENCE)
    ic = HeaderCollection(str(work_dir), keywords=LBC_KEYWORDS,
                          filenames=[name])
    _, fltr_dirs = make_targetdirectories(
        ic, image_directory=str(work_dir) + '/', verbose=False)
    assert (work_dir / 'J1419+4207-a' / 'g-SLOAN' / name).is_file()
    assert len(fltr_dirs) == 1
