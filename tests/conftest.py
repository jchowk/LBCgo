"""
Shared fixtures and helpers for the LBCgo test suite.

All tests use synthetic LBC-style FITS files built by make_lbc_hdu() rather
than real data.  Chip dimensions are kept small (100 cols × 50 rows) so the
suite runs in a few seconds.

BIASSEC = [81:100,1:50]  →  20 overscan columns (cols 81–100, 1-indexed)
TRIMSEC = [1:80,1:50]    →  80 science columns  (cols 1–80,  1-indexed)
After go_overscan the science region is 80 cols × 50 rows.
"""

import numpy as np
import pytest
from pathlib import Path
from astropy.io import fits
from ccdproc import ImageFileCollection

# ---------------------------------------------------------------------------
# Constants shared across tests
# ---------------------------------------------------------------------------
NX = 100          # total columns per chip (science + overscan)
NY = 50           # rows per chip
NX_SCIENCE = 80   # columns after overscan trim
N_CHIPS = 4
LBC_KEYWORDS = ['object', 'filter', 'exptime', 'imagetyp',
                'propid', 'lbcobnam', 'airmass', 'HA', 'objra', 'objdec']


# ---------------------------------------------------------------------------
# Core helper: build a synthetic LBC multi-extension FITS HDUList
# ---------------------------------------------------------------------------
def make_lbc_hdu(imagetyp='object',
                 filter_name='g-SLOAN',
                 object_name='NGC891',
                 lbcobnam='NGC891',
                 n_chips=N_CHIPS,
                 nx=NX,
                 ny=NY,
                 pixel_value=5000.0,
                 exptime=300.0,
                 airmass=1.2):
    """Return a synthetic LBC HDUList.

    Extension 0: primary header only (no data).
    Extensions 1–n_chips: chip images, each (ny, nx) float32.

    The overscan region occupies columns 81–100 (1-indexed FITS convention):
      BIASSEC = '[81:100,1:{ny}]'
      TRIMSEC = '[1:80,1:{ny}]'
    """
    primary_hdr = fits.Header()
    primary_hdr['IMAGETYP'] = imagetyp
    primary_hdr['FILTER']   = filter_name
    primary_hdr['OBJECT']   = object_name
    primary_hdr['LBCOBNAM'] = lbcobnam
    primary_hdr['EXPTIME']  = exptime
    primary_hdr['AIRMASS']  = airmass
    primary_hdr['OBJRA']    = 35.637
    primary_hdr['OBJDEC']   = 42.350
    primary_hdr['PROPID']   = 'TEST'
    primary_hdr['HA']       = 0.5
    primary_hdr['NEXTEND']  = n_chips
    primary_hdu = fits.PrimaryHDU(header=primary_hdr)

    hdul = fits.HDUList([primary_hdu])

    rng = np.random.default_rng(seed=42)
    for chip in range(1, n_chips + 1):
        # Science region: cols 0–79 (0-indexed); overscan: cols 80–99
        data = np.full((ny, nx), pixel_value, dtype='float32')
        # Add small Poisson-like noise to science region
        data[:, :NX_SCIENCE] += rng.normal(0, pixel_value * 0.01,
                                            size=(ny, NX_SCIENCE)).astype('float32')
        # Overscan region: low, near-constant value
        data[:, NX_SCIENCE:] = 50.0 + rng.normal(0, 2.0,
                                                   size=(ny, nx - NX_SCIENCE)).astype('float32')

        hdr = fits.Header()
        hdr['EXTNAME']  = f'chip{chip}'
        hdr['CHIPID']   = chip
        # Per-chip FILTER is truncated to 8 chars in real LBC headers
        hdr['FILTER']   = filter_name[:8]
        hdr['IMAGETYP'] = imagetyp
        hdr['BIASSEC']  = f'[81:100,1:{ny}]'
        hdr['TRIMSEC']  = f'[1:80,1:{ny}]'
        hdr['DATASEC']  = f'[1:80,1:{ny}]'
        hdr['ROTANGLE'] = 0.0
        hdr['PARANGLE'] = 45.0
        # Duplicate astrometric keyword that go_overscan must remove
        hdr['CTYPE1A']  = 'RA---TAN'
        hdr['BZERO']    = 32768.0
        hdr['BUNIT']    = 'adu'   # required by CCDData.read(..., unit=None)

        hdul.append(fits.ImageHDU(data=data, header=hdr))

    return hdul


def write_lbc_file(directory, filename, **kwargs):
    """Write a synthetic LBC FITS file; return its Path."""
    hdul = make_lbc_hdu(**kwargs)
    path = Path(directory) / filename
    hdul.writeto(str(path), overwrite=True)
    hdul.close()
    return path


# ---------------------------------------------------------------------------
# Directory fixtures
# ---------------------------------------------------------------------------
@pytest.fixture
def raw_dir(tmp_path):
    """A tmp raw/ subdirectory."""
    d = tmp_path / 'raw'
    d.mkdir()
    return d


@pytest.fixture
def work_dir(tmp_path):
    """The top-level tmp work directory (output goes here)."""
    return tmp_path


# ---------------------------------------------------------------------------
# File fixtures
# ---------------------------------------------------------------------------
@pytest.fixture
def one_object_file(raw_dir):
    """One synthetic object frame in raw_dir."""
    return write_lbc_file(raw_dir, 'lbcb.20230101.000001.fits',
                          imagetyp='object', filter_name='g-SLOAN',
                          object_name='NGC891', pixel_value=5000.0)


@pytest.fixture
def three_bias_files(raw_dir):
    """Three synthetic zero (bias) frames in raw_dir."""
    paths = []
    for i in range(1, 4):
        p = write_lbc_file(raw_dir, f'lbcb.20230101.00010{i}.fits',
                           imagetyp='zero', filter_name='g-SLOAN',
                           object_name='Bias', lbcobnam='Bias',
                           pixel_value=500.0)
        paths.append(p)
    return paths


@pytest.fixture
def three_flat_files(raw_dir):
    """Three synthetic sky-flat frames (g-SLOAN) in raw_dir."""
    paths = []
    for i in range(1, 4):
        p = write_lbc_file(raw_dir, f'lbcb.20230101.00020{i}.fits',
                           imagetyp='flat', filter_name='g-SLOAN',
                           object_name='SkyFlat', lbcobnam='SkyFlat',
                           pixel_value=20000.0)
        paths.append(p)
    return paths


# ---------------------------------------------------------------------------
# ImageFileCollection fixtures
# ---------------------------------------------------------------------------
@pytest.fixture
def ic_objects(raw_dir, one_object_file):
    """ImageFileCollection containing only the object file."""
    return ImageFileCollection(str(raw_dir), keywords=LBC_KEYWORDS)


@pytest.fixture
def ic_with_flats(raw_dir, one_object_file, three_flat_files):
    """ImageFileCollection containing object + flat files."""
    return ImageFileCollection(str(raw_dir), keywords=LBC_KEYWORDS)


@pytest.fixture
def ic_bias_only(raw_dir, three_bias_files):
    """ImageFileCollection containing only bias frames."""
    return ImageFileCollection(str(raw_dir), keywords=LBC_KEYWORDS)


# ---------------------------------------------------------------------------
# Synthetic master flat helper
# ---------------------------------------------------------------------------
def write_master_flat(directory, filter_name='g-SLOAN', pixel_value=20000.0,
                      n_chips=N_CHIPS, nx=NX_SCIENCE, ny=NY):
    """Write a ready-to-use master flat FITS file (no overscan region).

    Simulates the output of make_flatfield() with uniform, well-defined data.
    Uses NX_SCIENCE columns (80) since a real master flat has already been trimmed.
    """
    flat_path = Path(directory) / f'flat.{filter_name}.fits'

    primary_hdr = fits.Header()
    primary_hdr['FILTER']   = filter_name
    primary_hdr['IMAGETYP'] = 'flat'
    primary_hdr['NCOMBINE'] = 3
    primary_hdu = fits.PrimaryHDU(header=primary_hdr)
    hdul = fits.HDUList([primary_hdu])

    for chip in range(1, n_chips + 1):
        data = np.full((ny, nx), pixel_value, dtype='float32')
        hdr = fits.Header()
        hdr['FILTER']   = filter_name[:8]
        hdr['IMAGETYP'] = 'flat'
        hdr['BUNIT']    = 'adu'   # required by CCDData.read(..., unit=None)
        hdul.append(fits.ImageHDU(data=data, header=hdr))

    hdul.writeto(str(flat_path), overwrite=True)
    hdul.close()
    return flat_path


# ---------------------------------------------------------------------------
# Pipeline-output fixtures (run a stage and return the results)
# ---------------------------------------------------------------------------
@pytest.fixture
def overscan_files(ic_with_flats, work_dir, raw_dir):
    """Run go_overscan on object files; return list of _over.fits paths."""
    from LBCgo.lbcproc import go_overscan
    return go_overscan(ic_with_flats,
                       image_directory=str(work_dir) + '/',
                       raw_directory=str(raw_dir) + '/',
                       verbose=False,
                       return_files=True)
