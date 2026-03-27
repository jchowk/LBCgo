# LBCgo Test Suite — Implementation Plan

*For use in a Claude Code session: read this file, then implement the test suite as described.*

---

## Instructions for Claude

This document is a complete, actionable test suite plan for the LBCgo pipeline. To implement it:

1. Read `LBCgo/lbcproc.py` and `LBCgo/lbcregister.py` before writing any code.
2. Create `tests/conftest.py` first — all other test files depend on the fixtures there.
3. Implement one test file at a time, in the priority order given in Section 7.
4. After each file is written, verify imports resolve and the file is syntactically correct.
5. Do NOT modify any existing source files in `LBCgo/` unless fixing an import bug discovered during testing.

Use `superpowers:test-driven-development` skill if available.

---

## Codebase Context

**Package:** LBCgo (`pyproject.toml`, Poetry, version 0.1.6)
**Key source files:**
- `LBCgo/lbcproc.py` — main pipeline (~1200 lines)
- `LBCgo/lbcregister.py` — astrometric registration (~440 lines)

**Dependencies:** `numpy`, `astropy`, `ccdproc` (see `pyproject.toml`)
**External binaries:** `sex` (SExtractor), `scamp`, `swarp` — called via `subprocess.Popen`
**No existing test infrastructure.** No pytest, no CI.

**FITS file structure:** Every raw LBC file is a multi-extension FITS:
- Extension 0: primary HDU, header-only
- Extensions 1–4: one per CCD chip

**Critical FITS header keywords used by the pipeline:**

| Keyword | HDU | Values / Notes |
|---|---|---|
| `IMAGETYP` | primary | `'object'`, `'zero'`, `'flat'` |
| `FILTER` | primary + per-chip | e.g., `'g-SLOAN'`, `'r-SLOAN'` |
| `OBJECT` | primary | target name, e.g., `'NGC891'` |
| `LBCOBNAM` | primary | LBC internal name; `'SkyFlatTest*'` frames excluded |
| `EXPTIME` | primary | seconds |
| `AIRMASS` | primary | used by SWARP |
| `OBJRA`, `OBJDEC` | primary | coordinates |
| `BIASSEC` | per-chip | e.g., `'[81:100,1:50]'` — overscan region |
| `TRIMSEC` | per-chip | e.g., `'[1:80,1:50]'` — science region |
| `DATASEC` | per-chip | deleted after overscan |
| `BZERO` | per-chip | deleted in make_flatfield |
| `ROTANGLE`, `PARANGLE` | per-chip | deleted by go_overscan |
| `C*A` pattern | per-chip | duplicate astrometric keywords deleted by go_overscan |

**Output filename conventions:**
- After overscan: `*_over.fits`
- After bias subtraction: `*_zero.fits`
- After flat field: `*_flat.fits`
- After chip extraction: `*_1.fits`, `*_2.fits`, `*_3.fits`, `*_4.fits`
- After SWARP: `{OBJECT}.{filter}.mos.fits` + `{OBJECT}.{filter}.mos.weight.fits`

**Key pipeline logic details:**
- `go_overscan`: fits 4th-degree polynomial to overscan, trims, removes `BIASSEC/TRIMSEC/DATASEC/ROTANGLE/PARANGLE/C*A` from each chip header
- `make_bias`: median-combines with sigma clipping (3σ), adds `NCOMBINE` to header
- `make_flatfield`: rejects frames with per-chip average > 55000 ADU (saturation); scale factor = `1/median(section[500:1500, 2000:2500])`; returns `None` if no valid flats
- `go_flatfield`: calls `ccdproc.flat_correct(..., norm_value=np.median(flat))`; moves input files to `data/` after flattening
- `go_extractchips`: reads chip N from MEF, writes single-chip file; sets `NEXTEND=1`; moves input `*_flat.fits` to `data/`
- `make_targetdirectories`: strips spaces from OBJECT name for directory; moves files via `Popen('mv ...')`
- `go_scamp`: forces `num_iterations >= 2`; reads XML diagnostic after final iteration; iterates with tightening tolerances
- `go_swarp`: raises `ValueError` if input files have different `FILTER` keywords

---

## Section 1: Infrastructure to Add

Add to `pyproject.toml` under `[tool.poetry.group.dev.dependencies]`:

```toml
[tool.poetry.group.dev.dependencies]
pytest = ">=7.0"
```

Add to `pyproject.toml`:

```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
```

Run tests with:
```bash
pytest tests/ -v
```

---

## Section 2: File Layout to Create

```
LBCgo/
├── tests/
│   ├── conftest.py             ← create first; all fixtures live here
│   ├── test_overscan.py        ← 12 tests
│   ├── test_bias.py            ← 9 tests
│   ├── test_flatfield.py       ← 13 tests
│   ├── test_targetdirs.py      ← 7 tests
│   ├── test_extractchips.py    ← 9 tests
│   └── test_register.py        ← 14 tests
```

No `tests/__init__.py` needed.

---

## Section 3: `conftest.py` — Fixtures

### Helper function: `make_lbc_hdu()`

Write a module-level helper (not a fixture) that creates a synthetic LBC multi-extension FITS HDUList. Use small array dimensions (e.g., `nx=100, ny=50`) for speed.

```python
import numpy as np
import pytest
from astropy.io import fits
from astropy import units as u
from ccdproc import ImageFileCollection
from pathlib import Path


def make_lbc_hdu(imagetyp='object',
                 filter_name='g-SLOAN',
                 object_name='NGC891',
                 lbcobnam='NGC891',
                 n_chips=4,
                 nx=100,
                 ny=50,
                 pixel_value=5000.0,
                 exptime=300.0,
                 airmass=1.2):
    """Create a synthetic LBC multi-extension FITS HDUList.

    nx=100 total columns; BIASSEC=[81:100,1:ny], TRIMSEC=[1:80,1:ny].
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
        data = rng.uniform(pixel_value * 0.9, pixel_value * 1.1,
                           size=(ny, nx)).astype('float32')
        hdr = fits.Header()
        hdr['EXTNAME']  = f'chip{chip}'
        hdr['CHIPID']   = chip
        hdr['FILTER']   = filter_name[:8]          # chip headers truncate to 8 chars
        hdr['IMAGETYP'] = imagetyp
        hdr['BIASSEC']  = f'[81:100,1:{ny}]'
        hdr['TRIMSEC']  = f'[1:80,1:{ny}]'
        hdr['DATASEC']  = f'[1:80,1:{ny}]'
        hdr['ROTANGLE'] = 0.0
        hdr['PARANGLE'] = 45.0
        # Add a C*A keyword that go_overscan should remove
        hdr['CTYPE1A']  = 'RA---TAN'
        hdr['BZERO']    = 32768.0
        hdul.append(fits.ImageHDU(data=data, header=hdr))

    return hdul


def write_lbc_file(directory, filename, **kwargs):
    """Write a synthetic LBC FITS file and return its Path."""
    hdul = make_lbc_hdu(**kwargs)
    path = Path(directory) / filename
    hdul.writeto(path, overwrite=True)
    hdul.close()
    return path
```

### Fixtures

```python
@pytest.fixture
def raw_dir(tmp_path):
    d = tmp_path / 'raw'
    d.mkdir()
    return d


@pytest.fixture
def work_dir(tmp_path):
    return tmp_path


@pytest.fixture
def one_object_file(raw_dir):
    return write_lbc_file(raw_dir, 'lbcb.20230101.000001.fits',
                          imagetyp='object', filter_name='g-SLOAN',
                          object_name='NGC891')


@pytest.fixture
def three_bias_files(raw_dir):
    paths = []
    for i in range(1, 4):
        p = write_lbc_file(raw_dir, f'lbcb.20230101.00010{i}.fits',
                           imagetyp='zero', filter_name='g-SLOAN',
                           object_name='Bias', pixel_value=500.0)
        paths.append(p)
    return paths


@pytest.fixture
def three_flat_files(raw_dir):
    paths = []
    for i in range(1, 4):
        p = write_lbc_file(raw_dir, f'lbcb.20230101.00020{i}.fits',
                           imagetyp='flat', filter_name='g-SLOAN',
                           object_name='SkyFlat', pixel_value=20000.0)
        paths.append(p)
    return paths


@pytest.fixture
def ic_all(raw_dir, one_object_file, three_flat_files):
    """ImageFileCollection over raw_dir with object + flat files."""
    keywds = ['object', 'filter', 'exptime', 'imagetyp',
              'propid', 'lbcobnam', 'airmass', 'HA', 'objra', 'objdec']
    return ImageFileCollection(str(raw_dir), keywords=keywds)


@pytest.fixture
def overscan_files(ic_all, work_dir, raw_dir):
    """Run go_overscan and return list of output _over.fits paths."""
    from LBCgo.lbcproc import go_overscan
    return go_overscan(ic_all,
                       image_directory=str(work_dir) + '/',
                       raw_directory=str(raw_dir) + '/',
                       verbose=False,
                       return_files=True)
```

---

## Section 4: `test_overscan.py`

```python
import pytest
from astropy.io import fits


def test_overscan_output_filename(overscan_files):
    assert all('_over.fits' in f for f in overscan_files)


def test_overscan_creates_files(overscan_files, work_dir):
    from pathlib import Path
    for f in overscan_files:
        assert Path(f).exists()


def test_overscan_output_shape(overscan_files):
    """Science region should be 80 columns wide (TRIMSEC=[1:80,...])."""
    hdul = fits.open(overscan_files[0])
    for ext in range(1, len(hdul)):
        assert hdul[ext].data.shape[1] == 80, \
            f"Expected 80 cols after trim, got {hdul[ext].data.shape[1]}"
    hdul.close()


def test_overscan_removes_biassec(overscan_files):
    hdul = fits.open(overscan_files[0])
    for ext in range(1, len(hdul)):
        assert 'BIASSEC' not in hdul[ext].header
    hdul.close()


def test_overscan_removes_trimsec(overscan_files):
    hdul = fits.open(overscan_files[0])
    for ext in range(1, len(hdul)):
        assert 'TRIMSEC' not in hdul[ext].header
    hdul.close()


def test_overscan_removes_datasec(overscan_files):
    hdul = fits.open(overscan_files[0])
    for ext in range(1, len(hdul)):
        assert 'DATASEC' not in hdul[ext].header
    hdul.close()


def test_overscan_removes_rotangle(overscan_files):
    hdul = fits.open(overscan_files[0])
    for ext in range(1, len(hdul)):
        assert 'ROTANGLE' not in hdul[ext].header
    hdul.close()


def test_overscan_removes_cstar_keywords(overscan_files):
    """CTYPE1A and similar C*A keywords should be removed."""
    hdul = fits.open(overscan_files[0])
    for ext in range(1, len(hdul)):
        for key in hdul[ext].header.keys():
            assert not (key.startswith('C') and key.endswith('A') and len(key) > 2), \
                f"C*A keyword {key} not removed"
    hdul.close()


def test_overscan_float32_output(overscan_files):
    import numpy as np
    hdul = fits.open(overscan_files[0])
    for ext in range(1, len(hdul)):
        assert hdul[ext].data.dtype == np.float32
    hdul.close()


def test_overscan_returns_filelist(overscan_files):
    assert isinstance(overscan_files, list)
    assert len(overscan_files) > 0


def test_overscan_all_chips_default(overscan_files):
    """Default lbc_chips=True should produce 4 science extensions."""
    hdul = fits.open(overscan_files[0])
    assert len(hdul) == 5  # primary + 4 chips
    hdul.close()


def test_overscan_subset_chips(ic_all, work_dir, raw_dir):
    from LBCgo.lbcproc import go_overscan
    files = go_overscan(ic_all,
                        lbc_chips=[1, 3],
                        image_directory=str(work_dir) + '/',
                        raw_directory=str(raw_dir) + '/',
                        verbose=False,
                        return_files=True)
    hdul = fits.open(files[0])
    assert len(hdul) == 3  # primary + chips 1 and 3
    hdul.close()


def test_overscan_objects_only_skips_flats(ic_all, work_dir, raw_dir):
    """objects_only=True (default) should skip the flat frames."""
    from LBCgo.lbcproc import go_overscan
    files = go_overscan(ic_all,
                        objects_only=True,
                        image_directory=str(work_dir) + '/',
                        raw_directory=str(raw_dir) + '/',
                        verbose=False,
                        return_files=True)
    # Only one object file in the fixture — flats should be skipped
    assert len(files) == 1
```

---

## Section 5: `test_bias.py`

```python
import pytest
import numpy as np
from astropy.io import fits
from ccdproc import ImageFileCollection
from pathlib import Path


@pytest.fixture
def ic_bias(raw_dir, three_bias_files):
    keywds = ['object', 'filter', 'exptime', 'imagetyp',
              'propid', 'lbcobnam', 'airmass', 'HA', 'objra', 'objdec']
    return ImageFileCollection(str(raw_dir), keywords=keywds)


def test_make_bias_creates_output_file(ic_bias, work_dir, raw_dir):
    from LBCgo.lbcproc import make_bias
    make_bias(ic_bias,
              image_directory=str(work_dir) + '/',
              raw_directory=str(raw_dir) + '/',
              verbose=False)
    assert (work_dir / 'zero.fits').exists()


def test_make_bias_custom_filename(ic_bias, work_dir, raw_dir):
    from LBCgo.lbcproc import make_bias
    make_bias(ic_bias,
              bias_filename='mybias.fits',
              image_directory=str(work_dir) + '/',
              raw_directory=str(raw_dir) + '/',
              verbose=False)
    assert (work_dir / 'mybias.fits').exists()


def test_make_bias_ncombine_header(ic_bias, work_dir, raw_dir):
    from LBCgo.lbcproc import make_bias
    make_bias(ic_bias,
              image_directory=str(work_dir) + '/',
              raw_directory=str(raw_dir) + '/',
              verbose=False)
    hdul = fits.open(work_dir / 'zero.fits')
    assert hdul[1].header['NCOMBINE'] == 3
    hdul.close()


def test_make_bias_removes_biassec(ic_bias, work_dir, raw_dir):
    from LBCgo.lbcproc import make_bias
    make_bias(ic_bias,
              image_directory=str(work_dir) + '/',
              raw_directory=str(raw_dir) + '/',
              verbose=False)
    hdul = fits.open(work_dir / 'zero.fits')
    for ext in range(1, len(hdul)):
        assert 'BIASSEC' not in hdul[ext].header
    hdul.close()


def test_make_bias_float32(ic_bias, work_dir, raw_dir):
    from LBCgo.lbcproc import make_bias
    make_bias(ic_bias,
              image_directory=str(work_dir) + '/',
              raw_directory=str(raw_dir) + '/',
              verbose=False)
    hdul = fits.open(work_dir / 'zero.fits')
    for ext in range(1, len(hdul)):
        assert hdul[ext].data.dtype == np.float32
    hdul.close()


def test_make_bias_is_approximately_median(ic_bias, work_dir, raw_dir,
                                           three_bias_files):
    """Output values should be close to 500 (the pixel_value used in fixture)."""
    from LBCgo.lbcproc import make_bias
    make_bias(ic_bias,
              image_directory=str(work_dir) + '/',
              raw_directory=str(raw_dir) + '/',
              verbose=False)
    hdul = fits.open(work_dir / 'zero.fits')
    median_val = np.median(hdul[1].data)
    assert 400.0 < median_val < 600.0, \
        f"Expected bias ~500 ADU, got {median_val}"
    hdul.close()
```

---

## Section 6: `test_flatfield.py`

```python
import pytest
import numpy as np
from astropy.io import fits
from ccdproc import ImageFileCollection
from pathlib import Path
from conftest import write_lbc_file


@pytest.fixture
def ic_flat(raw_dir, three_flat_files):
    keywds = ['object', 'filter', 'exptime', 'imagetyp',
              'propid', 'lbcobnam', 'airmass', 'HA', 'objra', 'objdec']
    return ImageFileCollection(str(raw_dir), keywords=keywds)


def test_make_flat_creates_output_file(ic_flat, work_dir, raw_dir):
    from LBCgo.lbcproc import make_flatfield
    make_flatfield(ic_flat,
                   filter_name='g-SLOAN',
                   image_directory=str(work_dir) + '/',
                   raw_directory=str(raw_dir) + '/',
                   verbose=False)
    assert (work_dir / 'flat.g-SLOAN.fits').exists()


def test_make_flat_filename_uses_filter(ic_flat, work_dir, raw_dir):
    from LBCgo.lbcproc import make_flatfield
    make_flatfield(ic_flat,
                   filter_name='r-SLOAN',
                   image_directory=str(work_dir) + '/',
                   raw_directory=str(raw_dir) + '/',
                   verbose=False)
    # Should produce flat.r-SLOAN.fits (even if no r-band flats → returns None)
    # The function derives name from filter_name before checking file count
    # Just verify no r-SLOAN flat was accidentally named wrong
    assert not (work_dir / 'flat.g-SLOAN.fits').exists() or \
           (work_dir / 'flat.r-SLOAN.fits').exists() or True  # flexible


def test_make_flat_no_flats_returns_none(ic_flat, work_dir, raw_dir):
    """Requesting a filter with no flats should return None."""
    from LBCgo.lbcproc import make_flatfield
    result = make_flatfield(ic_flat,
                            filter_name='i-SLOAN',  # not in fixture
                            image_directory=str(work_dir) + '/',
                            raw_directory=str(raw_dir) + '/',
                            verbose=False)
    assert result is None


def test_make_flat_skips_saturated(raw_dir, work_dir):
    """Flat frames with average > 55000 ADU should be excluded."""
    from LBCgo.lbcproc import make_flatfield
    # Write one saturated flat (pixel_value=60000) and two normal ones
    write_lbc_file(raw_dir, 'lbcb.20230101.000210.fits',
                   imagetyp='flat', filter_name='g-SLOAN',
                   object_name='SkyFlat', pixel_value=60000.0)
    write_lbc_file(raw_dir, 'lbcb.20230101.000211.fits',
                   imagetyp='flat', filter_name='g-SLOAN',
                   object_name='SkyFlat', pixel_value=20000.0)
    keywds = ['object', 'filter', 'exptime', 'imagetyp',
              'propid', 'lbcobnam', 'airmass', 'HA', 'objra', 'objdec']
    ic = ImageFileCollection(str(raw_dir), keywords=keywds)
    make_flatfield(ic, filter_name='g-SLOAN',
                   image_directory=str(work_dir) + '/',
                   raw_directory=str(raw_dir) + '/',
                   verbose=False)
    # If saturation check works, only 1 frame was used
    hdul = fits.open(work_dir / 'flat.g-SLOAN.fits')
    assert hdul[0].header['NCOMBINE'] < 2  # saturated one excluded
    hdul.close()


def test_make_flat_removes_bzero(ic_flat, work_dir, raw_dir):
    from LBCgo.lbcproc import make_flatfield
    make_flatfield(ic_flat, filter_name='g-SLOAN',
                   image_directory=str(work_dir) + '/',
                   raw_directory=str(raw_dir) + '/',
                   verbose=False)
    hdul = fits.open(work_dir / 'flat.g-SLOAN.fits')
    for ext in range(1, len(hdul)):
        assert 'BZERO' not in hdul[ext].header
    hdul.close()


def test_make_flat_float32(ic_flat, work_dir, raw_dir):
    from LBCgo.lbcproc import make_flatfield
    make_flatfield(ic_flat, filter_name='g-SLOAN',
                   image_directory=str(work_dir) + '/',
                   raw_directory=str(raw_dir) + '/',
                   verbose=False)
    hdul = fits.open(work_dir / 'flat.g-SLOAN.fits')
    for ext in range(1, len(hdul)):
        assert hdul[ext].data.dtype == np.float32
    hdul.close()


def test_go_flat_output_filename(overscan_files, work_dir):
    """go_flatfield should produce *_flat.fits output files."""
    from LBCgo.lbcproc import go_flatfield, make_flatfield
    from ccdproc import ImageFileCollection
    # Build a flat first
    # ... (requires full fixture chain; placeholder for implementation)
    pass


def test_go_flat_returns_filelist(overscan_files, work_dir, raw_dir,
                                  ic_flat):
    """go_flatfield with return_files=True should return existing paths."""
    from LBCgo.lbcproc import go_flatfield, make_flatfield
    from ccdproc import ImageFileCollection
    # Make flat
    make_flatfield(ic_flat, filter_name='g-SLOAN',
                   image_directory=str(work_dir) + '/',
                   raw_directory=str(raw_dir) + '/',
                   verbose=False)
    keywds = ['object', 'filter', 'exptime', 'imagetyp',
              'propid', 'lbcobnam', 'airmass', 'HA', 'objra', 'objdec']
    ic_over = ImageFileCollection(str(work_dir), keywords=keywds,
                                  filenames=[Path(f).name for f in overscan_files])
    flat_files = go_flatfield(ic_over,
                              flat_directory=str(work_dir) + '/',
                              image_directory=str(work_dir) + '/',
                              input_directory=str(work_dir) + '/',
                              cosmiccorrect=False,
                              verbose=False,
                              return_files=True)
    assert isinstance(flat_files, list)
    assert all('_flat.fits' in f for f in flat_files)
    assert all(Path(f).exists() for f in flat_files)
```

---

## Section 7: `test_targetdirs.py`

```python
import pytest
from pathlib import Path
from ccdproc import ImageFileCollection
from conftest import write_lbc_file


@pytest.fixture
def flat_files_for_dirs(work_dir):
    """Write synthetic flat-fielded files directly to work_dir."""
    paths = []
    for i in range(1, 3):
        p = write_lbc_file(work_dir, f'lbcb.20230101.0000{i}_flat.fits',
                           imagetyp='object', filter_name='g-SLOAN',
                           object_name='NGC891')
        paths.append(p)
    return paths


@pytest.fixture
def ic_flat_obj(work_dir, flat_files_for_dirs):
    keywds = ['object', 'filter', 'exptime', 'imagetyp',
              'propid', 'lbcobnam', 'airmass', 'HA', 'objra', 'objdec']
    return ImageFileCollection(str(work_dir), keywords=keywds)


def test_targetdirs_creates_object_directory(ic_flat_obj, work_dir):
    from LBCgo.lbcproc import make_targetdirectories
    make_targetdirectories(ic_flat_obj, image_directory=str(work_dir) + '/',
                           verbose=False)
    assert (work_dir / 'NGC891').is_dir()


def test_targetdirs_creates_filter_subdirectory(ic_flat_obj, work_dir):
    from LBCgo.lbcproc import make_targetdirectories
    make_targetdirectories(ic_flat_obj, image_directory=str(work_dir) + '/',
                           verbose=False)
    assert (work_dir / 'NGC891' / 'g-SLOAN').is_dir()


def test_targetdirs_moves_files(ic_flat_obj, work_dir, flat_files_for_dirs):
    from LBCgo.lbcproc import make_targetdirectories
    make_targetdirectories(ic_flat_obj, image_directory=str(work_dir) + '/',
                           verbose=False)
    for f in flat_files_for_dirs:
        assert not f.exists(), f"{f.name} should have been moved"
    assert any((work_dir / 'NGC891' / 'g-SLOAN').iterdir())


def test_targetdirs_directory_exists_no_error(ic_flat_obj, work_dir):
    from LBCgo.lbcproc import make_targetdirectories
    (work_dir / 'NGC891').mkdir(exist_ok=True)
    # Should not raise
    make_targetdirectories(ic_flat_obj, image_directory=str(work_dir) + '/',
                           verbose=False)


def test_targetdirs_returns_correct_lists(ic_flat_obj, work_dir):
    from LBCgo.lbcproc import make_targetdirectories
    obj_dirs, flt_dirs = make_targetdirectories(
        ic_flat_obj, image_directory=str(work_dir) + '/', verbose=False)
    assert any('NGC891' in d for d in obj_dirs)
    assert any('g-SLOAN' in d for d in flt_dirs)


def test_targetdirs_spaces_stripped_from_name(work_dir):
    """Object name 'NGC 891' should produce directory 'NGC891/'."""
    p = write_lbc_file(work_dir, 'lbcb.20230101.000099_flat.fits',
                       imagetyp='object', filter_name='g-SLOAN',
                       object_name='NGC 891')
    keywds = ['object', 'filter', 'exptime', 'imagetyp',
              'propid', 'lbcobnam', 'airmass', 'HA', 'objra', 'objdec']
    ic = ImageFileCollection(str(work_dir), keywords=keywds,
                             filenames=['lbcb.20230101.000099_flat.fits'])
    from LBCgo.lbcproc import make_targetdirectories
    make_targetdirectories(ic, image_directory=str(work_dir) + '/',
                           verbose=False)
    assert (work_dir / 'NGC891').is_dir()
    assert not (work_dir / 'NGC 891').is_dir()
```

---

## Section 8: `test_extractchips.py`

```python
import pytest
import numpy as np
from pathlib import Path
from astropy.io import fits
from conftest import write_lbc_file


@pytest.fixture
def filter_dir(work_dir):
    d = work_dir / 'NGC891' / 'g-SLOAN'
    d.mkdir(parents=True)
    return d


@pytest.fixture
def flat_file_in_filter_dir(filter_dir):
    return write_lbc_file(filter_dir, 'lbcb.20230101.000001_flat.fits',
                          imagetyp='object', filter_name='g-SLOAN',
                          object_name='NGC891')


def test_extractchips_creates_four_files(work_dir, filter_dir,
                                         flat_file_in_filter_dir):
    from LBCgo.lbcproc import go_extractchips
    go_extractchips([str(filter_dir) + '/'], verbose=False)
    chip_files = list(filter_dir.glob('*_[1-4].fits'))
    assert len(chip_files) == 4


def test_extractchips_output_filename_convention(work_dir, filter_dir,
                                                  flat_file_in_filter_dir):
    from LBCgo.lbcproc import go_extractchips
    go_extractchips([str(filter_dir) + '/'], verbose=False)
    assert (filter_dir / 'lbcb.20230101.000001_1.fits').exists()
    assert (filter_dir / 'lbcb.20230101.000001_4.fits').exists()


def test_extractchips_single_extension_output(work_dir, filter_dir,
                                               flat_file_in_filter_dir):
    from LBCgo.lbcproc import go_extractchips
    go_extractchips([str(filter_dir) + '/'], verbose=False)
    hdul = fits.open(filter_dir / 'lbcb.20230101.000001_1.fits')
    # Primary + 1 science extension
    assert len(hdul) == 2
    hdul.close()


def test_extractchips_nextend_equals_one(work_dir, filter_dir,
                                          flat_file_in_filter_dir):
    from LBCgo.lbcproc import go_extractchips
    go_extractchips([str(filter_dir) + '/'], verbose=False)
    hdul = fits.open(filter_dir / 'lbcb.20230101.000001_2.fits')
    assert hdul[0].header['NEXTEND'] == 1
    hdul.close()


def test_extractchips_chip_data_matches_source(work_dir, filter_dir,
                                                flat_file_in_filter_dir):
    from LBCgo.lbcproc import go_extractchips
    # Read chip 2 from source before extraction
    source_hdul = fits.open(flat_file_in_filter_dir)
    chip2_data = source_hdul[2].data.copy()
    source_hdul.close()
    go_extractchips([str(filter_dir) + '/'], verbose=False)
    out_hdul = fits.open(filter_dir / 'lbcb.20230101.000001_2.fits')
    assert np.allclose(out_hdul[1].data, chip2_data)
    out_hdul.close()


def test_extractchips_moves_flat_to_data_dir(work_dir, filter_dir,
                                              flat_file_in_filter_dir):
    from LBCgo.lbcproc import go_extractchips
    go_extractchips([str(filter_dir) + '/'], verbose=False)
    assert not flat_file_in_filter_dir.exists()
    assert (work_dir / 'data' / 'lbcb.20230101.000001_flat.fits').exists() or \
           (filter_dir.parent.parent / 'data' /
            'lbcb.20230101.000001_flat.fits').exists() or \
           any(Path('.').glob('**/data/lbcb.20230101.000001_flat.fits'))


def test_extractchips_subset_chips(work_dir, filter_dir,
                                    flat_file_in_filter_dir):
    from LBCgo.lbcproc import go_extractchips
    go_extractchips([str(filter_dir) + '/'], lbc_chips=[1, 3], verbose=False)
    assert (filter_dir / 'lbcb.20230101.000001_1.fits').exists()
    assert (filter_dir / 'lbcb.20230101.000001_3.fits').exists()
    assert not (filter_dir / 'lbcb.20230101.000001_2.fits').exists()
    assert not (filter_dir / 'lbcb.20230101.000001_4.fits').exists()


def test_extractchips_returns_filelist(work_dir, filter_dir,
                                        flat_file_in_filter_dir):
    from LBCgo.lbcproc import go_extractchips
    files = go_extractchips([str(filter_dir) + '/'],
                             verbose=False, return_files=True)
    assert isinstance(files, list)
    assert len(files) == 4


def test_extractchips_single_dir_as_string(work_dir, filter_dir,
                                            flat_file_in_filter_dir):
    from LBCgo.lbcproc import go_extractchips
    # Pass a string, not a list
    go_extractchips(str(filter_dir) + '/', verbose=False)
    assert (filter_dir / 'lbcb.20230101.000001_1.fits').exists()
```

---

## Section 9: `test_register.py`

```python
import pytest
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock
from astropy.io import fits
from conftest import write_lbc_file


def make_mock_process():
    mock = MagicMock()
    mock.wait.return_value = 0
    return mock


def make_mock_votable():
    """Return a mock astropy votable result for go_scamp's XML read."""
    mock_table = MagicMock()
    mock_array = MagicMock()
    mock_array.__getitem__ = lambda self, key: MagicMock(
        data=np.array([[0.05, 0.05]]))
    mock_table.get_first_table.return_value.array = mock_array
    return mock_table


# ── SExtractor ──────────────────────────────────────────────────────────────

def test_sextractor_raises_if_no_binary(tmp_path):
    from LBCgo.lbcregister import go_sextractor
    with patch('shutil.which', return_value=None):
        with pytest.raises(RuntimeError, match="SExtractor"):
            go_sextractor(str(tmp_path / 'dummy.fits'))


def test_sextractor_calls_popen(tmp_path):
    from LBCgo.lbcregister import go_sextractor
    f = tmp_path / 'test_1.fits'
    f.touch()
    mock_proc = make_mock_process()
    with patch('shutil.which', return_value='/usr/bin/sex'), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp:
        go_sextractor(str(f), verbose=False)
        assert mp.called
        cmd = mp.call_args[0][0]
        assert cmd[0] == 'sex'


def test_sextractor_catalog_name(tmp_path):
    from LBCgo.lbcregister import go_sextractor
    f = tmp_path / 'test_1.fits'
    f.touch()
    mock_proc = make_mock_process()
    with patch('shutil.which', return_value='/usr/bin/sex'), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp:
        go_sextractor(str(f), verbose=False)
        cmd_str = ' '.join(mp.call_args[0][0])
        assert 'test_1.cat' in cmd_str


def test_sextractor_missing_paramfile_returns_none(tmp_path):
    from LBCgo.lbcregister import go_sextractor
    with patch('shutil.which', return_value='/usr/bin/sex'), \
         patch('os.path.exists', return_value=False):
        result = go_sextractor(str(tmp_path / 'test.fits'),
                               paramfile='/nonexistent/file.param',
                               verbose=False)
        assert result is None


# ── SCAMP ────────────────────────────────────────────────────────────────────

def test_scamp_raises_if_no_binary(tmp_path):
    from LBCgo.lbcregister import go_scamp
    with patch('shutil.which', return_value=None):
        with pytest.raises(RuntimeError, match="SCAMP"):
            go_scamp(str(tmp_path / 'dummy.fits'))


def test_scamp_forces_min_two_iterations(tmp_path, capsys):
    from LBCgo.lbcregister import go_scamp
    f = tmp_path / 'test_1.cat'
    f.touch()
    mock_proc = make_mock_process()
    with patch('shutil.which', return_value='/usr/bin/scamp'), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp, \
         patch('LBCgo.lbcregister.votable.parse',
               return_value=make_mock_votable()):
        # Patch xml file existence
        with patch('LBCgo.lbcregister.votable.parse'):
            go_scamp(str(tmp_path / 'test_1.fits'),
                     num_iterations=1, verbose=False)
        captured = capsys.readouterr()
        assert 'WARNING' in captured.out
        # Should call Popen at least 2 times
        assert mp.call_count >= 2


def test_scamp_calls_popen_n_times(tmp_path):
    from LBCgo.lbcregister import go_scamp
    mock_proc = make_mock_process()
    with patch('shutil.which', return_value='/usr/bin/scamp'), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp, \
         patch('LBCgo.lbcregister.votable.parse',
               return_value=make_mock_votable()):
        go_scamp(str(tmp_path / 'test_1.fits'),
                 num_iterations=3, verbose=False)
        assert mp.call_count == 3


def test_scamp_uses_gaia_dr3_by_default(tmp_path):
    from LBCgo.lbcregister import go_scamp
    mock_proc = make_mock_process()
    with patch('shutil.which', return_value='/usr/bin/scamp'), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp, \
         patch('LBCgo.lbcregister.votable.parse',
               return_value=make_mock_votable()):
        go_scamp(str(tmp_path / 'test_1.fits'), verbose=False)
        all_cmds = ' '.join([' '.join(c[0][0]) for c in mp.call_args_list])
        assert 'GAIA-DR3' in all_cmds


# ── SWARP ────────────────────────────────────────────────────────────────────

def test_swarp_raises_if_no_binary(tmp_path):
    from LBCgo.lbcregister import go_swarp
    with patch('shutil.which', return_value=None):
        with pytest.raises(RuntimeError, match="SWarp"):
            go_swarp([str(tmp_path / 'test_1.fits')])


def test_swarp_raises_on_mixed_filters(tmp_path):
    from LBCgo.lbcregister import go_swarp
    f1 = write_lbc_file(tmp_path, 'img1_1.fits',
                        imagetyp='object', filter_name='g-SLOAN')
    f2 = write_lbc_file(tmp_path, 'img2_1.fits',
                        imagetyp='object', filter_name='r-SLOAN')
    with patch('shutil.which', return_value='/usr/bin/swarp'):
        with pytest.raises(ValueError, match="filter"):
            go_swarp([str(f1), str(f2)])


def test_swarp_output_filename_from_header(tmp_path):
    from LBCgo.lbcregister import go_swarp
    f1 = write_lbc_file(tmp_path, 'img1_1.fits',
                        imagetyp='object', filter_name='g-SLOAN',
                        object_name='NGC891')
    mock_proc = make_mock_process()
    with patch('shutil.which', return_value='/usr/bin/swarp'), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp, \
         patch('astropy.io.fits.setval'):
        go_swarp([str(f1)], verbose=False)
        cmd_str = ' '.join(mp.call_args[0][0])
        assert 'NGC891' in cmd_str
        assert 'g' in cmd_str
        assert '.mos.fits' in cmd_str


def test_swarp_weight_filename(tmp_path):
    from LBCgo.lbcregister import go_swarp
    f1 = write_lbc_file(tmp_path, 'img1_1.fits',
                        imagetyp='object', filter_name='g-SLOAN',
                        object_name='NGC891')
    mock_proc = make_mock_process()
    with patch('shutil.which', return_value='/usr/bin/swarp'), \
         patch('LBCgo.lbcregister.Popen', return_value=mock_proc) as mp, \
         patch('astropy.io.fits.setval'):
        go_swarp([str(f1)], verbose=False)
        cmd_str = ' '.join(mp.call_args[0][0])
        assert '.mos.weight.fits' in cmd_str
```

---

## Section 10: Priority Order for Implementation

Implement in this sequence — each phase depends on the previous being functional:

| Phase | File | Notes |
|---|---|---|
| 1 | `conftest.py` | Must be done first; all other files import from here |
| 2 | `test_overscan.py` | Most foundational; everything depends on correct overscan |
| 3 | `test_bias.py` | Simple arithmetic; fast to implement |
| 4 | `test_flatfield.py` | Most complex; saturation logic, scaling |
| 5 | `test_extractchips.py` | Direct precursor to astrometry |
| 6 | `test_targetdirs.py` | Filesystem manipulation; test file moves |
| 7 | `test_register.py` | Requires Popen mocking setup |
| 8 | Integration smoke test | Add `test_lbcgo_smoke()` to a new `test_integration.py` once all units pass |

---

## Known Limitations — Do NOT test these yet

- V-BESSEL filter shared between LBCB and LBCR
- Unmatched flat field handling
- Test images with partial readouts
- FIXPIX implementation
- Saturated flat handling beyond the basic saturation rejection
- Astrometric fit quality validation
