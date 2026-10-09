# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What LBCgo Is

LBCgo is a data reduction pipeline for the Large Binocular Telescope's Large Binocular Camera (LBC). It processes multi-extension FITS files through overscan removal, bias subtraction, flat fielding, chip extraction, and optional astrometric registration using the astromatic.net suite (SExtractor, SCAMP, SWARP).

## Install / Development

```bash
# Install from source
pip install -e .

# Or via conda environment
conda env create -f lbcgo_environment.yaml
```

External tools required for astrometric processing:
```bash
brew install sextractor scamp swarp     # macOS
conda install astromatic-scamp astromatic-swarp sextractor  # conda
```

## Tests

pytest suite in `tests/` (synthetic 100x50-px LBC-style FITS built by `make_lbc_hdu()` in `tests/conftest.py`; no real data). Use the lbcgo env:

```bash
/opt/miniconda3/envs/lbcgo/bin/python -m pytest tests -q                  # all (~30 s); also ./py_runtests.sh
/opt/miniconda3/envs/lbcgo/bin/python -m pytest tests/test_detector.py -k <name>   # single test
```

Give `PYTHONPATH` as an **absolute** path if you set it: `tests/test_calibration_example.py` runs `calibration/gain_rdnoise_example/run.py` in a subprocess from a temp dir. `tests/TEST_PLAN.md` is the plan. Tests of estimators use tolerances derived from the expected noise, checked over several seeds.

## Pipeline Architecture

The reduction pipeline is driven by a single top-level call:

```python
from LBCgo import lbcproc
lbcproc.lbcgo(data_dir, ...)
```

**Processing stages in `lbcproc.py`:**

1. `go_overscan()` — fit and subtract overscan regions, trim images
2. `make_bias()` / `go_bias()` — create and apply master bias frames
3. `make_flat()` / `go_flat()` — create per-filter flats and apply them
4. `go_createobjectdirectories()` — organize output by target name and filter
5. `go_extractchips()` — split multi-extension FITS into per-chip files
6. `go_register()` (optional, in `lbcregister.py`) — astrometric pipeline:
   - `go_sextractor()` → source catalogs
   - `go_scamp()` → astrometric calibration against GAIA-DR3
   - `go_swarp()` → resample and co-add

**Key files:**
- `LBCgo/lbcproc.py` — main reduction pipeline (~1500 lines)
- `LBCgo/lbcregister.py` — astrometric registration (~1200 lines; astromatic path, still the default)
- `LBCgo/detector.py` — per-chip gain/read noise: table lookup (`conf/lbc_detector.ecsv`) plus the photon-transfer measurement (`measure_ptc`, `measure_gain_rdnoise_files`, `summarize_gain_rdnoise`)
- `LBCgo/masks.py` — pixel masks and inverse-variance weight maps, written as `<image>.mask.fits` / `<image>.weight.fits` sidecars (`.weight.fits` matches the SExtractor/SWarp `WEIGHT_SUFFIX` default). Saturation is flagged in `lbcproc.go_overscan` on raw ADU.
- `LBCgo/conf/` — SExtractor/SCAMP/SWARP configs and `lbc_detector.ecsv`

**Detector table semantics:** rows keyed by (channel, chip) with `mjd_start` (inclusive) / `mjd_end` (exclusive), NaN = open-ended; the matching row with the latest `mjd_start` wins. No matching row → fall back to the `GAIN`/`RDNOISE` image headers. `gain` is per-pixel (variance/weights); `gain_flux` is for Poisson errors of summed fluxes. `INSTRUME` is spelled inconsistently (`LBC_BLUE` vs `LBC-RED `); use `detector.lbc_channel`.

## Calibration provenance (`calibration/`)

Not part of the installed package. One subdirectory per calibration product *version* (`gain_rdnoise_lbc_<YYYYMM>/`); new epoch → new directory, never edit an old one. Convention in `calibration/README.md`; start a README from `NEW_PRODUCT_TEMPLATE.md`. Each holds `README.md`, `inputs.ecsv` (archive filenames, not paths), `run.py`, outputs, and `run_log.json` (LBCgo commit, parameters, SHA-256s).

- `run.py` only orchestrates calls into `LBCgo` (logic and tests live in the package). Raw-data path comes from `$LBCGO_RAW` or `--raw`; raw frames are not in git. Template: `gain_rdnoise_example/run.py`. The gain/read-noise drivers are near-identical; only `PARAMS['mjd_start']` differs.
- Workflow: `python calibration/make_inputs.py $LBCGO_RAW -o inputs_full.ecsv --levels` → edit down to `inputs.ecsv` (each `set` = 2 flats + 2 biases of one channel) → `LBCGO_RAW=... python run.py --box 0` → review `detector_rows.ecsv` → merge into `LBCgo/conf/lbc_detector.ecsv` as a separate, reviewed step (a test checks the merged rows match each product's `detector_rows.ecsv`).
- Raw frames may have invalid header cards (e.g. `PA_PNT = nan` in 2010-03 LBCB); `ccdproc.ImageFileCollection` fails on these, `make_inputs.py` reads only the needed keywords. `lbcproc.py` still has this problem (see `LBCgo/00ToDo.md`).
- Photon-transfer facts: gain rises with flat level (brighter-fatter effect), so product `gain` is the zero-level intercept of a linear fit; use 10,000–30,000 ADU flats; run `gain_sum` on the whole chip (`--box 0`).

## Planning docs

`docs/planning/registration_migration_plan.md` is the master plan (§2 = PI decisions not to re-litigate; §5–6 checkboxes = current state). `docs/planning/claude_handoff.md` holds working notes and verified LBC facts (e.g. chip 4 is rotated only in its WCS; the raw readout layout and `overscan_axis=1` are the same for all chips). Read both before registration or detector work.

## Known Limitations (from README/ToDo)

- V-BESSEL filter is shared between LBCB and LBCR; handling requires separate pipeline runs per camera
- Unmatched flat fields (object has flat, calibration data doesn't, or vice versa) need careful handling
- Test images with partial readouts cause failures
- FIXPIX implementation and saturated flat handling are pending
- Astrometric fit quality validation is not yet implemented
- Background subtraction and astrometric fits with extended objects (e.g., large galaxies).

## External Subprocess Calls

`lbcregister.py` calls SExtractor, SCAMP, and SWARP via `subprocess.Popen()` (recently converted from `call()`). Configuration files are bundled in `LBCgo/conf/` and their paths are resolved at runtime via `importlib.resources` or equivalent.

## Packaging

- Version is defined in `pyproject.toml` (currently 0.1.6) — update there for releases
- Published on PyPI as `LBCgo`
- Uses Poetry as build backend; Python ≥ 3.11, numpy ≥ 2

## Development

- References for further development in `./references` directory. 
- docs/whatidid.LBCgo.2016.02.py contains a sample of reducing a specific dataset, including one in which chip 4 was not considered. Useful for DOCUMENTATION goal. 
