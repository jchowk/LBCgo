# <Product name and version, e.g. "LBCB static distortion model, v1">

Copy this file to `calibration/<product_dir>/README.md` and fill it in.
See [`README.md`](README.md) for the convention.

## Summary

- **Product file(s):** e.g. `LBCgo/conf/distortion/lbcb_v1.fits`
- **SHA-256:** output of `sha256sum <file>`
- **Made by / date:**
- **LBCgo commit:** output of `git rev-parse HEAD` when `run.py` was run
- **Supersedes:** previous product directory, if any
- **Valid for:** channel(s), filter(s), MJD range

## Purpose

Why this product was made and what it is used for in the pipeline.

## Inputs

`inputs.ecsv`, one row per frame, with at least these columns (generate a
starting table with `python calibration/make_inputs.py <raw_dir> -o
inputs.ecsv`, filtering by `--channel`, `--imagetyp`, `--object`,
`--filter`, `--propid`):

| column | meaning |
|---|---|
| `filename` | archive filename, e.g. `lbcb.20141120.065509.fits` |
| `obs_id` | `OBS_ID` header value |
| `mjd_obs` | `MJD_OBS` |
| `channel` | `LBCB` or `LBCR` |
| `filter` | `FILTER` |
| `exptime` | `EXPTIME` [s] |
| `role` | e.g. `science`, `flat`, `bias`, `heldout` (test subset) |
| `notes` | free text (e.g. why excluded) |

Selection criteria used (refer to the plan section, e.g.
`docs/planning/registration_migration_plan.md` §6.3.2 "Calibration data").

## Method

What `run.py` does, in a few lines, and the parameters that matter
(e.g. polynomial order, clipping, magnitude range).

## Results

Key numbers (fit rms, residual maps, comparisons with earlier products or
published values) and figures, with caveats.

## Intermediate data

Where the intermediate catalogs are (Zenodo DOI / release asset / this
directory) and their format.

## How to reproduce

`run.py` follows [`gain_rdnoise_example/run.py`](gain_rdnoise_example/run.py):
data path from `$LBCGO_RAW` or `--raw`, outputs and `run_log.json` written
next to it.

```
export LBCGO_RAW=/path/to/raw/frames
python calibration/<product_dir>/run.py
```
