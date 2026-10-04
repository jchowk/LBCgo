# Calibration provenance

This directory records **how** each calibration product used by LBCgo was
made, so that others can reproduce it or check it. It is not part of the
installed package (`pyproject.toml` packages only `LBCgo/`).

## Where things live

| What | Where |
|---|---|
| Code that derives a calibration (e.g. `LBCgo.detector.measure_gain_rdnoise_files`, the future distortion fit `LBCgo/register/calibrate.py`) | The `LBCgo` package, with tests in `tests/` |
| Products the pipeline reads at run time (`conf/lbc_detector.ecsv`, the planned `conf/distortion/*.fits`) | `LBCgo/conf/` |
| How each product was made: inputs, parameters, driver script, code version | **Here**, one subdirectory per product |
| Intermediate catalogs (e.g. matched star lists) that let someone re-derive a product without the raw images | A Zenodo deposit or GitHub release asset, linked from the product README; commit to the product directory only if small (≲ 5 MB compressed) |
| Raw frames and full intermediates | Not in git. Local disk and the LBT archive; recorded by identifier in `inputs.ecsv` |

Raw LBT data may be proprietary or partner-restricted. That is why the
intermediate catalogs are worth publishing: they let people without archive
access re-run the fit.

## One subdirectory per product

Name it after the product and its version, e.g. `distortion_lbcb_v1/`,
`gain_rdnoise_lbcr_2027a/`. A new version gets a new directory; old ones
stay, so earlier reductions remain traceable. Copy
[`NEW_PRODUCT.md`](NEW_PRODUCT.md) as the starting `README.md`.

Contents:

- `README.md`: what the product is, why it was made, the date, who ran it,
  the LBCgo commit (`git rev-parse HEAD`), the product file(s) it wrote and
  their SHA-256 checksums, and the main results (e.g. fit rms) with any
  caveats.
- `inputs.ecsv`: one row per input frame (columns in `NEW_PRODUCT.md`).
  Archive filenames and `OBS_ID`s, not local paths.
- `run.py`: the driver. It only orchestrates calls into the `LBCgo`
  package and sets parameters; any real logic belongs in the package, with
  tests. Local data paths come from an environment variable or command-line
  argument, never hard-coded.
- Optional: small diagnostic figures, small intermediate catalogs.

The product file itself should point back here: e.g. the `source` column of
`conf/lbc_detector.ecsv`, or a `PROVENANCE` header keyword in a distortion
model file naming this directory.

## Products

| Directory | Product | Status |
|---|---|---|
| [`gain_rdnoise_lbcb_giallongo2008/`](gain_rdnoise_lbcb_giallongo2008/README.md) | LBCB rows of `conf/lbc_detector.ecsv` | Published values, seeded 2026-10-04 |
