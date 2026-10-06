# Example: measuring per-chip gain and read noise

**This is a worked example of a calibration driver, not a product.** The rows
in `inputs.ecsv` are illustrative, not real archive frames. To make a real
product, copy this directory to a new name (e.g. `gain_rdnoise_2027a/`),
replace `inputs.ecsv` with the real frames, rewrite this README from
[`../NEW_PRODUCT.md`](../NEW_PRODUCT.md), and run it.

## What `run.py` does

For each `set` in `inputs.ecsv` (two raw flats + two raw biases of one channel
and epoch), it calls `LBCgo.detector.measure_gain_rdnoise_files`, which
overscan-subtracts and trims each chip, uses a central region (`box`, default
1000 px), and applies the photon-transfer formula with the variances measured
in 50-px blocks (`cell`). Blocks make the result insensitive to small
illumination differences between the two flats, which otherwise bias the
gain low (−6 % to −12 % in the first real run, 2025-05). Each measurement
also records the flat level, the read noise in ADU, and the nearest-neighbour
correlations of the flat difference.

The apparent gain rises with flat level (brighter-fatter effect: +2 % per
10,000 ADU on LBCB, +3–5 % on LBCR in the 2025-05 data), so a median over
sets depends on which levels were observed. `LBCgo.detector.summarize_gain_rdnoise`
therefore fits gain = g0 + slope × level per chip and adopts the intercept
g0; read noise = g0 × median read noise in ADU. With fewer than
`PARAMS['min_sets']` (3) sets for a chip it falls back to the median;
`PARAMS['gain_model'] = 'median'` forces that. It writes, next to itself:

| file | content |
|---|---|
| `results_per_set.ecsv` | every measurement: one row per set and chip, with `level` [ADU], `rdnoise_adu`, `k`, `rho_x`, `rho_y`, `bias_rho_x`, `bias_rho_y`, `gain_nn`, `rho_sum`, `rho_sum_err`, `bias_rho_sum`, `gain_sum`, `gain_sum_err` |
| `gain_fit.ecsv` | one row per channel and chip: `model` (`linear`/`median`), `n`, `gain0` ± `gain0_err`, slope in % per 10,000 ADU ± error, residual `rms`, level range, median gain for comparison, read noise, and the brighter-fatter diagnostics `rho_slope_per_10k`, `gain_nn_slope_pct_per_10k`, `rho_sum_slope_per_10k`, `gain_sum_median`, `gain_sum_slope_pct_per_10k` ± `gain_sum_slope_pct_err` |
| `detector_rows.ecsv` | one row per channel and chip (gain = `gain0`), in the format of `LBCgo/conf/lbc_detector.ecsv`, with `source` naming this directory and the LBCgo commit |
| `run_log.json` | date, LBCgo version and commit (`+dirty` if the imported code has uncommitted changes), parameters, SHA-256 of every input and output |

Reading `gain_fit.ecsv`: if the brighter-fatter effect explains the slope,
`rho_slope_per_10k` is positive (correlations grow with level) and
`gain_nn_slope_pct_per_10k` (gain with the nearest-neighbour covariances
added back) is close to zero. `gain_nn` ignores longer-range covariances,
so it is a diagnostic, not the adopted gain.

`gain_sum` adds back the covariances of all lags up to `PARAMS['max_lag']`
(3 px). The brighter-fatter effect only moves charge between pixels, so it
leaves `gain_sum` flat; non-linearity creates no covariances and gives
`gain_sum` the same slope as the gain. Compare
`gain_sum_slope_pct_per_10k` ± `gain_sum_slope_pct_err` with
`slope_pct_per_10k`. Summing 48 lags is noisy (~1 % per set for the default
1000-px box), so for this test run on the whole chip: `run.py --box 0`.

It does **not** edit `LBCgo/conf/lbc_detector.ecsv`. Review
`detector_rows.ecsv` (e.g. against Giallongo et al. 2008 Table 1 for LBCB),
then merge the rows into the package table as a separate, reviewed change,
replacing or date-limiting older rows. For a real product, commit the three
output files with the directory; for this example, don't.

## Inputs

Build it with
`python calibration/make_inputs.py <raw_dir> --channel LBCR --imagetyp flat zero --levels -o inputs.ecsv`:
`--levels` adds each flat's overscan-subtracted level (`median_adu`, centre of
chip 2), which helps pick two flats at similar levels. Then keep exactly two
flats and two biases per `set` (sets start as one per channel and UT date),
with sets at several levels.

`inputs.ecsv` columns: `filename` (archive name), `obs_id`, `mjd_obs`,
`channel` (`LBCB`/`LBCR`), `filter`, `exptime`, `role` (`flat` or `bias`),
`set` (groups the four frames of one measurement), `notes`. `run.py` checks
that each set has exactly two flats and two biases and that the header
channel matches the `channel` column.

Choose each flat pair from **consecutive exposures of one sequence**: same
filter, same rotator angle (`LBCOBNAM` `..._pa0` vs `..._pa180` flips the
twilight gradient), similar exposure times, minutes apart at most. Avoid
z- and Y-band twilight flats: fringing changes through twilight on scales
the 50-px blocks do not remove. Use several pairs per channel at different
levels: the fit needs the spread, and its residuals show the uncertainty.

The two flats of a pair should be at similar levels; across pairs, spread
the levels over roughly 5,000–35,000 ADU above bias (at least three pairs
per chip, for the fit), well below saturation (`measure_gain_rdnoise_files` refuses flats above 70 %
of `SATURATE`). Several sets per channel at different epochs show whether the
values drift; if they do, write date-limited rows (`PARAMS['mjd_start']`,
`PARAMS['mjd_end']`, one product directory per epoch range) rather than one
open-ended median.

## Running

```
export LBCGO_RAW=/path/to/raw/frames
python calibration/<product_dir>/run.py          # or --raw /path
python calibration/<product_dir>/run.py --box 0  # whole chip
```

An existing product directory keeps its own copy of `run.py`; copy the
current example over it (keeping the product's `PARAMS`) to get new
outputs such as `gain_sum`.

`tests/test_calibration_example.py` runs this script on synthetic frames with
known gain and read noise, so it keeps working as the package changes.
