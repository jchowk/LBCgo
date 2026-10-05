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
1000 px), and applies the photon-transfer formula. It then writes, next to
itself:

| file | content |
|---|---|
| `results_per_set.ecsv` | every measurement: one row per set and chip |
| `detector_rows.ecsv` | one row per channel and chip (median over sets), in the format of `LBCgo/conf/lbc_detector.ecsv`, with `source` naming this directory and the LBCgo commit |
| `run_log.json` | date, LBCgo version and commit (`+dirty` if the imported code has uncommitted changes), parameters, SHA-256 of every input and output |

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
flats and two biases per `set` (sets start as one per channel and UT date).

`inputs.ecsv` columns: `filename` (archive name), `obs_id`, `mjd_obs`,
`channel` (`LBCB`/`LBCR`), `filter`, `exptime`, `role` (`flat` or `bias`),
`set` (groups the four frames of one measurement), `notes`. `run.py` checks
that each set has exactly two flats and two biases and that the header
channel matches the `channel` column.

Flats should be at similar levels, roughly 10,000–30,000 ADU above bias and
well below saturation (`measure_gain_rdnoise_files` refuses flats above 70 %
of `SATURATE`). Several sets per channel at different epochs show whether the
values drift; if they do, write date-limited rows (`PARAMS['mjd_start']`,
`PARAMS['mjd_end']`, one product directory per epoch range) rather than one
open-ended median.

## Running

```
export LBCGO_RAW=/path/to/raw/frames
python calibration/<product_dir>/run.py          # or --raw /path
```

`tests/test_calibration_example.py` runs this script on synthetic frames with
known gain and read noise, so it keeps working as the package changes.
