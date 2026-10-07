# <Product name and version, e.g. "LBCB static distortion model, v1">

Copy this file to `calibration/<product_dir>/README.md` and fill it in.
See [`README.md`](README.md) for the convention.

## Summary

- **Product file(s):** e.g. `LBCgo/conf/lbc_detector.ecsv`
- **SHA-256:** output of `00f1f5d883102fe009b55ede77660a96d0bfca2051b79ac64e52df627fb8792b  detector_rows.ecsv`
- **Made by / date:** `J.C. Howk / 2026.10.05`
- **LBCgo commit:** `edea30310147365df9ba4351b65b3baad93cdf59`

- **Valid for:** `LBCB+LBCR`, `MJD 60822–???`

## Purpose

Why this product was made and what it is used for in the pipeline.

## Inputs

`inputs.ecsv` contains input file information created using 
```
python ~/python/LBCgo/calibration/make_inputs.py ./ -o inputs.ecsv --overwrite --levels
```

The resulting file was edited to downselect to 6 "sets" of calibrations. The largest cut is on identifying unsaturated sky flats and spanning a wide range of flat intensities.

## Method

Uses `LBCgo.detector.measure_gain_rdnoise_files` to measure the gain and readnoise with procedures in place as of 2026.10.05. 

## Results

The gain results are based on an assumed "brighter-fatter" behavior, and the flats in this run confirm that interpretation for both sides. The gains for LBCB are lower than the values quoted in Giallongo et al. (2008). 

## Intermediate data

N/A (yet)

## How to reproduce

`run.py` follows [`gain_rdnoise_example/run.py`](gain_rdnoise_example/run.py):
data path from `$LBCGO_RAW` or `--raw`, outputs and `run_log.json` written
next to it.

```
export LBCGO_RAW=/Users/howk/Dropbox/Data/LBT/Raw/2025.05_calib/
python calibration/gain_rdnoise_lbc_202505/run.py
```
