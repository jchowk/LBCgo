# LBC per-chip gain and read noise, 2025-05-27 (LBCB + LBCR)

## Summary

- **Product file(s):** `detector_rows.ecsv`, merged into
  `LBCgo/conf/lbc_detector.ecsv` (8 rows, `mjd_start` = 60822)
- **SHA-256:** `5b02e5d97426c4d5ee043bbeb47793c9faeb0a4b1dd848b5e936e964e2730d4f  detector_rows.ecsv`
  (after setting `mjd_start` by hand, see below; `run_log.json` records
  `00f1f5d8…` for the file as `run.py` wrote it)
- **Made by / date:** J.C. Howk / 2026-10-06 (`run_log.json`)
- **LBCgo commit:** `e3f256c8e5b13e5ade06abc5c403128a1ed56ccc` (the code
  `run.py` imported, as in `run_log.json`; the outputs were committed in
  `edea303`)
- **Supersedes:** nothing for LBCR (header values before); for LBCB, the
  Giallongo et al. (2008) rows (`../gain_rdnoise_lbcb_giallongo2008/`) from
  MJD 60822 on. They remain in force for earlier dates.
- **Valid for:** LBCB and LBCR, all filters, MJD ≥ 60822 (open-ended until
  later epochs are measured)

## Purpose

Per-chip gain and read noise for the weight maps (`gain`, `rdnoise`), and
the flux gain (`gain_flux`) for Poisson errors of source fluxes, measured
with the photon-transfer method on raw twilight flats and biases from
2025-05-27 (UT).

## Inputs

`inputs_full.ecsv` lists every frame in the data directory:

```
python calibration/make_inputs.py $LBCGO_RAW -o inputs_full.ecsv --overwrite --levels
```

`inputs.ecsv` is that list edited down to 10 sets (5 per channel) of two
flats and two biases: unsaturated sky flats (`SkyFlat_BR_pa0`/`_pa180`,
B-BESSEL on LBCB and R-BESSEL on LBCR) spanning a wide range of levels
(about 5,000–27,000 ADU), plus bias pairs. `biascheck` is the PROPID of
every bias and does not mark a special subset.

## Method

`run.py` (identical to `../gain_rdnoise_example/run.py` at the commit
above), run on the whole trimmed chip (`--box 0`):

- `LBCgo.detector.measure_gain_rdnoise_files` per set: overscan-subtracted
  photon-transfer gain and read noise with variances in 50-px blocks, plus
  correlation diagnostics and `gain_sum` (covariances summed over lags
  ≤ 3 px).
- `LBCgo.detector.summarize_gain_rdnoise` per chip: `gain` = zero-level
  intercept of a straight-line fit of gain against level; `rdnoise` =
  `gain` × median read noise in ADU; `gain_flux` = median `gain_sum` of the
  sets whose two flats are ≤ 60 s apart (sets 1, 7, 8 on LBCB; 4, 9, 10 on
  LBCR; set 2 is 60.5 s apart and excluded).

## Results

- The apparent gain rises with level, by +1.9–2.0 % per 10,000 ADU on LBCB
  and +3.2–5.0 % on LBCR. With the covariances summed to 3 px the slope is
  consistent with zero on every chip (−0.78 to +0.05 % per 10,000 ADU,
  errors 0.25–0.93). The trend is the brighter-fatter effect on both
  channels (on LBCR with covariances beyond the nearest neighbours), not
  non-linearity (LBCR: ≲ 0.7 % at 22,000 ADU).
- `gain` and `gain_flux` agree within 1.3 % except on LBCR chip 1 (−4.4 %)
  and chip 2 (+3.4 %), where pixels are correlated along the readout
  direction even at zero signal.
- LBCB gains are 0.82–0.92 × the values of Giallongo et al. (2008, Table 1),
  while the read noise in ADU agrees for chip 2: the two electron scales
  differ. Not yet explained.
- LBCR chip 2 read noise is 12.4 e⁻ (7.3 ADU), against 9.3–9.8 e⁻ for the
  other LBCR chips.
- Uncertainties: `gain0_err` (0.1–0.3 %) is statistical only; the choice of
  fit model alone moves the intercept by up to ~1 %. The `gain_flux` errors
  come from three pairs per chip and are likely optimistic; pairs taken
  140–210 s apart read 0.5–0.8 % low.

## Intermediate data

None beyond the files here (`results_per_set.ecsv`, `gain_fit.ecsv`).

## How to reproduce

`run.py` follows [`../gain_rdnoise_example/run.py`](../gain_rdnoise_example/run.py):
data path from `$LBCGO_RAW` or `--raw`; outputs and `run_log.json` are
written next to it.

```
export LBCGO_RAW=/Users/howk/Dropbox/Data/LBT/Raw/2025.05_calib/
python calibration/gain_rdnoise_lbc_202505/run.py --box 0
```

`run.py` writes `mjd_start` = NaN (`PARAMS['mjd_start']`); 60822 was set
by hand in `detector_rows.ecsv` afterwards. Setting
`PARAMS['mjd_start'] = 60822.0` before a re-run reproduces the committed
file. See `whatidid.md` for the working notes.
