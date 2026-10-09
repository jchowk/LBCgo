# LBC per-chip gain and read noise, 2010-03-18/19 (LBCB + LBCR)

## Summary

- **Product file(s):** `detector_rows.ecsv`, merged into
  `LBCgo/conf/lbc_detector.ecsv` (8 rows, `mjd_start` = 55273, `mjd_end` open)
- **SHA-256:** `c251e5394fa222e8565b198f47a1d6e2ff1a8b443157946d256e8b54b743ca96  detector_rows.ecsv`
  (as `run.py` wrote it; `mjd_start` comes from `PARAMS`, no hand edit)
- **Made by / date:** J.C. Howk / 2026-10-09 (`run_log.json`)
- **LBCgo commit:** `91f3fd970dd32f3e6d6b2bc17a36a65d064c32f5+dirty` (as in
  `run_log.json`). The package directory then differed from the commit only
  in `LBCgo/00ToDo.md`. The lazy header scan in `calibration/make_inputs.py`
  (see Inputs) postdates this commit; it affects only how `inputs_full.ecsv`
  is built, not the measurement.
- **Supersedes:** nothing for LBCR (header values before); for LBCB, the
  Giallongo et al. (2008) rows (`../gain_rdnoise_lbcb_giallongo2008/`) from
  MJD 55273 on. They remain in force for earlier dates. These rows are in
  turn superseded by `../gain_rdnoise_lbc_202505/` from MJD 60822.
- **Valid for:** LBCB and LBCR, MJD ≥ 55273 (2010-03-18 UT). `mjd_end` is
  open: until further epochs are measured, the rows are applied up to MJD
  60822, on the assumption that the detectors did not change (see Results
  for how much the 2010 and 2025 values differ).

## Purpose

Per-chip gain and read noise for the weight maps (`gain`, `rdnoise`), and
the flux gain (`gain_flux`) for Poisson errors of source fluxes, measured
with the photon-transfer method on raw twilight flats and biases from the
nights of 2010-03-18 and 2010-03-19 (UT). This is the earliest epoch
measured by LBCgo for both channels; LBCB also has the published values of
Giallongo et al. (2008, 2006 commissioning data).

## Inputs

Raw frames: `/Users/howk/Dropbox/Data/LBT/Raw/2010.03/Calibration/`
(`$LBCGO_RAW`).

`inputs_full.ecsv` lists every raw frame in the directory that
`make_inputs.py` selects (83 frames):

```
python calibration/make_inputs.py $LBCGO_RAW -o inputs_full.ecsv --overwrite --levels
```

Twenty of the frames (the 10 LBCB and 10 LBCR biases of 2010-03-18
22:36–22:40 UT) have the invalid header card `PA_PNT = nan`, which
`ccdproc.ImageFileCollection` cannot parse (`VerifyError: Unparsable card
(PA_PNT)`). `make_inputs.py` was changed to read only the needed keywords
with `fits.getheader`; its output for the 2025.05 data is unchanged. The same
failure is expected in `lbcproc.py` (item in `LBCgo/00ToDo.md`).

`inputs.ecsv` is that list edited down to 12 sets (6 per channel) of two
flats and two biases. Times are UT on 2010-MMDD; flat levels are the chip-2
overscan-subtracted medians from `--levels`.

| set | ch | filter | exptime (s) | flats (MMDD hh:mm:ss, ADU) | flat dt (s) | biases (MMDD hh:mm) |
|---|---|---|---|---|---|---|
| 1 | LBCB | SDT_Uspec | 0.623 | 0318 13:02:48 (6617), 0318 13:03:20 (7577) | 32 | 0318 22:36, 0318 22:36 |
| 2 | LBCB | SDT_Uspec | 0.827 | 0318 01:52:54 (9159), 0318 01:53:28 (7923) | 35 | 0318 01:27, 0318 01:27 |
| 3 | LBCB | SDT_Uspec | 0.623 | 0318 13:03:54 (8754), 0318 13:04:28 (10065) | 34 | 0318 22:37, 0318 22:37 |
| 4 | LBCB | SDT_Uspec | 0.827 | 0318 01:51:45 (12257), 0318 01:52:21 (10572) | 35 | 0318 01:28, 0318 01:28 |
| 5 | LBCB | SDT_Uspec | 0.523 | 0318 13:06:40 (13396), 0318 13:07:14 (15247) | 33 | 0318 22:38, 0318 22:38 |
| 6 | LBCB | SDT_Uspec | 0.523 | 0318 13:07:47 (17444), 0318 13:08:21 (19805) | 34 | 0318 22:39, 0318 22:39 |
| 7 | LBCR | I-BESSEL | 2.232 | 0319 02:02:55 (14401), 0319 02:03:25 (12809) | 30 | 0318 22:36, 0318 22:37 |
| 8 | LBCR | I-BESSEL | 2.231 | 0319 02:01:54 (18223), 0319 02:02:25 (16143) | 30 | 0318 22:37, 0318 22:37 |
| 9 | LBCR | V-BESSEL | 0.531 | 0319 01:53:04 (19313), 0319 01:53:32 (16931) | 29 | 0318 22:38, 0318 22:38 |
| 10 | LBCR | I-BESSEL | 7.236 | 0318 02:06:39 (26505), 0318 02:07:18 (23134) | 39 | 0318 22:38, 0318 22:39 |
| 11 | LBCR | I-BESSEL | 5.228 | 0319 12:50:40 (22103), 0319 12:51:15 (25368) | 35 | 0318 22:39, 0318 22:39 |
| 12 | LBCR | V-BESSEL | 0.531 | 0319 01:52:07 (24878), 0319 01:52:36 (21745) | 29 | 0318 22:39, 0318 22:40 |

Selection:

- Flats are consecutive frames of one sequence with the same exposure time,
  each pair used once, at levels below 30,000 ADU (the 10,000–30,000 ADU the
  docstring of `measure_gain_rdnoise_files` recommends; the lowest LBCB sets
  reach 7,000 ADU).
- The six LBCB B-BESSEL flats (32,000–60,000 ADU) and the LBCR flats above
  30,000 ADU are excluded. Hence LBCB is measured through `SDT_Uspec` only.
- Biases are consecutive frames from the bias block nearest in time. The
  first four LBCR biases of 2010-03-18 01:27 (732, 413, 251 and 32 ADU
  after overscan subtraction, against about 20 ADU for the rest of the
  block) show a decaying level, probably a readout transient, and are not
  used.
- **Set 10 biases:** the settled 01:27 LBCR biases, which were the nearest in
  time to set 10, give a bias-difference noise of 6.6–8.1 ADU (chips 1–4),
  against 5.3–6.0 ADU in the 22:36 block, and set 10 then gave a read
  noise of 12–13.6 e⁻ against 9–10.5 e⁻ in the other LBCR sets. The overscan
  rms is the same in both blocks (5.4–5.9 ADU) and in the flats (5.7–6.5 ADU),
  so this is a bias pattern drifting in the image area, not read noise. Set
  10 therefore uses two 22:36 frames, which are also used by sets 9 and 11.
  The LBCB sets 2 and 4 use the 01:27 LBCB biases (nothing anomalous in their
  read noise).

## Method

`run.py` (identical to `../gain_rdnoise_lbc_202505/run.py` except for
`PARAMS['mjd_start'] = 55273.0`), run on the whole trimmed chip
(`--box 0`):

- `LBCgo.detector.measure_gain_rdnoise_files` per set: overscan-subtracted
  photon-transfer gain and read noise with variances in 50-px blocks, plus
  correlation diagnostics and `gain_sum` (covariances summed over lags
  ≤ 3 px).
- `LBCgo.detector.summarize_gain_rdnoise` per chip: `gain` = zero-level
  intercept of a straight-line fit of gain against level; `rdnoise` =
  `gain` × median read noise in ADU; `gain_flux` = median `gain_sum` of the
  sets whose two flats are ≤ 60 s apart (all 12 sets).

## Results

Fit per channel and chip (6 sets each; from `gain_fit.ecsv`):

| channel | chip | level range (ADU) | gain0 (e⁻/ADU) | slope (%/10k ADU) | gain_sum slope (%/10k ADU) | rdnoise (e⁻) | gain_flux (e⁻/ADU) |
|---|---|---|---|---|---|---|---|
| LBCB | 1 | 7871–20652 | 1.795 ± 0.006 | +1.81 ± 0.24 | −0.02 ± 0.66 | 11.72 | 1.741 ± 0.004 |
| LBCB | 2 | 6988–18348 | 1.973 ± 0.026 | +2.65 ± 1.10 | +0.53 ± 0.74 | 11.04 | 1.950 ± 0.005 |
| LBCB | 3 | 7056–18532 | 1.932 ± 0.004 | +2.04 ± 0.17 | +0.19 ± 0.52 | 10.94 | 1.900 ± 0.003 |
| LBCB | 4 | 7312–19190 | 1.816 ± 0.002 | +2.22 ± 0.10 | +0.14 ± 0.45 | 10.36 | 1.798 ± 0.003 |
| LBCR | 1 | 13334–24194 | 1.637 ± 0.030 | +3.96 ± 0.90 | −1.36 ± 3.14 | 9.75 | 1.633 ± 0.019 |
| LBCR | 2 | 13582–24781 | 1.705 ± 0.033 | +4.40 ± 0.95 | +0.53 ± 3.34 | 9.59 | 1.696 ± 0.020 |
| LBCR | 3 | 14721–27206 | 1.594 ± 0.053 | +1.99 ± 1.53 | +2.98 ± 2.90 | 8.86 | 1.562 ± 0.018 |
| LBCR | 4 | 12690–22846 | 1.743 ± 0.040 | +4.32 ± 1.21 | −2.32 ± 4.12 | 9.26 | 1.740 ± 0.026 |

- The apparent gain rises with level on both channels (+1.8 to +2.7 % per
  10,000 ADU on LBCB, +2.0 to +4.4 % on LBCR), as in 2025. With the
  covariances summed to 3 px, the LBCB slopes are consistent with zero
  (−0.02 to +0.53 %, errors 0.45–0.74), the brighter-fatter signature. On
  LBCR the `gain_sum` slopes have 3–4 % errors and cannot separate
  brighter-fatter from non-linearity.
- `gain_flux` is below `gain` by 1–3 % on LBCB and by 0.2–2 % on LBCR.
- **Against Giallongo et al. (2008, Table 1), LBCB:** gain 0.92, 0.94, 0.94,
  0.92 × the published values (`gain_flux`: 0.89–0.93 ×); read noise 1.03,
  0.95, 0.94, 0.93 × (11.4, 11.6, 11.6, 11.2 e⁻ published). Their flat
  levels are not recorded here: at 14,000–16,000 ADU (set 5) the per-set
  gains are 0.94–0.98 × theirs, so part of the offset may be the level at
  which they measured.
- **Against 2025 (`../gain_rdnoise_lbc_202505/`):** LBCB gains are 1.00,
  1.15, 1.11 and 1.10 × the 2025 values, and the read noise is 15–35 %
  higher in 2010 (8.7–9.6 e⁻ in 2025). The 2010 chip-to-chip gain pattern
  (chip 2 highest, chip 1 lowest) matches Giallongo's, not 2025's. LBCR
  gains agree to within about 6 % (0.94, 1.00, 0.94, 1.00 ×), and the read
  noise to within about 6 % except chip 2 (9.6 e⁻ here, 12.4 e⁻ in 2025). Applying these
  rows up to MJD 60822 is therefore better supported for LBCR than for
  LBCB.
- **Not resolved:** four of the six LBCR sets are I-band, where fringing in
  twilight flats may bias the photon-transfer variances (see the docstring of
  `measure_gain_rdnoise_files`). At about 23,500 ADU the V-band set 12 reads
  2 % higher on chip 1 than the I-band set 11, which could be fringing or
  noise. The LBCR flats span only 13,000–27,000 ADU, so the zero-level
  intercept is an extrapolation with 2–3 % errors; LBCB covers 7,000–20,000
  ADU.
- Uncertainties: `gain0_err` is statistical only; the choice of fit model
  alone moves the intercept by up to ~1 % (the 2025 estimate, not repeated
  here). The `gain_flux` errors come from six pairs per chip.

## Intermediate data

None beyond the files here (`results_per_set.ecsv`, `gain_fit.ecsv`).
`inputs_full.ecsv` is the unedited frame list.

## How to reproduce

`run.py` follows [`../gain_rdnoise_example/run.py`](../gain_rdnoise_example/run.py):
data path from `$LBCGO_RAW` or `--raw`; outputs and `run_log.json` are
written next to it. It takes about 7 minutes.

```
export LBCGO_RAW=/Users/howk/Dropbox/Data/LBT/Raw/2010.03/Calibration/
python calibration/make_inputs.py $LBCGO_RAW -o inputs_full.ecsv --overwrite --levels
python calibration/gain_rdnoise_lbc_201003/run.py --box 0
```

Run these in the `lbcgo` environment. `inputs.ecsv` is edited by hand from
`inputs_full.ecsv` as described under Inputs, and is checked in. Its SHA-256
is `35bf69495fba069acbfbf031897b00931976c7b079dfbde88f8de573249a52e2`.
Merging into `LBCgo/conf/lbc_detector.ecsv` is a separate step, covered by
`tests/test_detector.py::test_packaged_table_has_2010_rows_from_product`.
See `whatidid.md` for the working notes.
