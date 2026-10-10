# LBCB per-chip gain and read noise, 2009-06 (LBCB only)

## Summary

- **Product file(s):** `detector_rows.ecsv`, merged into
  `LBCgo/conf/lbc_detector.ecsv` (4 LBCB rows, `mjd_start` = 55003)
- **SHA-256:** `ee49fe992cad690bc9646e0ed0bc47a85913ee1d57c7cd8aa0097bc261d5d02a  detector_rows.ecsv`
  (after setting `mjd_start` by hand, see below; `run_log.json` records
  `e504c71e…` for the file as `run.py` wrote it)
- **Made by / date:** J.C. Howk / 2026-10-10 (`run_log.json`)
- **LBCgo commit:** `89091c38174a0d97ebd1a78c14b50ab7b3676140` (the code
  `run.py` imported, as in `run_log.json`)
- **Supersedes:** for LBCB, the Giallongo et al. (2008) rows
  (`../gain_rdnoise_lbcb_giallongo2008/`) from MJD 55003 on. They remain in
  force for earlier dates. These rows are in turn superseded by
  `../gain_rdnoise_lbc_201003/` from MJD 55273.
- **Valid for:** LBCB, all filters, MJD ≥ 55003 (2009-06-21 UT, the date of
  the flats used); the 2010-03 rows take over at MJD 55273. **No LBCR
  rows**: the LBCR data of this epoch could not be used (see Results).

## Purpose

Per-chip gain and read noise for the weight maps (`gain`, `rdnoise`), and
the flux gain (`gain_flux`) for Poisson errors of source fluxes, measured
with the photon-transfer method on raw twilight flats (2009-06-21 UT) and
biases (2009-06-22 UT). This is the earliest epoch measured with
`LBCgo.detector`; it covers the period between the Giallongo et al. (2008)
commissioning values and the 2010-03 product.

## Inputs

`inputs_full.ecsv` lists every frame in the data directory (155 frames: 67
flats, 50 biases, 38 science frames):

```
python calibration/make_inputs.py $LBCGO_RAW -o inputs_full.ecsv --overwrite --levels
```

`inputs.ecsv` is that list edited down to 10 LBCB sets of two flats and two
biases:

- **Flats (all 2009-06-21):** sets 1–6 are r-SLOAN sky flats (`SkyFlat_ri_5`,
  1.5–6.5 s) at about 10–17k ADU; sets 7–10 are g-SLOAN sky flats
  (`SkyFlat_gR_5`, 0.63–0.84 s) at about 22–33k ADU. Each set is two
  consecutive flats of one OB (equal exposure time) 33–40 s apart; no frame
  is used in two sets. The g-SLOAN flats at 64–65k ADU and the 60k frame at
  the start of the g sequence are saturated and were not used.
- **Biases (2009-06-22):** consecutive pairs of the `25Bias_Bino` sequence
  (`propid` `biascheck`), skipping the first two frames, a different pair
  for each set. They were taken about 22 h after the flats; the bias level
  and read noise are assumed constant over that interval.

`inputs_lbcr_rejected.ecsv` and `results_lbcr_rejected.ecsv` hold the 8 LBCR
sets (two z-SLOAN, six i-SLOAN) that were tried and rejected, with their
per-set measurements.

## Method

`run.py` (identical to `../gain_rdnoise_example/run.py`, see the commit
above), run on the whole trimmed chip (`--box 0`):

- `LBCgo.detector.measure_gain_rdnoise_files` per set: overscan-subtracted
  photon-transfer gain and read noise with variances in 50-px blocks, plus
  correlation diagnostics and `gain_sum` (covariances summed over lags
  ≤ 3 px).
- `LBCgo.detector.summarize_gain_rdnoise` per chip: `gain` = zero-level
  intercept of a straight-line fit of gain against level; `rdnoise` =
  `gain` × median read noise in ADU; `gain_flux` = median `gain_sum` of the
  sets whose two flats are ≤ 60 s apart (all ten sets).

## Results

| chip | `gain` | `rdnoise` (e⁻) | `gain_flux` |
|---|---|---|---|
| 1 | 1.789 ± 0.002 | 10.20 | 1.744 ± 0.006 |
| 2 | 1.938 ± 0.002 | 9.76 | 1.909 ± 0.008 |
| 3 | 1.930 ± 0.002 | 9.85 | 1.903 ± 0.006 |
| 4 | 1.802 ± 0.002 | 9.48 | 1.801 ± 0.006 |

- The apparent gain rises with level, by +1.9 to +2.4 % per 10,000 ADU.
  With the covariances summed to 3 px the slope is consistent with zero on
  chips 1, 3 and 4 (errors 0.4–0.55 % per 10k); chip 2 retains +1.4 ± 0.4 %.
  The trend is the brighter-fatter effect, with at most a small remainder
  on chip 2.
- **Agreement with 2010-03** (`../gain_rdnoise_lbc_201003/`): `gain` within
  −0.1 to −1.8 % and `gain_flux` within −2.1 to +0.2 %; the read noise is
  8.5–13 % lower (e⁻) here. The two epochs are nine months apart.
- **Disagreement with 2014 and 2025** (`../gain_rdnoise_lbc_201406/`,
  `../gain_rdnoise_lbc_202505/`): on chips 2–4, `gain` is 8–13 % and
  `gain_flux` 8–10 % higher than in those epochs; chip 1 agrees to about
  1–3 %. The read noise agrees within 0.2–6 % on chips 2 and 4, but is
  8–9 % (chip 3) and 17 % (chip 1) higher in 2009 than in 2014 and 2025.
  Together with the 2010 agreement, this puts the LBCB gain change between
  2010-03 (MJD 55273) and 2014-06 (MJD 56830).
- **Against Giallongo et al. (2008, Table 1):** `gain` is 6–9 % lower and
  the read noise 10–16 % lower.
- **LBCR is not usable.** Every LBCR set (both z-SLOAN and i-SLOAN) has
  `rho_x` = −0.27 to −0.42 in the flat difference, an rms about 2.7 times
  the shot-noise expectation (431 against 157 ADU, set 13 chip 1), and
  power in the difference rising towards the Nyquist frequency in x. The
  per-set gains (0.15–0.74) are not credible. The same chips show
  `rho_x` ≈ −0.01 in 2010 and ≈ +0.01 in 2025. The headers (geometry,
  amplifier sections) match those of LBCB and of later LBCR data, so the
  cause is not known; the biases look normal. These data should not be used
  for LBCR.
- Uncertainties: `gain0_err` (0.1 %) is statistical only; the choice of fit
  model alone moves the intercept by up to ~1 %. The `gain_flux` errors come
  from ten sets; the quoted `rms` of the fit is 0.2–0.3 %.

## Intermediate data

None beyond the files here (`results_per_set.ecsv`, `gain_fit.ecsv`; the
rejected LBCR measurements are in `results_lbcr_rejected.ecsv`).

## How to reproduce

`run.py` follows [`../gain_rdnoise_example/run.py`](../gain_rdnoise_example/run.py):
data path from `$LBCGO_RAW` or `--raw`; outputs and `run_log.json` are
written next to it.

```
export LBCGO_RAW=/Users/howk/Dropbox/Data/LBT/Raw/2009.06/
python calibration/gain_rdnoise_lbc_200906/run.py --box 0
```

`run.py` writes `mjd_start` = NaN (`PARAMS['mjd_start']`); 55003 was set
by hand in `detector_rows.ecsv` afterwards. Setting
`PARAMS['mjd_start'] = 55003.0` before a re-run reproduces the committed
file. The run used the LBCB sets only: the LBCR rows of an earlier run on
all 18 sets were discarded (see Results). See `whatidid.md` for the working
notes.
