# LBCB per-chip gain and read noise, 2008-07 (LBCB only)

## Summary

- **Product file(s):** `detector_rows.ecsv`, merged into
  `LBCgo/conf/lbc_detector.ecsv` (4 LBCB rows, `mjd_start` = 54652)
- **SHA-256:** `99bf576dbbc58f349ecf0629438a383ae12e3bdd6a497b98a6dc091442d7bb82  detector_rows.ecsv`
  (after setting `mjd_start` by hand, see below; `run_log.json` records
  `d40d9b31…` for the file as `run.py` wrote it)
- **Made by / date:** J.C. Howk / 2026-10-10 (`run_log.json`)
- **LBCgo commit:** `18ed5ba0a851f206c6874a3f7a0189a0780048cb` (the code
  `run.py` imported, as in `run_log.json`)
- **Supersedes:** for LBCB, the Giallongo et al. (2008) rows
  (`../gain_rdnoise_lbcb_giallongo2008/`) from MJD 54652 on. They remain in
  force for earlier dates. These rows are in turn superseded by
  `../gain_rdnoise_lbc_200906/` from MJD 55003.
- **Valid for:** LBCB, all filters, MJD ≥ 54652 (2008-07-05 UT, the date of
  the flats used); the 2009-06 rows take over at MJD 55003. **No LBCR
  rows**: the LBCR data of this epoch could not be used (see Results).

## Purpose

Per-chip gain and read noise for the weight maps (`gain`, `rdnoise`), and
the flux gain (`gain_flux`) for Poisson errors of source fluxes, measured
with the photon-transfer method on raw twilight flats and biases of
2008-07-05 (UT). It is the earliest epoch measured with `LBCgo.detector`
and covers the period between the Giallongo et al. (2008) commissioning
values and the 2009-06 product.

## Inputs

Only part of the calibration data of 2008-07-05 was available (observing
log `20080705.log.txt`, Raw/2008.07/): archive tar files in
`Raw/2008.07/Calibrations/` (`IA2_LBT-SDT_12092008_*.tar`), extracted with
`tar -xf` into a scratch directory (`$LBCGO_RAW`), plus eight LBCR frames
(R-BESSEL) copied later as loose `.fits.gz` files. `inputs_full.ecsv`
lists the 40 frames used by `make_inputs.py`: 20 biases (10 LBCB, 10 LBCR)
and 20 flats:

```
python calibration/make_inputs.py $LBCGO_RAW -o inputs_full.ecsv --overwrite --levels
```

The 96 calibration frames of the log are not all present: no LBCB V- or
B-BESSEL flats, no LBCR flats of the `SkyFlatTest*` OBs (the five present
have a single extension, i.e. partial readouts, and were skipped by
`make_inputs.py`), and 15 of the 25 biases per channel are missing.

`inputs.ecsv` has 4 LBCB sets of two flats and two biases:

- **Flats:** the five `SkyFlat_Usr_rot180` SDT_Uspec flats (1.48 s,
  2008-07-05 11:52–11:55 UT, 25–43k ADU depending on chip, 40 s apart).
  With only five usable flats the four consecutive pairs **share frames**
  (only sets 1 and 3 are disjoint), so the sets are not independent. The
  levels are above the recommended 10–30k ADU, up to 0.66 × SATURATE.
- **Biases:** consecutive pairs of the `25Bias_Bino` sequence (`propid`
  `biascheck`, 06:29–06:34 UT), skipping the first two frames, a different
  pair for each set; taken about 5 h before the flats.

`inputs_lbcr_rejected.ecsv` and `results_lbcr_rejected.ecsv` hold the 7
LBCR sets that were tried (sets 5–8: i-SLOAN, 10.2 s; sets 9–11: R-BESSEL,
3.2 s; the sets of one filter share frames) and their per-set
measurements. LBCR r-SLOAN flats (42–65k ADU) are saturated or above the
0.7 × SATURATE limit of `measure_gain_rdnoise_files` and were not used.

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
  sets whose two flats are ≤ 60 s apart (all four).

## Results

| chip | `gain` | `rdnoise` (e⁻) | `gain_flux` |
|---|---|---|---|
| 1 | 1.792 ± 0.005 | 10.04 | 1.758 ± 0.004 |
| 2 | 1.957 ± 0.007 | 9.61 | 1.937 ± 0.003 |
| 3 | 1.921 ± 0.006 | 9.98 | 1.903 ± 0.004 |
| 4 | 1.811 ± 0.002 | 9.37 | 1.801 ± 0.003 |

- The apparent gain rises with level, by +2.0 to +2.4 % per 10,000 ADU.
  With the covariances summed to 3 px the slope is consistent with zero
  (−0.13 to +0.68 %, errors 0.3–0.5 %): the brighter-fatter effect.
  `rho_x` is +0.01 and `rho_sum` 0.06–0.11 on all chips.
- **Agreement with 2009-06** (`../gain_rdnoise_lbc_200906/`, independent
  flats and biases a year later): `gain` +1.0 to −0.4 %, `gain_flux` 0 to
  +1.4 %, read noise within 1.5 %. **Agreement with 2010-03**
  (`../gain_rdnoise_lbc_201003/`): `gain` within 0.8 %, `gain_flux` within
  1.0 %; the read noise is 9–14 % lower here and in 2009. The LBCB gains are
  therefore stable from 2008-07 to 2010-03.
- **Disagreement with 2014 and 2025** (`../gain_rdnoise_lbc_201406/`,
  `../gain_rdnoise_lbc_202505/`): on chips 2–4, `gain` is 9–14 % and
  `gain_flux` 8–11 % higher than in those epochs; chip 1 agrees to 0.4–2 %.
  Read noise agrees within 0.2–4 % on chips 2 and 4 and is 10 % (chip 3)
  and 15 % (chip 1) higher in 2008.
- **Against Giallongo et al. (2008, Table 1):** `gain` is 6–9 % lower and
  the read noise 12–17 % lower.
- **LBCR is not usable.** The flat differences have `rho_x` = −0.19 to
  −0.32 on chips 1, 3 and 4 (−0.03 to −0.05 on chip 2), in both the i-SLOAN
  and the R-BESSEL sets and independent of level, so it is not fringing.
  The per-pixel gains are 0.63–1.13 on chips 1, 3, 4 (1.65–1.73 on chip 2),
  inconsistent with any other epoch, and the fitted intercepts (0.79–1.65)
  are not credible. The bias differences do not show the effect
  (`bias_rho_x` ≥ −0.04), so it is signal-dependent. The same behaviour is
  present in 2009-06 and absent in 2010-03 and later (`rho_x` ≈ −0.01). The
  three usable R-BESSEL pairs (29–40k ADU) give `gain_sum` = 1.68, 1.61,
  1.64, 1.60 e⁻/ADU on chips 1–4, constant to ≲ 1 % over all levels and
  within ~8 % of the later `gain_flux` values; the i-SLOAN sets give only
  1.25–1.4, probably biased low. These flux gains are recorded here but not
  merged: the per-pixel `gain` and `rdnoise` columns cannot be filled
  credibly.
- Uncertainties: `gain0_err` (0.1–0.4 %) is statistical only, from four
  non-independent sets whose levels (25–43k ADU) lie far above zero, so the
  intercept is an extrapolation; the agreement with 2009-06 suggests ~1 %.
  A two-pair median (disjoint sets 1 and 3 only, the fallback for fewer
  than three sets) would give gains ~6 % higher (1.91, 2.09, 2.06, 1.94),
  because it samples the gain at ~30k ADU instead of at zero level.
- The biases come from 5 h before the flats (assumed stable).

## Intermediate data

None beyond the files here (`results_per_set.ecsv`, `gain_fit.ecsv`; the
rejected LBCR measurements are in `results_lbcr_rejected.ecsv`).

## How to reproduce

`run.py` follows [`../gain_rdnoise_example/run.py`](../gain_rdnoise_example/run.py):
data path from `$LBCGO_RAW` or `--raw`; outputs and `run_log.json` are
written next to it. Extract the tars first (see `whatidid.md`).

```
export LBCGO_RAW=<directory with the extracted .fits.gz files>
python calibration/gain_rdnoise_lbc_200807/run.py --box 0
```

`run.py` writes `mjd_start` = NaN (`PARAMS['mjd_start']`); 54652 was set
by hand in `detector_rows.ecsv` afterwards. Setting
`PARAMS['mjd_start'] = 54652.0` before a re-run reproduces the committed
file. The final run used the LBCB sets only: an earlier run with the LBCR
i-SLOAN sets included gave identical LBCB numbers, and its LBCR rows were
discarded (see Results). See `whatidid.md` for the working notes.
