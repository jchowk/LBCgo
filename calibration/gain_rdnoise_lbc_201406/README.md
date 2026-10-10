# LBC per-chip gain and read noise, 2014-06/07 (LBCB + LBCR)

## Summary

- **Product file(s):** `detector_rows.ecsv`, merged into
  `LBCgo/conf/lbc_detector.ecsv` (8 rows, `mjd_start` = 56830)
- **SHA-256:** `56e22dae4de28f87b5edf31a9538e1c9fcc84f612aab284339ddcc8d0528e4c7  detector_rows.ecsv`
  (after setting `mjd_start` by hand, see below; `run_log.json` records
  `ccdd292b…` for the file as `run.py` wrote it)
- **Made by / date:** J.C. Howk / 2026-10-09 (`run_log.json`)
- **LBCgo commit:** `e0c2e98e8d65ff61d19a8878c44703e66e45e318` (the code
  `run.py` imported, as in `run_log.json`)
- **Supersedes:** the 2010-03 rows (`../gain_rdnoise_lbc_201003/`) from
  MJD 56830 on. They remain in force for earlier dates. These rows are in
  turn superseded by `../gain_rdnoise_lbc_202505/` from MJD 60822.
- **Valid for:** LBCB and LBCR, all filters, MJD ≥ 56830 (2014-06-22 UT,
  the date of the earliest flats used). `mjd_end` is open; the 2025 rows
  take over at MJD 60822. Epochs between 2014 and 2025 are not measured.

## Purpose

Per-chip gain and read noise for the weight maps (`gain`, `rdnoise`), and
the flux gain (`gain_flux`) for Poisson errors of source fluxes, measured
with the photon-transfer method on raw twilight flats and biases from the
2014-06 (LBCB, LBCR) and 2014-07-01 (LBCR) calibration data. This is the
third measured epoch, after 2010-03 and before 2025-05.

## Inputs

`inputs_full.ecsv` lists every frame in the data directory (188 frames:
105 flats, 59 biases, 24 science frames). The flats and science frames
were already local; the biases were fetched from the LBT archive (59 files,
`archive.lbto.org`; the directory has no biases otherwise):

```
python calibration/make_inputs.py $LBCGO_RAW -o inputs_full.ecsv --overwrite --levels
```

`inputs.ecsv` is that list edited down to 11 sets (8 LBCB, 3 LBCR) of two
flats and two biases:

- **LBCB (sets 1–8):** unsaturated sky flats, 10–28k ADU, taken 2014-06-22
  (3 sets, g-SLOAN), 2014-06-24 (4 sets, g-SLOAN) and 2014-06-29 (1 set,
  r-SLOAN).
- **LBCR (sets 9–11):** 2014-07-01 r-SLOAN sky flats at 11k, 17k and 30k
  ADU.
- Each set is two consecutive flats of one OB (`lbcobnam`; same rotator
  angle), equal exposure time and ≤ 45 s apart; no frame is used in two
  sets. Exposure times are 1.2–6.2 s.
- Biases: consecutive pairs of the `25Bias_Bino` sequence of 2014-06-24
  (`propid` `biascheck`), skipping the first two frames of each channel; a
  different pair for each set. The same sequence supplies the biases for
  the LBCR flats of 2014-07-01 (7 days later) and for the LBCB flats of 06-22
  and 06-29.

Selection was iterative (the first two passes are not in the repository):

1. The first pass used 19 sets (10 LBCB, 9 LBCR), including three LBCB
   pairs of 3.28 s exposures from 2014-06-25. Their flat differences were
   strongly correlated spatially (`rho_sum` 1.5–2.3, against 0.03–0.06 for
   the other sets), so `gain_sum` read 0.5–0.7 against ~1.75 and the
   per-pixel gain was low. Replacement 1.28 s pairs from the same night
   (06-25) showed the same problem (`rho_sum` 0.4–0.45), so **all LBCB
   flats of 2014-06-25 were excluded**. The cause is unknown.
2. LBCR pairs with flats more than 50 s apart (sky level changing by 18–24 %
   between the two) gave low `gain_sum` (1.4–1.68 on chip 2, against 1.71–1.76
   for pairs ≤ 38 s apart) and were excluded. Two of the remaining LBCR
   pairs shared a flat and were reduced to one. That leaves three LBCR sets,
   the minimum for `summarize_gain_rdnoise` (`min_sets` = 3).

Saturated flats (> 40k ADU) were not used.

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
  sets whose two flats are ≤ 60 s apart (all sets here).

## Results

| channel, chip | sets | `gain` | `rdnoise` (e⁻) | `gain_flux` |
|---|---|---|---|---|
| LBCB 1 | 8 | 1.809 ± 0.010 | 8.75 | 1.796 ± 0.005 |
| LBCB 2 | 8 | 1.744 ± 0.008 | 9.78 | 1.749 ± 0.006 |
| LBCB 3 | 8 | 1.755 ± 0.009 | 9.07 | 1.756 ± 0.004 |
| LBCB 4 | 8 | 1.667 ± 0.010 | 9.43 | 1.660 ± 0.004 |
| LBCR 1 | 3 | 1.739 | 8.25 | 1.681 ± 0.006 |
| LBCR 2 | 3 | 1.758 | 8.09 | 1.759 ± 0.017 |
| LBCR 3 | 3 | 1.731 | 8.41 | 1.665 ± 0.006 |
| LBCR 4 | 3 | 1.778 | 7.85 | 1.778 ± 0.021 |

- The apparent gain rises with level, by +2.2–2.3 % per 10,000 ADU on LBCB
  and +3.6–4.9 % on LBCR. With the covariances summed to 3 px the slope is
  consistent with zero (LBCB −0.4 to +0.9 %, errors 0.5–0.75; LBCR −0.6 to
  −2.0 %, errors 0.3–0.8): the brighter-fatter effect, not non-linearity.
- **LBCB agrees with 2025** (`../gain_rdnoise_lbc_202505/`): `gain` +0.5 to
  +1.3 %, `gain_flux` +0.1 to +0.3 %, read noise +0.5 to +2.0 % (chip 4:
  +5 %), and the read noise in ADU is the same (e.g. chip 2: 5.6 ADU). Two
  products 11 years apart, from different flats and biases, agree to
  about 1 %.
- **LBCB disagrees with 2010** (`../gain_rdnoise_lbc_201003/`) on chips 2–4:
  `gain` −8 to −12 %, `gain_flux` −8 to −10 %, read noise −9 to −17 %
  (chip 1: `gain` +0.8 %, `gain_flux` +3 %, read noise −25 %). Either the
  detectors or electronics changed between 2010-03 and 2014-06, or the 2010
  measurement was biased. Not resolved here.
- **LBCR `gain_flux` agrees with 2025** within 1 % on every chip; `gain`
  agrees to within 0.2–3.2 %. Against 2010 `gain` is +2 to +9 %.
- **LBCR read noise is 7.9–8.4 e⁻ (4.4–4.9 ADU)**, 11–16 % below 2025
  (chip 2: 35 % below, since 2025 has 12.4 e⁻ there) and 5–16 % below 2010.
  The same method reproduces the LBCB read noise from 2025, so this does not
  look like a measurement artifact, but a change in the LBCR electronics
  between epochs has not been checked.
- Uncertainties: LBCB `gain0_err` (0.8–1.0 %) is statistical; the LBCB fit
  has 8 sets and scatter 0.6–0.8 %. **LBCR has three sets, so the straight
  line is nearly unconstrained (1 degree of freedom) and the quoted LBCR
  errors on `gain0` (≤ 0.010; two chips ≈ 0) are not meaningful.** The slope
  is set largely by the 30k-ADU set, which also has the largest
  correlation (`rho_sum` 0.18, against 0.05–0.08). The agreement of
  `gain_flux` with 2025 suggests the LBCR gains are good to 1–3 %, but this is
  a comparison, not an error estimate. Treat the LBCR rows as provisional
  until a deeper LBCR set is measured.
- Bias frames and flats are from different nights (see Inputs); the
  bias level and read noise are assumed constant over the 7 days.

## Intermediate data

None beyond the files here (`results_per_set.ecsv`, `gain_fit.ecsv`).

## How to reproduce

`run.py` follows [`../gain_rdnoise_example/run.py`](../gain_rdnoise_example/run.py):
data path from `$LBCGO_RAW` or `--raw`; outputs and `run_log.json` are
written next to it. The biases must be downloaded from the LBT archive first
(see `whatidid.md`).

```
export LBCGO_RAW=/Users/howk/Dropbox/Data/LBT/Raw/2014.06/LBC/
python calibration/gain_rdnoise_lbc_201406/run.py --box 0
```

`run.py` writes `mjd_start` = NaN (`PARAMS['mjd_start']`); 56830 was set
by hand in `detector_rows.ecsv` afterwards. Setting
`PARAMS['mjd_start'] = 56830.0` before a re-run reproduces the committed
file. See `whatidid.md` for the working notes.
