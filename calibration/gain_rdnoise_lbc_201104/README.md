# LBCB per-chip gain and read noise, 2011-04 (LBCB only)

## Summary

- **Product file(s):** `detector_rows.ecsv`, merged into
  `LBCgo/conf/lbc_detector.ecsv` (4 LBCB rows, `mjd_start` = 55677)
- **SHA-256:** `9fb87760f7239b9aaf72b31993bd934773ff535cebd661cc9df00b4771734ba1  detector_rows.ecsv`
  (as written by `run.py`; `PARAMS['mjd_start']` = 55677.0 is set in `run.py`,
  no hand edit)
- **Made by / date:** J.C. Howk / 2026-10-11 (`run_log.json`)
- **LBCgo commit:** `edcbd0d32294c15794808acdbe5c6e83a6b3bc7f` (the code
  `run.py` imported, as in `run_log.json`; LBCgo 0.2.0.dev0)
- **Supersedes:** for LBCB, the 2010-03 rows (`../gain_rdnoise_lbc_201003/`)
  from MJD 55677 on; those remain in force from MJD 55273. These rows are in
  turn superseded by `../gain_rdnoise_lbc_201406/` from MJD 56830.
- **Valid for:** LBCB, all filters, MJD ≥ 55677 (2011-04-26 UT, the date of
  the first flats used); the 2014-06 rows take over at MJD 56830. **No LBCR
  rows**: the LBCR data of this epoch were not usable (see Results); LBCR
  stays on the 2010-03 rows.

## Purpose

Per-chip gain and read noise for the weight maps (`gain`, `rdnoise`), and
the flux gain (`gain_flux`) for Poisson errors of source fluxes, measured
with the photon-transfer method on raw twilight flats (2011-04-26, -28 and
-29 UT) and biases (2011-04-27 and -28 UT). It adds an LBCB epoch between
the 2010-03 product and the 2014-06 product, which differ in gain on chips
2–4 by about 10 %, and so brackets when that change did not yet happen.

## Inputs

`inputs_full.ecsv` lists every frame in the data directory (279 frames: 100
biases, 40 flats, 87 science and 52 unclassified science frames; 8 default
sets):

```
python calibration/make_inputs.py $LBCGO_RAW -o inputs_full.ecsv --overwrite --levels
```

The biases were not in the original download of the directory and were
fetched afterwards from the LBT archive; the file was regenerated after that.

`inputs.ecsv` is that list edited down to 5 LBCB sets of two flats and two
biases:

- **Flats:** set 1, B-BESSEL, 2011-04-26 (12.8 s, ~11k ADU); sets 2 and 3,
  B-BESSEL, 2011-04-28 (15.3 s, ~20k and ~15k ADU); sets 4 and 5, U_spec
  (`SDT_Uspec`), 2011-04-29 (27.3 s, ~30k and ~20k ADU). Each set is two
  consecutive flats of one OB (equal exposure time) 45–61 s apart (sets 4
  and 5: 60.5 s, so they are outside the 60 s cut for `gain_flux`); no
  frame is used in two sets. The remaining ramp frames (the 42k ADU first U_spec flat and the fainter B-BESSEL frames) were not
  used. The flats are twilight ramps, so the level changes by 10–15 % between the two
  frames of a set.
- **Biases:** consecutive pairs of the `25Bias_Bino` sequences (`propid`
  `biascheck`) of 2011-04-27 (set 1) and 2011-04-28 (sets 2–5), skipping
  the first two frames of each block, a different pair for each set. They
  were taken 2 h (sets 2, 3) to 22 h (set 1, biases after the flats) and
  26 h (sets 4, 5, biases before the flats) from the flats; the bias level
  and read noise are assumed constant over that interval.

`inputs_lbcr_rejected.ecsv` and `results_lbcr_rejected.ecsv` hold the 5 LBCR
sets (I-BESSEL sets 6–8, V-BESSEL sets 9–10) that were tried and rejected,
with their per-set measurements (chips 1, 2, 4).

## Method

`run.py` (the 2008-07 driver, with `PARAMS['chips']` added to skip the blank
LBCR chip 3; no LBCB-relevant change), run on the whole trimmed chip
(`--box 0`):

- `LBCgo.detector.measure_gain_rdnoise_files` per set: overscan-subtracted
  photon-transfer gain and read noise with variances in 50-px blocks, plus
  correlation diagnostics and `gain_sum` (covariances summed over lags
  ≤ 3 px).
- `LBCgo.detector.summarize_gain_rdnoise` per chip: `gain` = zero-level
  intercept of a straight-line fit of gain against level; `rdnoise` =
  `gain` × median read noise in ADU; `gain_flux` = median `gain_sum` of the
  sets whose two flats are ≤ 60 s apart (sets 1–3).

## Results

| chip | `gain` | `rdnoise` (e⁻) | `gain_flux` |
|---|---|---|---|
| 1 | 1.793 ± 0.004 | 10.32 | 1.721 ± 0.004 |
| 2 | 1.951 ± 0.002 | 9.82 | 1.914 ± 0.005 |
| 3 | 1.928 ± 0.003 | 10.20 | 1.885 ± 0.006 |
| 4 | 1.822 ± 0.006 | 9.61 | 1.802 ± 0.004 |

- The apparent gain rises with level, by +1.5 to +2.1 % per 10,000 ADU
  (rms about the fit 0.1–0.4 %). With the covariances summed to 3 px the
  slope is −0.8 to +0.4 % per 10k (errors 0.2–0.3 %): consistent with zero
  on chip 2, and slightly negative on chips 1, 3 and 4 (−0.5 to −0.8 %,
  2–3σ). The trend is the brighter-fatter effect.
- **Agreement with 2008-07, 2009-06 and 2010-03** (`../gain_rdnoise_lbc_200807/`,
  `../gain_rdnoise_lbc_200906/`, `../gain_rdnoise_lbc_201003/`): `gain`
  within −1.1 to +1.1 % and `gain_flux` within −2.1 to +0.3 %. The read noise
  is 0.6–3.6 % higher than in 2008-07 and 2009-06, and 6.7–12 % lower than
  in 2010-03.
- **Disagreement with 2014 and 2025** (`../gain_rdnoise_lbc_201406/`,
  `../gain_rdnoise_lbc_202505/`): on chips 2–4, `gain` is 9–13 % and
  `gain_flux` 7–10 % higher than in those epochs; chip 1 agrees to about
  1 % (`gain`) and 4 % (`gain_flux`). The read noise agrees within 0.4–2.4 %
  on chip 2 and within 2–7 % on chip 4, but is 12–13 % (chip 3) and 18–19 %
  (chip 1) higher in 2011 than in 2014 and 2025. Together with the 2010
  agreement, this puts the LBCB gain change between 2011-04 (MJD 55677) and
  2014-06 (MJD 56830).
- **LBCR is not usable here.**
  - Chip 3 is blank in this epoch: every pixel is 2 ADU in all LBCR flats,
    biases and science frames, so it cannot be measured (the first run
    stopped with "Flats must be brighter than the biases").
  - V-BESSEL sets 9 and 10 (04-29 dawn flats) have `gain_sum` ≈ 0.3 against
    a per-pixel gain of 1.6–1.7, `rho_sum` ≈ 4 and `rho_x`, `rho_y` of
    0.03–0.08. The flat differences are not shot-noise limited; the cause
    was not identified (sky structure, stars or the 14–17 % flux change
    between the frames are candidates).
  - I-BESSEL sets 6–8 are better, but chips 1 and 4 have `rho_x` of −0.02 to
    −0.04 (the sign of the 2009-06 LBCR problem, at a much lower level) and
    set 7 (~30k ADU) has `gain_sum` 12 % below `gain`. The `gain` fit over
    all five sets gives 1.58 ± 0.06, 1.66 ± 0.09 and 1.68 ± 0.09 (chips 1,
    2, 4), 3–4 % below 2010-03 with large errors; read noise 9.2–9.4 e⁻.
    These numbers are not recommended.
- Uncertainties: `gain0_err` (0.1–0.3 %) is statistical only; the choice of
  fit model alone moves the intercept by up to ~1 %. The `gain_flux` errors
  come from only three sets.

## Intermediate data

None beyond the files here (`results_per_set.ecsv`, `gain_fit.ecsv`; the
rejected LBCR measurements are in `results_lbcr_rejected.ecsv`, from an
earlier run of all ten sets).

## How to reproduce

`run.py` follows [`../gain_rdnoise_example/run.py`](../gain_rdnoise_example/run.py):
data path from `$LBCGO_RAW` or `--raw`; outputs and `run_log.json` are
written next to it.

```
export LBCGO_RAW=/Users/howk/Dropbox/Data/LBT/Raw/2011.04/
python calibration/gain_rdnoise_lbc_201104/run.py --box 0
```

The committed run used the LBCB sets only (`inputs.ecsv`). The LBCR sets
were run once together with them (all ten sets) and moved to
`inputs_lbcr_rejected.ecsv`; their per-set results are the file in this
directory. See `whatidid.md` for the working notes.
