# LBCB per-chip gain and read noise, 2013-06 (LBCB only)

## Summary

- **Product file(s):** `detector_rows.ecsv`, merged into
  `LBCgo/conf/lbc_detector.ecsv` (4 LBCB rows, `mjd_start` = 56449)
- **SHA-256:** `12d653c45cd72d7376be9e2af6f417801c9a75bde71801f9bdd46f185e6595c2  detector_rows.ecsv`
  (as written by `run.py`; `PARAMS['mjd_start']` = 56449.0 is set in `run.py`,
  no hand edit)
- **Made by / date:** J.C. Howk / 2026-10-11 (`run_log.json`)
- **LBCgo commit:** `4213c7584d3a0aba1e6da29dc24a3ed71e2b0606` (the code
  `run.py` imported, as in `run_log.json`; LBCgo 0.2.0.dev0)
- **Supersedes:** for LBCB, the 2011-04 rows (`../gain_rdnoise_lbc_201104/`)
  from MJD 56449 on; those remain in force from MJD 55677. These rows are in
  turn superseded by `../gain_rdnoise_lbc_201406/` from MJD 56830.
- **Valid for:** LBCB, all filters, MJD ≥ 56449 (2013-06-06 UT, the date of
  the flats used); the 2014-06 rows take over at MJD 56830. **No LBCR
  rows**: the LBCR data of this epoch were not usable (see Results); LBCR
  stays on the 2010-03 rows.

## Purpose

Per-chip gain and read noise for the weight maps (`gain`, `rdnoise`), and
the flux gain (`gain_flux`) for Poisson errors of source fluxes, measured
with the photon-transfer method on raw twilight flats and biases of
2013-06-06 UT. It bridges the LBCB gain change on chips 2–4 (about 10 %)
seen between the 2011-04 and 2014-06 products, and dates it to before
MJD 56449.

## Inputs

`inputs_full.ecsv` lists every frame in the data directory (110 frames: 48
biases, 50 flats, 12 science frames):

```
python calibration/make_inputs.py $LBCGO_RAW -o inputs_full.ecsv --overwrite --levels
```

`inputs.ecsv` is that list edited down to 7 LBCB sets of two flats and two
biases (all frames 2013-06-06 UT):

- **Flats (dusk, 02:47–03:05 UT):** sets 1–3 are U_spec (`SDT_Uspec`) sky
  flats (1.28 s at ~26–27k ADU; 1.78 s at ~22k and ~16.5k ADU); sets 4–7 are
  r-SLOAN sky flats (2.28 s at ~22–25k and ~16–17k ADU; 8.27 s at ~27–30k
  and ~19–21k ADU, as the mean level of the chips). Each set is two
  consecutive flats of one OB (equal exposure time) 39–48 s apart; no frame
  is used in two sets. The twilight ramps are falling, so the level changes
  by 10–15 % between the two frames of a set. Frames above 34k ADU were not
  used.
- **Biases:** consecutive pairs of the `25Bias_Bino` sequence (`propid`
  `biascheck`) taken at 11:50–12:04 UT, skipping the first two frames, a
  different pair for each set. They were taken about 9 h after the dusk
  flats; the bias level and read noise are assumed constant over that
  interval.

Rejected, with the per-set measurements kept for the record:

- `inputs_lbcb_rejected.ecsv`, `results_lbcb_rejected.ecsv`: LBCB sets 8–10,
  g-SLOAN **dawn** flats (11:19–11:30 UT).
- `inputs_lbcr_rejected.ecsv`, `results_lbcr_rejected.ecsv`: all 8 LBCR sets
  (r-SLOAN dusk sets 11–14, i-SLOAN dawn sets 15–18).

## Method

`run.py` (the 2011-04 driver with `mjd_start` and the chip list changed),
run on the whole trimmed chip (`--box 0`):

- `LBCgo.detector.measure_gain_rdnoise_files` per set: overscan-subtracted
  photon-transfer gain and read noise with variances in 50-px blocks, plus
  correlation diagnostics and `gain_sum` (covariances summed over lags
  ≤ 3 px).
- `LBCgo.detector.summarize_gain_rdnoise` per chip: `gain` = zero-level
  intercept of a straight-line fit of gain against level; `rdnoise` =
  `gain` × median read noise in ADU; `gain_flux` = median `gain_sum` of the
  sets whose two flats are ≤ 60 s apart (all seven sets).

## Results

| chip | `gain` | `rdnoise` (e⁻) | `gain_flux` |
|---|---|---|---|
| 1 | 1.820 ± 0.006 | 9.91 | 1.794 ± 0.004 |
| 2 | 1.736 ± 0.008 | 11.12 | 1.747 ± 0.005 |
| 3 | 1.758 ± 0.009 | 10.45 | 1.756 ± 0.007 |
| 4 | 1.662 ± 0.003 | 10.12 | 1.657 ± 0.005 |

- The apparent gain rises with level, by +1.8 to +2.2 % per 10,000 ADU
  (rms about the fit 0.2–0.4 %; levels 16–30k ADU). With the covariances
  summed to 3 px the slope is consistent with zero (−0.6 to +0.8 %, errors
  0.4–0.9 % per 10k). The trend is the brighter-fatter effect.
- **Agreement with 2014 and 2025** (`../gain_rdnoise_lbc_201406/`,
  `../gain_rdnoise_lbc_202505/`): `gain` within −0.4 to +1.1 % on all four
  chips and `gain_flux` within −0.2 to +0.3 %. The read noise is 7–16 %
  higher here than in those epochs.
- **Disagreement with 2008-07 through 2011-04** (`../gain_rdnoise_lbc_200807/`,
  `../gain_rdnoise_lbc_200906/`, `../gain_rdnoise_lbc_201003/`,
  `../gain_rdnoise_lbc_201104/`): on chips 2–4, `gain` is 8–12 % and
  `gain_flux` 7–10 % lower here; chip 1 agrees to 1.4–1.8 % in `gain`
  (2–4 % in `gain_flux`). Together with the 2011-04 and 2014-06 products,
  this puts the LBCB gain change between 2011-04 (MJD 55677) and 2013-06
  (MJD 56449).
- Read noise (chip 2: 11.1 e⁻; 9.8 in 2011-04 and 2014-06) varies between
  epochs by up to ~15 % in a way that does not follow the gain change.
- **Dawn flats are unusable.** The g-SLOAN sets 8–10 (11:19–11:30 UT, rising
  sky level) have `gain_sum` of 0.86–1.24 e⁻/ADU against a per-pixel gain of
  1.7–1.9, `rho_sum` of 0.5–1.1, and `rho_x`, `rho_y` of 0.01–0.02
  (several times the dusk values). The same failure affected the dawn flats of
  2011-04. Including them in the fit leaves the central gains almost
  unchanged but inflates the errors 3–5 times and makes the `gain_sum` slope
  meaningless. The cause was not identified.
- **LBCR is not usable.**
  - The dusk r-SLOAN sets (11–14) have `rho_sum` of 0.09–0.18 (about three
    times the LBCB value) and `gain_sum` 8–13 % below `gain`. Their `gain`
    (1.86–2.01) is above both the 2010-03 (1.59–1.74) and 2014-06
    (1.73–1.78) values, and chips 1 and 2 have a read noise of 11.2–12.6 e⁻
    against 9.4–10.0 e⁻ on chips 3 and 4. Their levels span only 20–32k ADU.
  - The dawn i-SLOAN sets (15–18) have `rho_sum` of 1.2–2.1, `rho_x` of
    0.02–0.06 and `gain_sum` of 0.55–0.85 e⁻/ADU.
- Uncertainties: `gain0_err` (0.2–0.5 %) is statistical only; the choice of
  fit model alone moves the intercept by up to ~1 %. The fit spans only
  16–30k ADU.

## Intermediate data

None beyond the files here (`results_per_set.ecsv`, `gain_fit.ecsv`; the
rejected measurements are in `results_lbcb_rejected.ecsv` and
`results_lbcr_rejected.ecsv`, from an earlier run of all eighteen sets).

## How to reproduce

`run.py` follows [`../gain_rdnoise_example/run.py`](../gain_rdnoise_example/run.py):
data path from `$LBCGO_RAW` or `--raw`; outputs and `run_log.json` are
written next to it.

```
export LBCGO_RAW=/Users/howk/Dropbox/Data/LBT/Raw/2013.06/
python calibration/gain_rdnoise_lbc_201306/run.py --box 0
```

The committed run used the seven dusk LBCB sets only (`inputs.ecsv`). The
rejected sets were run once together with them and moved to the
`inputs_*_rejected.ecsv` files.
