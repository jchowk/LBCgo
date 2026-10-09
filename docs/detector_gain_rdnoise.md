# LBC gain and read noise: method, results and open work

Split out of `docs/planning/registration_migration_plan.md` (§3.3 and §5.2)
on 2026-10-09 so that the plan holds planning only. This file is the
running record of the detector (gain / read-noise) work: what the code does,
what the measurements show across epochs, and what is still open. Each
measurement itself is documented in its own `calibration/` directory
(convention: `calibration/README.md`); this file summarizes across them.

Why it matters for registration and coaddition: the per-chip gain and read
noise set the inverse-variance weight maps (`LBCgo/masks.py`). In the
sky-limited regime a gain error rescales all exposures of a chip alike and
barely changes the coadd weights; it matters where read noise is not
negligible (LBCB U band) and for absolute flux errors. A ~1 % gain is
therefore more than enough for the registration work; further refinement is
for its own sake (epoch tracking, flux errors).

## 1. Where things are

| What | Where |
|---|---|
| Table read by the pipeline | `LBCgo/conf/lbc_detector.ecsv` |
| Lookup, photon-transfer measurement and summary | `LBCgo/detector.py` (`lookup_detector_params`, `lookup_gain_flux`, `measure_ptc`, `measure_gain_rdnoise_files`, `summarize_gain_rdnoise`) |
| Tests (simulations of brighter-fatter, non-linearity, CTI) | `tests/test_detector.py` |
| Products (one directory per epoch) | `calibration/gain_rdnoise_lbcb_giallongo2008/`, `calibration/gain_rdnoise_lbc_201003/`, `calibration/gain_rdnoise_lbc_202505/` |
| Input-list builder | `calibration/make_inputs.py` |

Table semantics: rows keyed by (channel, chip) with `mjd_start` (inclusive)
/ `mjd_end` (exclusive), NaN = open-ended; the matching row with the latest
`mjd_start` wins. No matching row → `GAIN`/`RDNOISE` header keywords →
nominal defaults. Weight headers record `GAINSRC` (`table`/`header`/
`default`). `gain` is per pixel (variance, weights); `gain_flux` is for
Poisson errors of summed fluxes (nothing in the pipeline uses it yet).

## 2. Current table (2026-10-09)

| Rows | Source | Valid |
|---|---|---|
| LBCB chips 1–4 | Giallongo et al. (2008, Table 1), 2006 commissioning | open (applies before MJD 55273 and when the date is unknown) |
| LBCB + LBCR chips 1–4 | `calibration/gain_rdnoise_lbc_201003/` | MJD 55273 → (superseded at 60822) |
| LBCB + LBCR chips 1–4 | `calibration/gain_rdnoise_lbc_202505/` | MJD ≥ 60822 |

LBCR before MJD 55273 uses header values (`GAIN = 1.75`, `RDNOISE = 12`,
nominal and identical on every chip of both cameras).

## 3. Results across epochs

| | chip 1 | chip 2 | chip 3 | chip 4 |
|---|---|---|---|---|
| LBCB gain 2006 (Giallongo) | 1.96 | 2.09 | 2.06 | 1.98 |
| LBCB gain 2010-03 | 1.80 | 1.97 | 1.93 | 1.82 |
| LBCB gain 2025-05 | 1.80 | 1.72 | 1.74 | 1.66 |
| LBCB read noise 2006 (e⁻) | 11.4 | 11.6 | 11.6 | 11.2 |
| LBCB read noise 2010-03 (e⁻) | 11.7 | 11.0 | 10.9 | 10.4 |
| LBCB read noise 2025-05 (e⁻) | 8.7 | 9.6 | 9.0 | 9.0 |
| LBCR gain 2010-03 | 1.64 | 1.71 | 1.59 | 1.74 |
| LBCR gain 2025-05 | 1.74 | 1.70 | 1.69 | 1.75 |
| LBCR read noise 2010-03 (e⁻) | 9.7 | 9.6 | 8.9 | 9.3 |
| LBCR read noise 2025-05 (e⁻) | 9.8 | 12.4 | 9.5 | 9.3 |

(Gains in e⁻/ADU, zero-level intercepts; see the product READMEs for errors
and `gain_flux`.)

- **LBCB changed between 2010 and 2025.** Gains 2010/2025 = 1.00, 1.15,
  1.11, 1.10; read noise 15–35 % higher in 2010. The 2010 chip-to-chip
  pattern (chip 2 highest, chip 1 lowest) matches 2006, not 2025. The 2010
  values are 0.92–0.94 × the 2006 ones, partly explainable by the level at
  which Giallongo et al. measured (§4, brighter-fatter). Read noise in ADU
  in 2025 agrees with 2006 for chip 2 while the gains differ by 16–18 %:
  the electron scales differ, not the ADC conversion alone.
- **LBCR is roughly stable:** gains agree within ~6 %, read noise within
  ~6 % except chip 2 (9.6 → 12.4 e⁻).
- No known hardware or controller change (PI, 2026-10-09). The PI is
  measuring intermediate epochs to see whether the LBCB change is a step or
  a drift (§5).

## 4. Method and what was learned

- Photon transfer with unequal flat levels (k = μ1/μ2); variances measured
  in 50-px blocks. Whole-region variances bias the gain low by 6–12 % on
  real twilight pairs (illumination differences between the two flats: −6 %
  for a 0.5 % peak-to-peak mismatch over 1000 px in simulation, < 0.5 % with
  blocks). Pick consecutive flats of one sequence at the same rotator angle;
  avoid z/Y-band flats (fringing); use 10,000–30,000 ADU flats.
- **Brighter-fatter.** The apparent gain rises linearly with level: +1.8–2.7
  % per 10k ADU on LBCB, +2–5 % on LBCR (both epochs), several times what
  the < 1 % linearity residuals of Giallongo et al. (2008) allow. Adopted
  `gain` = zero-level intercept g0 of gain vs level (median if fewer than 3
  sets or no spread in level); `rdnoise` = g0 × median read noise in ADU.
  References: Antilogus et al. 2014, JINST 9, C03048; Astier et al. 2019,
  A&A 629, A36.
- **Diagnostics.** `measure_ptc` reports `rho_x`, `rho_y` (lag-1 correlation
  of the flat difference), `gain_nn` (lag-1 covariances added back) and,
  with `max_lag=3`, `rho_sum` / `gain_sum` (covariances summed over
  |dx|, |dy| ≤ 3, plane removed per block, −3(1 + S)/n bias corrected per
  lag). Charge conservation makes `gain_sum` free of brighter-fatter within
  the lag range; non-linearity creates no covariances and survives.
  Simulations: half the charge shared at 2 px → gain slope 5.6 %/10k ADU,
  `gain_nn` 2.4, `gain_sum` 0.1 ± 0.5; sublinear response → gain slope 4.9
  stays 4.9 ± 0.6 in `gain_sum`. Noise ~1 % per set on a 1000² box, ~0.3 %
  on the whole chip: run with `run.py --box 0`.
- **2025 whole-chip result:** `gain_sum` slopes −0.78 to +0.05 %/10k ADU on
  every chip, so the LBCR trend is brighter-fatter with covariances beyond
  lag 1 (`rho_sum` grows 0.035–0.054 per 10k ADU on LBCR, 0.020–0.026 on
  LBCB), not non-linearity: LBCR non-linearity ≲ 0.3 % at 10k ADU and ≲ 0.7
  % at 22k ADU. g0 changes by ≤ 0.13 % when widely spaced pairs are
  dropped. (The 2010 LBCR `gain_sum` slopes have 3–4 % errors and cannot
  separate the two.)
- **Two gains.** Correlations present at zero signal make the per-pixel
  and flux gains differ: 2025 median `gain_sum`/g0 = −4.6 % on LBCR chip 1
  (positive serial correlation at low level, falling with level: CTI-like),
  +2.8 % on LBCR chip 2 (serial anti-correlation ρ_x ≈ −0.020, electronic;
  the chip with the high read noise and `bias_rho_x` ≈ −0.05), −1.2 to
  +1.2 % elsewhere; both reproduced by 1/(1 + S) with S the summed
  correlation at zero signal. For a linear readout kernel with weights
  summing to H, an aperture sum's mean scales as H and its variance as H²,
  so `gain_sum` (once `max_lag` covers the kernel) is the electrons per ADU
  of a flux. Hence `gain` for per-pixel variance and weights, `gain_flux`
  (median `gain_sum`) for flux errors. `gain_flux` is an optional column
  (`detector.OPTIONAL_TABLE_COLUMNS`; NaN = unknown; `lookup_gain_flux`
  falls back to `gain` and says so).
- **Correlated read noise:** `bias_rho_sum` 0.03–0.24 (LBCR chip 2
  negative), so read noise in an aperture is up to ~11 % above the
  independent-pixel value; per-pixel read noise is unaffected.
- **Flat spacing.** `gain_sum` depends on the time between the two flats
  (2025: pairs 41–60 s apart +0.3 to +0.6 %, 143–212 s apart −0.5 to
  −0.8 %); g0 is insensitive. `summarize_gain_rdnoise(flux_max_dt=60)` uses
  only pairs ≤ 60 s apart for `gain_flux`; `run.py` records `flat_dt` and
  warns when the flats come from different OBs.
- **Uncertainty.** The intercept moves by up to 1.1 % between a linear and
  a quadratic fit: realistic `gain` uncertainty ~1 %, not `gain0_err`
  (0.1–0.4 %).
- **Other 2025 findings:** a level-independent lag-1 anti-correlation along
  x (ρ_x down to −0.019, LBCR chip 2) pushes `gain_nn` 1–8 % above the
  intercept; LBCB chips 2–3 read noise rises ~5 % through the 01:28–01:41 UT
  bias sequence (cause unknown; `biascheck` is the PROPID of every LBC bias,
  PI 2026-10-06).
- **2010 findings:** the first LBCR biases of a block show a decaying level
  (readout transient; drop them); one LBCR bias block had an image-area
  pattern drift (higher bias-difference noise with normal overscan rms) and
  was replaced; four of six LBCR sets are I band (fringing may bias the
  variances); LBCR flats span only 13–27k ADU, so the intercept has 2–3 %
  errors.
- Raw frames may carry invalid header cards (`PA_PNT = nan`, 2010-03);
  `make_inputs.py` reads only the needed keywords.

## 5. Open work

- [ ] **Bracket the LBCB change** (PI, in progress): measure epochs between
      2010 and 2025. If a step is found, give the 2010 rows a finite
      `mjd_end` and add rows for the new state; if a drift, add rows per
      epoch. Until then, the 2010 rows apply up to MJD 60822, which is
      well supported for LBCR and doubtful for LBCB.
- [ ] Compare with the `GAIN`/`RDNOISE` keywords of the 2025 headers and any
      values the LBT/LBC team publishes (the 2010 commissioning table seen
      only in search snippets: LBCB gains as Giallongo, LBCR 2.08–2.14
      e⁻/ADU, LBCB read noise 4.8–5.2 ADU; unverified).
- [ ] LBCR: add V/r-band flat pairs to check the I-band fringing concern of
      2010, and lower-level pairs (1–4k ADU) for a shorter extrapolation.
- [ ] More data generally: consecutive pairs ≤ 60 s apart at the same
      rotator angle (for `gain_flux`); biases from later in the night
      (LBCB chips 2–3 read-noise drift).
- [ ] Use `gain_flux` in the pipeline once catalogs need flux errors
      (registration plan, catalog stage).

Procedure for a new epoch: `calibration/README.md` and `CLAUDE.md`
(make_inputs → edit `inputs.ecsv` → `run.py --box 0` → review
`detector_rows.ecsv` → merge into `conf/lbc_detector.ecsv` as a separate
step; a test checks the merged rows).

## 6. References

- Giallongo, E. et al. 2008, A&A, 482, 349 — LBCB per-chip gain and read
  noise (Table 1), linearity, full well.
- Speziali, R. et al. 2008, Proc. SPIE 7014, 70144T — LBCR read noise
  "< 10 e⁻ @500 Kpix/s/ch".
- Antilogus, P. et al. 2014, JINST, 9, C03048; Astier, P. et al. 2019, A&A,
  629, A36 — brighter-fatter effect and covariances in flat pairs.
