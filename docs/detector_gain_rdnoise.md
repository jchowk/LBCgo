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
| Products (one directory per epoch) | `calibration/gain_rdnoise_lbcb_giallongo2008/`, `calibration/gain_rdnoise_lbc_200807/`, `calibration/gain_rdnoise_lbc_200906/`, `calibration/gain_rdnoise_lbc_201003/`, `calibration/gain_rdnoise_lbc_201406/`, `calibration/gain_rdnoise_lbc_202505/` |
| Input-list builder | `calibration/make_inputs.py` |

Table semantics: rows keyed by (channel, chip) with `mjd_start` (inclusive)
/ `mjd_end` (exclusive), NaN = open-ended; the matching row with the latest
`mjd_start` wins. No matching row → `GAIN`/`RDNOISE` header keywords →
nominal defaults. Weight headers record `GAINSRC` (`table`/`header`/
`default`). `gain` is per pixel (variance, weights); `gain_flux` is for
Poisson errors of summed fluxes (nothing in the pipeline uses it yet).

## 2. Current table (2026-10-10)

| Rows | Source | Valid |
|---|---|---|
| LBCB chips 1–4 | Giallongo et al. (2008, Table 1), 2006 commissioning | open (applies before MJD 54652 and when the date is unknown) |
| LBCB chips 1–4 (LBCR unusable) | `calibration/gain_rdnoise_lbc_200807/` (four non-independent pairs, 25–43k ADU) | MJD 54652 → (superseded at 55003) |
| LBCB chips 1–4 (LBCR unusable) | `calibration/gain_rdnoise_lbc_200906/` | MJD 55003 → (superseded at 55273) |
| LBCB + LBCR chips 1–4 | `calibration/gain_rdnoise_lbc_201003/` | MJD 55273 → (superseded at 56830) |
| LBCB + LBCR chips 1–4 | `calibration/gain_rdnoise_lbc_201406/` (LBCR provisional: three flat pairs) | MJD 56830 → (superseded at 60822) |
| LBCB + LBCR chips 1–4 | `calibration/gain_rdnoise_lbc_202505/` | MJD ≥ 60822 |

LBCR before MJD 55273 (no LBCR rows exist for 2008 or 2009) uses header values (`GAIN = 1.75`, `RDNOISE = 12`,
nominal and identical on every chip of both cameras).

## 3. Results across epochs

| | chip 1 | chip 2 | chip 3 | chip 4 |
|---|---|---|---|---|
| LBCB gain 2006 (Giallongo) | 1.96 | 2.09 | 2.06 | 1.98 |
| LBCB gain 2008-07 | 1.79 | 1.96 | 1.92 | 1.81 |
| LBCB gain 2009-06 | 1.79 | 1.94 | 1.93 | 1.80 |
| LBCB gain 2010-03 | 1.80 | 1.97 | 1.93 | 1.82 |
| LBCB gain 2014-06 | 1.81 | 1.74 | 1.75 | 1.67 |
| LBCB gain 2025-05 | 1.80 | 1.72 | 1.74 | 1.66 |
| LBCB read noise 2006 (e⁻) | 11.4 | 11.6 | 11.6 | 11.2 |
| LBCB read noise 2008-07 (e⁻) | 10.0 | 9.6 | 10.0 | 9.4 |
| LBCB read noise 2009-06 (e⁻) | 10.2 | 9.8 | 9.8 | 9.5 |
| LBCB read noise 2010-03 (e⁻) | 11.7 | 11.0 | 10.9 | 10.4 |
| LBCB read noise 2014-06 (e⁻) | 8.8 | 9.8 | 9.1 | 9.4 |
| LBCB read noise 2025-05 (e⁻) | 8.7 | 9.6 | 9.0 | 9.0 |
| LBCR gain 2010-03 | 1.64 | 1.71 | 1.59 | 1.74 |
| LBCR gain 2014-07 (provisional) | 1.74 | 1.76 | 1.73 | 1.78 |
| LBCR gain 2025-05 | 1.74 | 1.70 | 1.69 | 1.75 |
| LBCR read noise 2010-03 (e⁻) | 9.7 | 9.6 | 8.9 | 9.3 |
| LBCR read noise 2014-07 (e⁻) | 8.3 | 8.1 | 8.4 | 7.9 |
| LBCR read noise 2025-05 (e⁻) | 9.8 | 12.4 | 9.5 | 9.3 |

(Gains in e⁻/ADU, zero-level intercepts; see the product READMEs for errors
and `gain_flux`.)

- **LBCB: 2014 agrees with 2025, 2010 does not.** 2014 and 2025 agree to
  ~1 % (`gain` +0.5 to +1.3 %, `gain_flux` ≤ 0.3 %, read noise in ADU the
  same), from different flats and biases 11 years apart. The 2010 gains of
  chips 2–4 are 8–12 % higher than 2014 (chip 1 the same); 2010/2025 =
  1.00, 1.15, 1.11, 1.10, read noise 15–35 % higher in 2010. So whatever
  differs lies between 2010-03 and 2014-06: a change in the detectors or
  electronics, or a biased 2010 measurement. The 2010 chip-to-chip
  pattern (chip 2 highest, chip 1 lowest) matches 2006, not 2025. The 2010
  values are 0.92–0.94 × the 2006 ones, partly explainable by the level at
  which Giallongo et al. measured (§4, brighter-fatter). Read noise in ADU
  in 2025 agrees with 2006 for chip 2 while the gains differ by 16–18 %:
  the electron scales differ, not the ADC conversion alone.
- **Read noise in ADU is the stable quantity.** It comes from the biases
  alone, independent of the flats. LBCB chips 2 and 4 read 5.6–5.7 ADU in
  all four epochs (2006, 2010, 2014, 2025; chip 4 5.4 in 2025), while
  their flat-based gains fall by 12–17 % between 2010 and 2014. Either the
  dominant read noise arises after the stage whose gain changed, or the
  2006 and 2010 flat-based gains are biased high in the same way. Chip 1 is
  the exception: same gain in every epoch, but 6.5 ADU in 2010 against
  4.8 in 2014 and 2025. Worth checking against the raw 2010 flats before
  attributing the difference to the instrument.
- **LBCR gains are roughly stable:** 2014 `gain_flux` agrees with 2025
  within 1 % on every chip, `gain` within 0.2–3.2 %; 2010 is 2–9 % lower.
  **LBCR read noise varies:** 2014 reads 4.4–4.9 ADU (7.9–8.4 e⁻), 11–16 %
  below 2025 (chip 2: 35 %, since 2025 has 12.4 e⁻ there) and 5–16 % below
  2010. The same method reproduces the LBCB read noise between 2014 and
  2025, so this is probably real; not checked.
- No known hardware or controller change (PI, 2026-10-09).

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
- **2014 findings** (selection lessons for future epochs):
  - All LBCB flats of 2014-06-25 (3.28 s and 1.28 s pairs) had strongly
    spatially correlated differences (`rho_sum` 0.4–2.3, against 0.03–0.08
    for good pairs), giving `gain_sum` ≈ 0.5–0.7; they were excluded.
    Cause unknown. `rho_sum` is a useful screen for bad pairs.
  - LBCR pairs more than 50 s apart (sky level changing 18–24 % between
    the two) gave low `gain_sum`; pairs ≤ 38 s apart were consistent. This
    supports, and may argue for tightening, the 60 s `flux_max_dt` cut.
  - Nights without biases: biases can be fetched from the LBT archive
    (`archive.lbto.org`); the 2014 product used a 2014-06-24 bias sequence
    for flats up to 7 days away, assuming bias level and read noise are
    constant.
  - Three sets is the minimum for `summarize_gain_rdnoise` (`min_sets`);
    with three the line has one degree of freedom and the quoted `gain0`
    errors are not meaningful.
- **2008 findings** (tars of the archive calibration data; only 37 of the 96
  logged frames, later 8 more loose LBCR files):
  - LBCB (five SDT_Uspec flats at 25–43k ADU, so four overlapping pairs)
    agrees with 2009-06 (`gain` +1.0 to −0.4 %, `gain_flux` 0 to +1.4 %) and
    2010-03 (`gain` within 0.8 %): the LBCB gains are stable from 2008-07 to
    2010-03, and the change to the 2014/2025 values lies in 2010-03 → 2014-06.
  - LBCR 2008-07 shows the 2009 behaviour: `rho_x` −0.19 to −0.32 on chips
    1, 3, 4 (−0.04 on chip 2) in i-SLOAN (10.2 s) and R-BESSEL (3.2 s) flats
    alike, independent of level; bias differences are clean
    (`bias_rho_x` ≥ −0.04 in 2008 and 2009), so the effect is
    signal-dependent. Per-pixel gains are 0.63–1.13 (chips 1, 3, 4). The
    R-BESSEL `gain_sum` (1.68, 1.61, 1.64, 1.60; constant over 29–40k ADU)
    is plausible as a flux gain, the i-SLOAN `gain_sum` (1.25–1.4) is not.
    LBCR before 2010-03 is therefore unusable for the per-pixel gain; two
    epochs, two filters and two exposure times show the same signature, so
    more pre-2010 LBCR flats are unlikely to help. Whether the excess
    noise appears at sky levels in science frames is untested.
  - Frames needed for the LBCR measurement above 0.7 × SATURATE (r-SLOAN
    42–65k ADU) cannot be used by `measure_gain_rdnoise_files`.
  - Extract archive tars outside Dropbox (they are 1.5 GB); biases were 5 h
    before the flats.
- **2009 findings:**
  - LBCB (10 sets, 2009-06-21 flats, 2009-06-22 biases) agrees with 2010-03
    (`gain` −0.1 to −1.8 %, `gain_flux` −2.1 to +0.2 %), not with 2014/2025
    (chips 2–4 8–13 % higher). An independent epoch, with different flats
    and biases, thus reproduces the 2010 gains: a bias of the 2010
    measurement is unlikely, and the change lies between 2010-03 and
    2014-06. Read noise is 8.5–13 % below 2010 and chip 1 is 17 % above 2014/2025.
  - LBCR 2009-06 is unusable: all eight sets (z- and i-SLOAN) have
    `rho_x` −0.27 to −0.42 in the flat difference and an rms 2.7 × the
    shot-noise expectation, with power rising towards the x Nyquist
    frequency; the headers match those of LBCB and later LBCR data. Cause
    unknown. The rejected sets are in `inputs_lbcr_rejected.ecsv` and
    `results_lbcr_rejected.ecsv` of the product directory.
  - Biases come from 22 h after the flats (assumed stable).
- Raw frames may carry invalid header cards (`PA_PNT = nan`, 2010-03);
  `make_inputs.py` reads only the needed keywords.

## 5. Open work

- [x] 2014-06/07 epoch measured and merged (`calibration/gain_rdnoise_lbc_201406/`,
      rows from MJD 56830).
- [x] 2008-07 epoch measured and merged for LBCB (`calibration/gain_rdnoise_lbc_200807/`,
      rows from MJD 54652); LBCR unusable (as in 2009-06).
- [x] 2009-06 epoch measured and merged for LBCB (`calibration/gain_rdnoise_lbc_200906/`,
      rows from MJD 55003); LBCR unusable.
- [ ] **Bracket the LBCB change** (PI, in progress), within 2010-03 →
      2014-06 (MJD 55273–56830). Measure epochs in that window. A step →
      the 2010 rows end at the step (they already yield to the 2014 rows
      from 56830); no step anywhere in the window, and a re-check of the
      2010 flats finds a bias → replace the 2010 LBCB rows.
- [ ] Re-check the 2010 LBCB measurement (single filter, `SDT_Uspec`;
      chip 1 read noise in ADU 35 % above other epochs): look for flat
      structure or bias problems as in the 2014-06-25 flats. The 2009-06
      epochs of 2008-07 and 2009-06 reproduce the 2010 gains, which argues
      against a bias in them.
- [ ] LBCR read noise: 2014 is 11–16 % below 2025 and 5–16 % below 2010.
      Check the bias frames of each epoch (bias-difference rms vs overscan
      rms, as done for 2010 set 10) before treating it as an instrument
      change.
- [ ] LBCR 2014 rows are provisional (three flat pairs): replace with a
      deeper LBCR set from the same period if one exists.
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
