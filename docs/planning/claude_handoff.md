# Working notes for Claude Code sessions on LBCgo

(`CLAUDE.md` is git-ignored in this repository. To have Claude Code load
these notes automatically, copy this file to `CLAUDE.md` at the repository
root, or start a session by asking it to read this file.)

LBCgo reduces Large Binocular Camera (LBT) imaging: overscan, bias, flats,
chip extraction, registration and coaddition. PI and owner: J. C. Howk
(astronomer; codes in Python). This file holds what is not obvious from the
code, the plan or the commit history. Read it before working here.

## Where things are

- `docs/planning/registration_migration_plan.md` is the master plan. §2
  records PI decisions that are not to be re-litigated without the PI. The
  checkboxes in §5–§6 are the current state. (The "Status: not yet
  implemented" line at the top is stale; Phase 0 work and the
  detector/gain work are done.) §3.3 holds the verified facts about LBC
  data, §9.1 the calibration provenance rules, §11 the open items for the
  PI, and §12 the references.
- `LBCgo/detector.py`: the detector table (`conf/lbc_detector.ecsv`), its
  lookups, and the photon-transfer measurement and summary.
  `LBCgo/masks.py`: weight and mask maps. `LBCgo/lbcproc.py`: the
  reduction steps. `LBCgo/lbcregister.py`: the astromatic path
  (SExtractor/SCAMP/SWarp), still the default.
- `calibration/`: one directory per calibration product (README,
  `inputs.ecsv`, `run.py`, outputs, `run_log.json`). See
  `calibration/README.md` and the README skeleton,
  `NEW_PRODUCT_TEMPLATE.md` (renamed from `NEW_PRODUCT.md` on 2026-10-07).
  `make_inputs.py` builds `inputs.ecsv` with ccdproc
  `ImageFileCollection`. `gain_rdnoise_lbc_202505/` is the first real
  product.
- These notes describe `main` after branch `202505_calibration` is merged.
  If `calibration/gain_rdnoise_lbc_202505/` is missing on `main`, that
  branch is still unmerged: the 2025 detector rows and the template rename
  are there.

## How the PI works

- **Style.** No flattery, and never "you're absolutely right." Be direct
  and concise. Verify claims, including the PI's, and say plainly when
  something is wrong. Give references (paper, table, section) so results
  can be checked.
- **Branches.** Each task goes on a new `claude/<topic>` branch made from
  the current `origin/main`. Starting a fresh branch after a PR merges is
  the practice, not a PI requirement. Ask before opening a PR; once the PI
  agrees, open it. `git checkout -B` on an existing branch has been
  refused as destructive.
- **Real data.** Real-data runs happen on the PI's machine; raw data are
  not in the repo (e.g. `LBCGO_RAW=/Users/howk/Dropbox/Data/LBT/Raw/2025.05_calib/`).
  The PI pushes the outputs to a branch (calibration runs:
  `202505_calibration`), and Claude reviews them from GitHub. Several
  threads sometimes work in parallel, which can produce merge conflicts.
- **Commits.** Commit messages and docs carry no model identifiers. PR
  bodies keep the repository's own style.

## Testing

- Cloud sessions: `pytest` on PATH may lack numpy. Use a venv with the
  package's dependencies, and give the repo as an **absolute**
  `PYTHONPATH`:
  `PYTHONPATH=$PWD python -m pytest tests -q`. With a relative path,
  `tests/test_calibration_example.py` fails, because it runs `run.py` in a
  subprocess from a temporary directory.
- Practice used for all detector work: check an estimator on a
  simulation before writing its tests; give tests tolerances derived from
  the expected noise and check them over several seeds; then break the
  code deliberately (e.g. median instead of intercept, drop a correction)
  and confirm that a test fails each time.
- The PI's local environment is `/opt/miniconda3/envs/lbcgo`
  (`.claude/settings.json`).

## Verified facts and PI corrections (do not re-derive or contradict)

- LBC: 4 E2V 42-90 CCDs per channel (LBCB blue, LBCR red). **Chip 4 is
  rotated 90° on the sky only**: in the raw readout all chips have the
  same layout with the overscan along x (`overscan_axis=1` is correct for
  every chip). Chip 4's rotation lives in its WCS. (PI correction.)
- `INSTRUME` is spelled inconsistently (`LBC_BLUE` vs `LBC-RED `);
  `detector.lbc_channel` normalizes it.
- The `SDT_Uspec` filter is the U-LBC interference filter (PI).
- `biascheck` is the PROPID of **every** LBC bias, not a special subset
  (PI, 2026-10-06).
- PI: no hardware changes expected over time apart from natural drift;
  the preferred coadd flux unit is ADU/s.
- drizzle 3.0, with its defaults, **preserves surface brightness, not
  counts** (verified by test; an earlier claim otherwise was wrong).
  Sky-flattened frames are in surface-brightness units.
- Environment, D5: Python ≥ 3.11, numpy ≥ 2. Keep SCAMP and SWarp as an
  optional back-end; SourceXtractor++ is not in the core path.

## Gain and read noise: state at 2026-10-08

Details are in plan §5.2 and `calibration/gain_rdnoise_lbc_202505/README.md`.

- Method (`detector.measure_ptc`): photon transfer with unequal flat
  levels (k = μ1/μ2), variances measured in 50-px blocks (whole-region
  variances bias the gain low by 6–12 % on real twilight pairs). Also
  reports nearest-neighbour correlations, `gain_nn`, and `gain_sum` (the
  covariances summed over lags ≤ 3 px and added back). Plane removal per
  block and a −3(1+S)/n bias correction per lag are both required.
- The apparent gain rises linearly with level: +2 % per 10k ADU on LBCB
  and +3–5 % on LBCR. This is the **brighter-fatter effect**, confirmed:
  `gain_sum` is flat on every chip. LBCR non-linearity is ≲ 0.7 % at
  22k ADU.
- **Two gains.** `gain` is the zero-level intercept, per pixel: it sets
  the per-pixel variance and weight maps. `gain_flux` is the median
  `gain_sum` of flat pairs ≤ 60 s apart: use it for Poisson errors of
  fluxes. They differ by −4.4 % (LBCR chip 1, CTI-like positive serial
  correlation) and +3.4 % (LBCR chip 2, electronic anti-correlation;
  this chip also has the high read noise, 12.4 e⁻).
- `gain_sum` is noise-limited: run on the whole chip (`run.py --box 0`).
  It reads 0.5–0.8 % low for pairs taken 140–210 s apart.
- `conf/lbc_detector.ecsv`: the 2006 LBCB rows (Giallongo et al. 2008,
  Table 1; open-ended) plus 8 rows measured on 2025-05-27, valid from MJD
  60822. The latest `mjd_start` wins. They are copied unchanged from
  `calibration/gain_rdnoise_lbc_202505/detector_rows.ecsv` (a test
  checks this). The 2006 rows still apply before MJD 60822 and for an
  unknown date; LBCR before 60822 uses header values.
- Unresolved: the 2025 LBCB gains are 0.82–0.92 × the 2006 values, while
  the read noise in ADU agrees for chip 2, so the electron scales differ.
  Next: check the `GAIN` keywords in 2025 headers and the LBT team's
  values. More epochs are planned by the PI: move `mjd_start` earlier if
  they agree.
- Realistic uncertainty of `gain` is about 1 % (fit model), not
  `gain0_err`.

## Other open work

- Phase 0 astromatic improvements, plan §5.3–5.5: done, except
  extended-target mode for SExtractor (§5.3, partial).
- Not started: the validation harness (§5.6) and the native back-end (§6;
  §6.6 is the switch to the native default).
- The V1–V6 validation datasets will come from the PI (§11). Real-data registration work needs
network access to Gaia/VizieR, which cloud sessions have blocked (§10).
