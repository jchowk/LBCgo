# Working notes for Claude Code sessions on LBCgo

(`CLAUDE.md` at the repository root is tracked and loaded automatically;
it points here. Read this file before working on LBCgo.)

LBCgo reduces Large Binocular Camera (LBT) imaging: overscan, bias, flats,
chip extraction, registration and coaddition. PI and owner: J. C. Howk
(astronomer; codes in Python). This file holds what is not obvious from the
code, the plan or the commit history. Read it before working here.

## Where things are

- `docs/planning/registration_migration_plan.md` is the master plan
  (planning only). §1.1 gives the order of work agreed on 2026-10-09; §2
  records PI decisions that are not to be re-litigated without the PI
  (D8 is the order of work). §3.3 holds the verified facts about LBC data,
  §6 the native back-end design (stage D, last), §9.1 the calibration
  provenance rules, §11 the open items for the PI, and §12 the references.
- `docs/planning/astromatic_path_plan.md`: the current work. Record of
  Phase 0 (§2), then stage A (robustness to known failure modes, items
  A1–A14), stage B (validation harness and baseline on V1–V6), stage C
  (extended-target mode on the astromatic path). Its checkboxes are the
  current state.
- `docs/detector_gain_rdnoise.md`: the gain/read-noise work (method,
  results by epoch, open measurements), split out of the plan on
  2026-10-09.
- `LBCgo/detector.py`: the detector table (`conf/lbc_detector.ecsv`), its
  lookups, and the photon-transfer measurement and summary.
  `LBCgo/masks.py`: weight and mask maps. `LBCgo/lbcproc.py`: the
  reduction steps. `LBCgo/lbcregister.py`: the astromatic path
  (SExtractor/SCAMP/SWarp), the working path until stage D.
- `calibration/`: one directory per calibration product (README,
  `inputs.ecsv`, `run.py`, outputs, `run_log.json`). See
  `calibration/README.md` and the README skeleton,
  `NEW_PRODUCT_TEMPLATE.md`. `make_inputs.py` builds `inputs.ecsv`, reading
  only the needed header keywords (raw frames can carry invalid cards).
  Measured products: `gain_rdnoise_lbc_201003/`, `gain_rdnoise_lbc_202505/`.

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

## Gain and read noise: state at 2026-10-10

Details are in `docs/detector_gain_rdnoise.md` and the product READMEs.

- `gain` (zero-level intercept of gain vs flat level; per pixel; weight
  maps) and `gain_flux` (median `gain_sum`; Poisson errors of fluxes).
  The level dependence is the brighter-fatter effect. Realistic `gain`
  uncertainty ~1 %. Run measurements on the whole chip (`run.py --box 0`).
- `conf/lbc_detector.ecsv`: 2006 LBCB rows (Giallongo et al. 2008, open),
  then LBCB-only 2009-06 (MJD 55003; LBCR 2009 unusable), 2010-03 (MJD 55273), 2014-06/07 (MJD 56830) and 2025-05 (MJD 60822)
  rows for both channels. Latest `mjd_start` wins; tests check each block
  against its product's `detector_rows.ecsv`.
- LBCB: 2014 agrees with 2025 to ~1 %; 2010 differs (chips 2–4 gain 8–12 %
  higher, read noise higher), and 2009-06 reproduces the 2010 gains. So the
  change lies between 2010-03 and 2014-06; no known hardware change (PI,
  2026-10-09).
- LBCR: 2014 rows rest on three flat pairs (provisional); 2014 read noise
  is 11–16 % below 2025, cause unchecked.
- For registration, the gain is good enough; further gain work is epoch
  tracking, not a blocker.

## Other open work

Order of work (PI, 2026-10-09; plan §1.1, D8): astromatic path first,
distortion last.
- Phase 0 astromatic improvements: done on synthetic tests.
- Stage A (robustness), B (harness + baseline, needs the PI's V1–V6
  datasets and Gaia/VizieR access, i.e. local runs; plan §10, §11) and C
  (extended-target mode): open; see `docs/planning/astromatic_path_plan.md`.
- Stage D (static distortion model, native back-end, plan §6): deferred;
  to be carved out into its own plan when started.
