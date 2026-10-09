# LBCgo astromatic path: robustness, validation and extended targets

Status: working plan (written 2026-10-09). Companion to the master plan,
`docs/planning/registration_migration_plan.md` ("plan" below), whose §2
decisions apply here. This document plans the work on the **astromatic
path** (SExtractor → SCAMP → SWarp, `LBCgo/lbcregister.py`), which stays
the working path until the native back-end is built.

## 1. Scope and order of work

PI direction (2026-10-09): make the current astromatic approach robust to
the errors we know about and able to combine extended-source fields; the
static distortion model is the last item, because it is the most time
consuming (plan §1.1, D8).

| Stage | Content | Section | Status |
|---|---|---|---|
| Phase 0 | Astromatic path brought to its best state on synthetic tests | §2 | done (record below) |
| A | Robustness to known failure modes; first real-data checks | §3 | open |
| B | Validation harness and baseline numbers on V1–V6 | §4 | open |
| C | Extended-target mode (science sky model feeding SExtractor and SWarp) | §5 | open (SExtractor hook done) |
| later | Static distortion model, native back-end | plan §6 | deferred |

A and B interleave: the harness is how A's fixes are checked on real data.
Order within them: the items that stop a real-data run (A1–A4), then the
harness on V3 (typical field), then the remaining A items as V3–V6 expose
them, then C on V4.

Execution: real-data runs happen on the PI's machine (data volume; cloud
sessions cannot reach Gaia/VizieR, plan §10). Code and synthetic tests can
be written in cloud sessions.

## 2. Phase 0 record (done)

### 2.1 Housekeeping
- [x] `lbcproc.py:20` → `from .lbcregister import *`; confirm `import LBCgo`
      works without `PYTHONPATH` hacks; run tests.
- [x] `pyproject.toml`: `requires-python = ">=3.11"`, `numpy>=2`,
      `astropy>=6.1.4`, `ccdproc>=2.5`; add `scipy`. (New deps added in later
      phases.) Update README/`docs/installation.rst`.
- [x] Accept both astromatic executable names (`sex`/`source-extractor`,
      `swarp`/`SWarp`; `lbcregister.find_astromatic_tool`, used by
      `go_sextractor`, `go_scamp`, `go_scamp_joint`, `go_swarp` and
      `lbcproc.check_external_dependencies`). Ubuntu/Debian apt installs
      `source-extractor` and `SWarp`; SCAMP is `scamp` everywhere.

### 2.2 Weight and mask maps (shared by both back-ends)
- [x] In `go_flatfield`/`go_extractchips`, write per-chip mask + weight
      (implemented in `LBCgo/masks.py`; the bad-pixel list is an optional
      user `badpix_file`, none is shipped yet):
  - saturation: raw ADU ≥ 0.9 × `SATURATE` (must be flagged *before* flat
    division, i.e. carry a mask through `go_overscan` → `go_flatfield`);
  - bad columns: static per-chip bad-pixel list in `conf/` (start from the
    FIXPIX ToDo; derive from flats: pixels deviating > Nσ in the normalized
    master flat);
  - vignetting: normalized flat < 0.5 (tunable) → mask; else weight ∝ flat²
    / sky variance (background-limited inverse variance in flattened units).
- [x] Unit tests with synthetic chips (known bad column, saturated star):
      `tests/test_masks.py`.
- Gain/read noise for the weights come from `LBCgo/detector.py`
  (`conf/lbc_detector.ecsv`, header fallback; weight headers record
  `GAINSRC`). The detector work is recorded separately in
  `docs/detector_gain_rdnoise.md`.

### 2.3 SExtractor
- [x] Pass `-WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE <weight>`; `-FLAG_IMAGE`
      from the mask. (`go_sextractor` picks up the `<base>.weight.fits` /
      `<base>.mask.fits` sidecars automatically; absent sidecars → unweighted
      run. With a flag image, `IMAFLAGS_ISO`/`NIMAFLAGS_ISO` are added to a
      staged copy of the param file: SExtractor fails if they are requested
      without a `FLAG_IMAGE`, so they cannot be in the default param file.)
- [x] Alignment run (default): `-BACK_SIZE 32 -BACK_FILTERSIZE 3`,
      `-DEBLEND_MINCONT 1e-4`, `-DETECT_THRESH 5`, drop the odd
      `-ANALYSIS_THRESH 8`. Make all of these function arguments.
      (`ANALYSIS_THRESH` now follows `DETECT_THRESH` unless given; the conf
      file defaults were changed to match; `go_register(sextractor_args=…)`
      forwards overrides.)
- [~] Extended-target mode: `go_sextractor(subtracted_image=…)` runs on a
      supplied background-subtracted image with `BACK_TYPE MANUAL`,
      `BACK_VALUE 0`. **Pending:** the producer of that image (Stage C, §5).
- Found while implementing:
  - SExtractor does not honour quotes and splits option values at spaces;
    with a package or data path containing a space (e.g. a Dropbox folder)
    `-c` was silently dropped ("not found, using internal defaults") or the
    run failed. `go_sextractor` now stages symlinks in a temp directory in
    that case.
  - The conv/nnw/param files were validated but never passed to `sex`
    (the config's relative `default.conv` only resolved if the cwd held it);
    now passed explicitly.
  - Executable names: `sex`/`source-extractor`, `swarp`/`SWarp` accepted
    via `find_astromatic_tool` (SExtractor, SCAMP, SWarp call sites and
    `check_external_dependencies`).
  - Synthetic check (SExtractor 2.28.2): with the old settings
    (`DEBLEND_MINCONT 0.005`, mesh 64) two stars 40–50 px from a galaxy core
    stay merged into the galaxy segment; with the new defaults both are
    recovered. A zero-weight bad column yields no detections.

### 2.4 SCAMP: one joint run per filter directory
SCAMP's focal-plane modes need one catalog per **exposure** with one
extension per chip. Implemented in `lbcregister.py` (`go_scamp_joint`,
default in `go_register`; `scamp_joint=False` restores per-chip solutions).
**Status: done and validated end to end on a simulated Gaia field with
SCAMP 2.15.0 (conda-forge); the 2.14.1 binary in `/usr/local/bin` is broken
(see below).**
- [x] `merge_ldac(chip_cats, output) -> exposure_cat`: concatenates the
      (`LDAC_IMHEAD`, `LDAC_OBJECTS`) pairs in chip order (astropy `fits`);
      written as `<base>_exp.cat`.
- [x] Run SCAMP **once** on all exposure catalogs of the filter directory:
      iteration 1 `MOSAIC_TYPE LOOSE`, later `FIX_FOCALPLANE`;
      `STABILITY_TYPE INSTRUMENT`; `ASTRINSTRU_KEY FILTER` (drop CFHT's
      `QRUNID`); `-MOSAIC_TYPE` restored (also in the per-chip `go_scamp`);
      the no-op `replace()` deleted (`astrometric_method` is ignored, kept for
      compatibility); `DISTORT_DEGREES 3` kept. SCAMP runs with the catalog
      directory as cwd and relative names (spaces in paths break its option
      parser too); a non-zero exit status raises `RuntimeError` (it was
      silently ignored before).
- [x] `split_head(exposure_head, chip_heads)`: splits the multi-section
      `.head` (sections separated by `END`) into `<base>_<chip>.head`.
- [x] Parse SCAMP output into `astrometry_qa.ecsv` (one row per exposure
      chip: internal/reference rms in arcsec, `FLXSCALE`, `XY_Contrast`,
      reference-match dof, `bad`/`reason`; SCAMP version, reference catalog
      and epoch mode in the table metadata). Per-chip rms comes from the
      `.head` (`ASTIRMS`/`ASTRRMS`, deg), per-exposure numbers from the XML
      `Fields` table — SCAMP's XML has **no per-chip rows**, and
      `AstromSigma_*`/`AstromNDets_*` exist only per field group. Thresholds
      (`max_ref_rms=0.2″`, `min_xy_contrast=2`) are untuned placeholders.
- [x] Proper motions (SCAMP 2.15.0, recorded in the QA metadata):
      `ASTREFEPOCH_TYPE FIELDS_AVERAGE` (what `go_scamp_joint` sets) applies
      Gaia DR3 proper motions from VizieR at the header epoch. Test: 3
      simulated exposures with real Gaia DR3 stars moved by their PM to
      2012.0 (rms PM 12.7, max 50 mas/yr): residual vs truth with
      `FIELDS_AVERAGE` median (−1, +2) mas, same rms as the zero-PM case
      (23/70 mas, limited by the simulation); with `ORIGINAL` the
      high-PM stars sit at a median (−48, −35) mas. **Still to check** (Stage A, §3): that
      real chip headers carry `DATE-OBS` (the simulation had both it and
      `MJD-OBS`; LBC headers have `MJD_OBS`, underscore).
- Validation numbers (2016.0 epoch, no PM, 3 exposures × 2 chips, ~100
  Gaia stars/chip): median residual vs truth (−1, 0) mas; rms 23 (RA) / 69
  (Dec) mas, dominated by the simulation's centroiding of ~90 faint stars.
  QA caveat: SCAMP reports internal rms per instrument and reference rms per
  exposure, so the per-chip rows repeat those values; they are not
  independent per-chip fits.
- **Bugs found while validating:**
  - *Merged catalogs crashed SCAMP* (SIGSEGV/SIGBUS). Cause 1: chip
    catalogs all carry `FITSEXT = 1`, `FITSNEXT = 1`; SCAMP needs them to
    index extensions in the merged file. Cause 2: an astropy round trip of
    `LDAC_IMHEAD` rewrites SExtractor's space-padded 80-character cards as
    NUL-padded, which SCAMP also faults on. `merge_ldac` now works on raw
    FITS bytes and renumbers `FITSFILE/FITSEXT/FITSNEXT` only.
  - Any code that rewrites an LDAC catalog with astropy must preserve the
    space padding.
- **Environment:** `/usr/local/bin/scamp` 2.14.1 (arm64) lacks an
  `LC_RPATH` for `libcurl` and was unstable; it is superseded by the
  conda-forge `astromatic-scamp` 2.15.0 (stable: 12/12 repeated runs and the
  full joint run). See the install notes in the session report.
- SWarp 2.41.5 accepts the SCAMP heads (`CTYPE TAN` + 20 `PV` terms); output
  scale 0.2251″/px for a 0.5 % scale error in the simulation. `astropy.wcs`
  ignores `PV` on `TAN`; astropy consumers must rewrite `CTYPE` to `TPV`
  (as done in the validation script). `go_swarp` path handling was fixed in §2.5.

### 2.5 SWarp
Implemented in `go_swarp` (`go_register(swarp_args=dict(...))` forwards
overrides); the packaged `swarp.lbc.conf` defaults were changed to match.
- [x] Background: `go_swarp(subtract_back=True, back_size=1024)` by default
      (standalone use, since `register/sky.py` (§5) does not exist yet);
      `subtract_back=False` gives `-SUBTRACT_BACK N` for use once the sky is
      removed upstream or for extended targets.
- [x] `-WEIGHT_TYPE MAP_WEIGHT` with the §2.2 `<base>.weight.fits` sidecars
      (SWarp finds them by suffix and fails if only some exist, so weights
      are used only when *every* input has one; otherwise unweighted with a
      message); `-FSCALE_KEYWORD FLXSCALE` when any SCAMP `.head` carries it
      (SWarp merges the `.head` into the header first); `-COMBINE_TYPE
      CLIPPED` default (`CLIP_SIGMA 4`, `CLIP_AMPFRAC 0.3`), `combine_type=
      'MEDIAN'` (also `WEIGHTED`, `AVERAGE`) optional.
- [x] Paths with spaces: SWarp's option parser splits on whitespace, so
      `go_swarp` now builds an argument list and, if any path has a space,
      runs in a temporary directory with symlinked inputs/`.head`/weights/
      config and moves the products back. A non-zero SWarp exit raises
      `RuntimeError`.
- Check (SWarp 2.38.0, synthetic 400² field with a galaxy, a star, a
  half-flux exposure with `FLXSCALE 2`, a cosmic ray, a zero-weight column,
  space in the path): flux scaling reproduces the true galaxy profile
  (`SUBTRACT_BACK N`: r≈120 px level 374 vs 380 true); the cosmic ray is
  rejected by both `CLIPPED` and `MEDIAN`; `BACK_SIZE 128` removes more
  galaxy light than 1024 (centre 1784 vs 1847, and both below the unsubtracted
  2018, as the mesh also removes the galaxy's own pedestal on this small
  frame). Real-data comparison belongs to the harness (§4).
- **Not verified:** that real SCAMP `FLXSCALE` values are sensible for LBC
  data (Stage A, §3).

## 3. Stage A — robustness to known failure modes

Collected from `LBCgo/00ToDo.md`, the README limitations, and the open
checks of Phase 0. Each item gets a synthetic unit test where possible and
is checked on the dataset named.

Stops a real-data run:
- [ ] **A1. Invalid header cards.** `ccdproc.ImageFileCollection` parses
      every card, so one bad card (e.g. `PA_PNT = nan` in the 2010-03 LBCB
      frames) kills the scan. Read only the needed keywords with
      `fits.getheader`, as `calibration/make_inputs.py` (`_scan_headers`)
      does, wherever `lbcproc.py` or `lbcregister.py` (`go_swarp` also
      builds a collection) scan frames. Test with a 2010-03 frame.
- [ ] **A2. Missing or switched-off chips.** Co-pointing or single-chip
      files raise `IndexError` (the code assumes all files have the same
      chips); 2011 data have chip 3 off (V5). Use `LBCCHIP1..4` / the
      extensions actually present per file; carry the chip list per
      exposure through overscan, flat, extraction, `merge_ldac` and the QA
      table; QA flags the absent chip.
- [ ] **A3. Saturated flats.** A filter whose flats are all saturated stops
      the whole run, even if other filters are fine. Skip that filter with
      a warning (and its science frames, or process them unflattened only
      on request).
- [ ] **A4. Flat normalization box.** `make_flatfield` normalizes on a
      hard-coded `data[500:1500, 2000:2500]`; derive the region from the
      chip shape (also needed for partial readouts, A6).

Correctness on real data:
- [ ] **A5. Epoch for proper motions.** SCAMP's `FIELDS_AVERAGE` needs the
      epoch from `DATE-OBS`/`MJD-OBS`; LBC headers have `MJD_OBS`
      (underscore). Check what real chip headers carry after
      `go_extractchips`; if needed write `MJD-OBS` (and `DATE-OBS`) into the
      chip headers. Verify on V3 that high-proper-motion stars match.
- [ ] **A6. Partial readouts / test images.** Detect frames whose size or
      `TRIMSEC` differs from the full 2304 × 4608 layout; skip them with a
      message (or process consistently with A4).
- [ ] **A7. Unmatched flats.** Object frames with no flat in their filter,
      or flats with no object frames, need explicit handling (skip with a
      message, or accept a user-supplied master flat).
- [ ] **A8. V-BESSEL shared by both cameras.** Key flats and output
      directories on (channel, filter) using `detector.lbc_channel`, so one
      run can hold LBCB and LBCR data with the same filter name.
- [ ] **A9. Bad astrometric fits.** `astrometry_qa.ecsv` sets `bad`/
      `reason`, but `go_register` only prints the count. Decide the action
      (default proposal: exclude bad chips from the SWarp input with a
      warning; option to keep them) and tune `max_ref_rms`/`min_xy_contrast`
      on V3/V6. Resolves the ToDo "Auto-identify bad astrometric fits".
- [ ] **A10. Flux scale and unit.** Check that SCAMP `FLXSCALE` values are
      sensible on V3 (compare with the ratio of matched-star fluxes and
      with `EXPTIME`), and make the coadd unit ADU/s (PI preference, plan
      §11 item 4), recorded as `BUNIT`.
- [ ] **A11. Mask thresholds.** Tune `badpix_threshold`/`vignette_threshold`
      on real LBCB/LBCR flats (defaults 0.2 / 0.5 are untested on real
      data). Tune per camera: the LBCR corrector's larger field was designed
      to remove the vignetting seen in LBCB (plan §3.3), so expect the
      threshold to matter mainly for LBCB. Optionally ship a per-chip
      bad-pixel list derived from flats (the FIXPIX ToDo).
- [ ] **A12. Few reference stars (V6).** Check how SCAMP and the QA behave
      with few Gaia stars per chip; the joint solution should carry
      star-poor chips.
- [ ] **A13. `go_scrub` is a no-op** (`rm data/` lacks `-r`, no glob
      expansion, works relative to the cwd). Fix with `shutil`/`pathlib`
      on `image_directory`, so `lbcgo(clean=True)` does what it says.
- [ ] **A14. Chip 3/4 readout layout** (low priority; PI expects it to
      hold). Confirm from raw headers of
      both cameras that chips 3 and 4 have the same readout layout as chips
      1–2 (`NAXIS1/2 = 2304/4608`, `TRIMSEC [51:2098,…]`, `BIASSEC
      [2099:2304,…]`), i.e. that chip 4's 90° rotation is on the sky only
      and `go_overscan`'s `overscan_axis=1` is right for all chips. Record
      chip 4's CD matrix for plan §6.3.2. Also verify the positions of all
      four chips relative to chip 2 from their headers (plan §11 item 2).

## 4. Stage B — validation harness and baseline

The harness (`LBCgo/register/qa.py`) is written against the astromatic
outputs first and is reused unchanged for the native back-end (plan §7).
It needs its own Gaia DR3 access for the absolute residuals:
build `register/refcat.py` (plan §6.3.1: query, cache, epoch propagation,
quality cuts) now rather than in the native phase. The cached catalog can
also be handed to SCAMP as a file (`ASTREF_CATALOG FILE`), which makes runs
reproducible and lets them run without network access.

Datasets (PI to provide paths; see plan §10):
| ID | Content | Purpose |
|----|---------|---------|
| V1 | Many (≥ 20–50) moderately rich, galaxy-free LBCB exposures; selection criteria in plan §6.3.2 "Calibration data" | distortion calibration + astrometry accuracy |
| V2 | Same for LBCR | as V1 |
| V3 | Typical science field (e.g. J1419+4207 OB: 3 dithers × U, g / r, i) | end-to-end regression |
| V4 | NGC 891 (LBCB, `SDT_Uspec` and others) | extended-target mode |
| V5 | Data with chip 3 off (2011) | missing-chip handling |
| V6 | A field with few Gaia stars (high latitude, short exposures) | failure modes |

Metrics written to an ECSV + PDF/PNG figures:
- **Absolute astrometry**: rms and median residual vs Gaia (mas) per
  chip/exposure; residual vector map across the focal plane (systematics).
- **Internal astrometry**: rms of positions of sources matched across
  exposures (including non-Gaia, fainter sources) after registration.
- **Photometry**: per-exposure flux scale; star flux ratio coadd vs
  per-exposure average; agreement native vs SWarp coadd (target < 0.5 %).
- **Image quality**: stellar FWHM in coadd vs median input (target
  ≤ 3 % broadening).
- **Noise**: sky rms in coadd vs prediction from weights; pixel-to-pixel
  correlation (resampling kernel effect).
- **Sky**: mean level and gradient in masked empty regions; for V4, minor-
  and major-axis surface-brightness profiles of NGC 891 compared between
  methods; no negative "moat" around the galaxy.
- **Runtime / peak memory** per stage.

Run order: V3 (typical field; shakes out Stage A), V5 (missing chip), V6
(few stars), V1/V2 (astrometric accuracy over many exposures; the same
data later calibrate the distortion), V4 (NGC 891, before and after
Stage C).

- [ ] `register/refcat.py` with a mocked query in tests; add `astroquery`
      to the dependencies.
- [ ] `register/qa.py`: metrics above from the chip images, the SCAMP
      `.head` files and the SWarp coadd; ECSV + figures.
- [ ] Real-data tests marked `@pytest.mark.realdata`, skipped unless
      `LBCGO_TESTDATA` is set (plan §9).
- [ ] Baseline: the astromatic path runs end-to-end on V1–V6 (or fails with
      a clear message and a QA flag); numbers recorded in
      `docs/planning/baseline_results.md`. Note which plan §7 targets the
      astromatic path already meets (internal ≤ 0.05″, absolute < 0.1″, no
      coherent residual pattern > 10 mas): that sets the priority of the
      native astrometry.

## 5. Stage C — extended-target mode on the astromatic path

Goal: combine fields with a large galaxy (worst case NGC 891, V4) without
over-subtracting its light and without losing the point sources needed
for alignment. Switched on by the user (D4). The sky model
(`LBCgo/register/sky.py`) is shared with the native back-end; here it
feeds SExtractor and SWarp.

Sky model design (moved from plan §6.2):

Default mode (most fields):
- Source mask: `sep.extract` at 1.5σ with a 2–3× dilated segmentation map
  (iterate twice: mask → re-estimate rms → re-detect).
- Model per chip: 2-D polynomial of order ≤ 2 fitted to unmasked pixels
  (sigma-clipped, on a 16×-binned image), or `photutils.Background2D` with a
  very large box (≥ 512 px) as an alternative. Subtract before resampling.

Extended-target mode (`extended_target=True|dict`):
- Galaxy mask: user ellipse, else HyperLEDA (`astroquery.vizier` or
  HyperLEDA query) D25 ellipse scaled × 2 (configurable).
- Model **per exposure across the focal plane**: one 2-D polynomial
  (order ≤ 2) in focal-plane coordinates shared by all chips, plus one
  additive offset per chip, fitted to unmasked pixels of all chips together.
  (Rationale: in the NGC 891 case chips 1, 3, 4 and the outer parts of chip 2
  are galaxy-free; dithers are too small to sample the sky under the galaxy.)
- Diagnostics: fraction of each chip masked; warn if > 60 % of the focal
  plane is masked (sky then unconstrained → recommend offset sky frames).
- Optional: inter-exposure additive offset matching in overlaps (Montage-like
  rectification) to remove residual exposure-to-exposure sky differences.
- Ghosts (LBCB, Giallongo et al. 2008): with the U-LBC interference filter
  (header `FILTER = 'SDT_Uspec'`, so this applies to the NGC 891 U data),
  mask each bright star's ghost (ring 75 px + diffuse 200 px component,
  shifted radially outward, 2.8 % of the star's flux) before fitting the
  sky; the ~0.15 % sky ghost near the field centre is part of the sky model
  or flat, not a source. Not needed for Bessel U, B, V or G, R.

References: Watkins et al. 2024 (masking + parametric modelling vs dithered
stacking); Borlaff et al. 2019 (over-subtraction of extended outskirts);
Trujillo & Fliri 2016; Akhlaghi & Ichikawa 2015.

Wiring into the astromatic path:
- `go_register(extended_target=None | True | dict(ra, dec, a, b, pa))`.
  Per exposure: build the masks, fit the focal-plane sky model, write
  sky-subtracted chips (`<base>_<chip>.skysub.fits`, with the weight/mask
  sidecars and the SCAMP `.head` named to match, since SWarp finds them by
  suffix).
- Two images, two purposes: the **coadd input** is science − sky only; the
  **alignment detection image** additionally has smooth galaxy light
  flattened (median filter, box ≈ 5 × FWHM, as plan §6.1 step 1) so stars on
  the galaxy body are detected. Feed the latter to
  `go_sextractor(subtracted_image=…)` (`BACK_TYPE MANUAL`, `BACK_VALUE 0`;
  hook already in place, §2.3).
- SWarp on the sky-subtracted chips with `subtract_back=False`.
- The ellipse: user-supplied first (works offline); HyperLEDA D25 × 2 as a
  cached network lookup.
- Default (non-extended) mode keeps SWarp's 1024-px background; whether
  `sky.py`'s default per-chip model should replace it is decided from the
  harness numbers (V3).

Tasks:
- [ ] `register/sky.py`: masks (ellipse, iterative source mask), per-chip
      and focal-plane polynomial models, diagnostics (masked fraction,
      warning > 60 %). Add `sep` to the dependencies (plan §8; LGPL note
      in the README).
- [ ] Synthetic tests: four-chip focal plane with a tilted sky, a Sérsic
      galaxy across chip 2 and stars on it; the model recovers the sky
      without a moat; stars on the galaxy body are detected.
- [ ] Wire into `go_register`, `go_sextractor`, `go_swarp` as above.
- [ ] Ghost masks for `SDT_Uspec` (above).
- [ ] V4 acceptance (plan §7 row): no negative moat; sky in a masked
      annulus consistent with 0 within 1σ of the pixel-noise-limited mean;
      ≥ 90 % of Gaia stars projected onto the galaxy body detected and
      matched; major/minor-axis profiles compared with the default mode.
- [ ] (Optional, low-surface-brightness work.) Electronic cross-talk
      ~3 × 10⁻⁵ (Giallongo et al. 2008): a saturated star imprints ~2 ADU
      ghosts in the other chips/channels. Either correct it (needs the
      coefficient matrix) or mask those positions in extended-target mode.

Resolves the ToDo items "Alignment in presence of extended sources" and
"Background estimation in presence of extended sources" for the astromatic
path.

## 6. Later: distortion model and native back-end

Plan §6.3.2 (static per-chip distortion, D2) and the rest of plan §6. To be
carved out into its own planning document when that work starts; the V1/V2
data and the Stage B baseline are its inputs.
