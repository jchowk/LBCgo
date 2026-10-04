# LBCgo registration & coaddition: migration plan

Status: **agreed plan, not yet implemented** (written 2026-10-04).
Audience: a later implementation session (human or Claude). Read this whole
document before writing code; section 2 records decisions that are not to be
re-litigated without the PI (J. C. Howk).

---

## 1. Goal

Replace the hard dependency on the astromatic tools (SExtractor, SCAMP, SWarp)
in `LBCgo/lbcregister.py` with an in-process Python implementation, while
keeping the astromatic tools as an optional, *improved* back-end used for
validation. The new path must:

1. Register LBC chip images to Gaia DR3 with **internal (exposure-to-exposure)
   registration ≤ 0.05″ rms** and **absolute residuals vs Gaia < 0.1″ rms**.
2. Produce a coadd with correct flux scaling, weight maps and outlier
   rejection.
3. Handle fields containing a large galaxy (worst case: NGC 891, D25 ≈ 13.5′,
   lying diagonally across chip 2) without over-subtracting galaxy light and
   without losing the point sources needed for alignment.
4. (Secondary) Produce final source catalogs from the coadd.

## 2. Decisions already made

| # | Decision |
|---|----------|
| D1 | Python implementation is the **default** back-end. SExtractor/SCAMP/SWarp remain as an **optional** back-end (`backend='astromatic'`), used for validation and as a fallback. |
| D2 | Astrometry uses a **static per-chip distortion model** (per channel LBCB/LBCR, per filter if needed, per epoch if needed) fitted once from the PI's calibration datasets and shipped in the package, plus a **per-exposure linear correction fitted jointly across all chips** against Gaia DR3. |
| D3 | Default combine = **per-exposure flux-scaled, weighted, sigma-clipped mean**; median available as an option. |
| D4 | Large-galaxy handling is an **optional "extended-target" mode** switched on by the user (optionally with a user-supplied ellipse, else HyperLEDA D25 × ~2). |
| D5 | Minimum environment: **Python ≥ 3.11, numpy ≥ 2**. |
| D6 | Phase 0 first: bring the existing astromatic path to its best state so it is a fair benchmark for the Python path. |
| D7 | SourceXtractor++ is **not** used in the core path (conda-only, its docs recommend a separate environment; it replaces only SExtractor). It may later be an optional back-end for science catalogs. |

## 3. Current state (verified 2026-10-04, commit `4b9efbc`)

### 3.1 Pipeline flow
`lbcgo()` → `go_overscan` → `go_flatfield` → `make_targetdirectories` →
`go_extractchips` (writes `<base>_<chip>.fits`, single-extension, one per chip;
moves the MEF `*_flat.fits` into `data/`) → `go_register(fltr_dirs, …)`.

`go_register` (`LBCgo/lbcregister.py:421`) loops over chip files and, **for each
chip file individually**, runs `go_sextractor` then `go_scamp`; then runs
`go_swarp` on all chip files of the filter directory.

### 3.2 Defects / limitations in the astromatic usage
| Location | Issue |
|----------|-------|
| `lbcregister.py:537` | SCAMP is invoked **once per chip catalog**, so every chip is solved independently with a 3rd-order polynomial (`scamp.lbc.conf:66`, `DISTORT_DEGREES 3`). No focal-plane or cross-exposure constraint. Under-constrained on star-poor chips → likely source of the "bad astrometric fits" ToDo item. |
| `lbcregister.py:248` | `-MOSAIC_TYPE` flag commented out; the per-iteration `mosaic_type` values are dead code. `scamp.lbc.conf:48` → `UNCHANGED`. |
| `lbcregister.py:251` | `cmd_flags.replace('INSTRUMENT','EXPOSURE')` discards its result (and the string never contains `INSTRUMENT`): `astrometric_method` is a no-op. |
| `lbcregister.py:119` | Command-line `-DETECT_THRESH 5.0 -ANALYSIS_THRESH 8.0` overrides the config's 1.5/1.5. |
| `sextractor.lbc.conf:23,54` | `DEBLEND_MINCONT 0.005`, `BACK_SIZE 64` (≈14″). Stars on a bright galaxy hold < 0.5 % of the galaxy segment's flux → merged into the galaxy segment, not deblended; galaxy structure inflates the RMS map → raised thresholds. |
| `swarp.lbc.conf:77,84` | `SUBTRACT_BACK Y`, `BACK_SIZE 128` (≈29″), not overridden in `go_swarp`. **SWarp itself removes extended galaxy light in the coadd**, independent of SExtractor. |
| `swarp.lbc.conf:13,24`; `lbcregister.py:380` | `WEIGHT_TYPE NONE`, `COMBINE_TYPE MEDIAN`, `FSCALE_KEYWORD NONE`: no weights (vignetted, noisy field corners get full weight), no transparency/exposure scaling (SCAMP's `FLXSCALE` ignored). |
| `go_scamp` | XML diagnostics are read (`astrometric_dispersion`) but never used. |
| `lbcproc.py:20` | `from lbcregister import *` (absolute, not relative) → `import LBCgo` fails when installed; tests only pass with `LBCgo/LBCgo` on `PYTHONPATH`. |
| pipeline | No weight/mask images are produced anywhere (bad columns, saturation, vignetting). |

### 3.3 Verified facts about LBC data (real headers: LBCB NGC 891 `lbcb.20141120.065509.fits`; LBCR sky flat `lbcr.20141229.132703.fits`, chips 1–2 of each)
- Raw chip: 2304 × 4608, `TRIMSEC [51:2098,1:4608]`, `BIASSEC [2099:2304,…]`.
  BITPIX 16, `SATURATE = 65536`, `GAIN = 1.75`, `RDNOISE = 12`.
- Per-chip WCS: `RA---TAN`/`DEC--TAN`, **common `CRVAL`** (= telescope pointing,
  `TELRA/TELDEC`), chip layout encoded in `CRPIX` (chip 1 `CRPIX1 = −1087`,
  chip 2 `CRPIX1 = 1035`, both `CRPIX2 = 2924`) → chip 1/2 gap = 74 px ≈ 16.6″.
  `CD` ≈ 0.2240″/px with ≈ 0.186° rotation (`PA_PNT = 360.186`). **No
  distortion terms.**
- `ccdproc.trim_image` (used in `go_overscan`) correctly shifts `CRPIX1` by the
  50-px prescan (tested: 1035 → 985); it converts `CD` to `PC` + `CDELT=1`.
- Alternate WCS `…A` keywords (AZ/EL) are stripped by `go_overscan`.
- Useful keywords: `MJD_OBS` (epoch for proper motions), `LBCFWHM` (seeing,
  arcsec), `DITHSEQ/DITHOFFX/DITHOFFY`, `INSTRUME`, `DETECTOR`,
  `FILTER`, `AIRMASS`, `EXPTIME`, `LBCCHIP1..4` (chip on/off), `DETSEC`.
- **LBCR vs LBCB layout** (chips 1–2 compared; chips 3–4 not yet seen):
  | | LBCB | LBCR |
  |---|---|---|
  | `CRPIX1` chip 1 / chip 2 | −1087 / 1035 | −1044 / 1078 |
  | chip 1 − chip 2 offset | 2122 px | 2122 px |
  | `CRPIX2` | 2924 | 2913 |
  | `DETSEC`, `TRIMSEC`, `BIASSEC` | identical | identical |
  | header scale (`CD`) | 6.222e-5° = 0.224″/px | same |
  | header rotation vs `PA_PNT` | +0.186° vs 360.186 | −0.156° vs −0.156 |
  - Same chip spacing and readout format → trimming/mask code is
    channel-independent.
  - Reference point (optical axis/rotator centre) differs by (+43, −11) px
    ≈ (9.6″, −2.5″) between cameras → distortion models must be per channel
    (as planned in §6.3.2).
  - Header scale 0.224″/px is nominal for both; the published LBC scale is
    ~0.2254″/px (+0.6 %), i.e. ~18 px (~4″) at ~2900 px from the reference
    point *before* optical distortion. Start the per-exposure solve from the
    fitted static model, not the bare header WCS (§6.3.3).
  - **Chip-gap inconsistency:** CRPIX (with the 50-px prescan) implies a
    74-px gap between chip 2's last and chip 1's first data column; `DETSEC`
    implies 49 px (cols 4451–4499), in both cameras. `DETSEC` is probably a
    nominal readout layout. Do not use `DETSEC` geometrically; the
    distortion calibration fits chip placement (the 25-px ≈ 5.6″ difference
    is inside the ±60″ offset search).
- **Header traps:**
  - `INSTRUME` is spelled inconsistently: `'LBC_BLUE'` (underscore) vs
    `'LBC-RED '` (hyphen, trailing space); the header comment says
    `'LBC-BLUE' or 'LBC-RED'`. Identify the channel with
    `LBCgo.detector.lbc_channel` (normalizes `INSTRUME`, then `DETECTOR`
    `EEV-BLUE`/`EEV-RED`, then the `lbcb`/`lbcr` filename prefix).
  - Sentinel values: `LBCFWHM = -3600.00` and `LBCBACK = -1.0` appear when the
    trackers did not measure them (seen in the LBCR sky flat). Treat
    `LBCFWHM <= 0` as missing (§6.1).
  - `TELESCOP` is `LBT-SX` (LBCB) / `LBT-DX` (LBCR).
- **Gain/read noise:** `GAIN = 1.75` e⁻/ADU and `RDNOISE = 12` e⁻ appear in
  the primary and chip headers of **both** cameras, identical on every chip
  seen → nominal values, not per-chip measurements. Published values differ:
  a per-chip table attributed to April 2010 commissioning (LBTO/Arizona LBC
  pages; *not verified*, pages unreachable from the planning session) gives
  LBCB 1.96–2.09 e⁻/ADU and LBCR 2.08–2.14 e⁻/ADU, read noise 4.8–5.3 ADU
  (≈ 10–11 e⁻); LBC papers quote ~2.02 e⁻/ADU & 5.0 ADU and ~1.75 e⁻/ADU &
  ~9 ADU (arXiv:1703.09874, arXiv:2305.10516; attribution not checked).
  Impact: in the sky-limited regime a gain error rescales all exposures of a
  chip alike (coadd weights barely change); it matters where read noise is
  not negligible (LBCB U band: at a 150-ADU sky the header and published
  values give variances ~35 % apart) and for absolute flux errors.
  Handled by `LBCgo/detector.py`: per-chip table
  `conf/lbc_detector.ecsv` (ships empty) overrides headers; header
  values are the fallback (§5.2).
- Typical observing pattern (from OB `j1419.ob`): `NDIT = 3` dither positions,
  offsets (0,0), (−40,−80), (−20,+60)″; **one exposure per filter per dither
  position** → ~3 exposures per filter per OB (repeated OBs add more).
  Dithers (20–90″) fill chip gaps but are ≪ a large galaxy.
- LBC-Blue optical distortion: pincushion, ≤ 1.75 % at field edge
  (Giallongo et al. 2008). A static average correction gives ~15 mas relative
  precision in B, V (Bellini & Bedin 2010). A pure-TAN header is therefore off
  by several arcsec (order 4–13″, depending on how the 1.75 % is defined) near
  the field edge.

### 3.4 Package landscape (verified by installing from PyPI, 2026-10-04)
| Package | Version | Notes |
|---------|---------|-------|
| `sep` | 1.4.1 | C core of SExtractor; wheels; `Background`, `extract`, `winpos`, `flux_radius`, `kron_radius`, aperture sums. LGPL. |
| `photutils` | 3.0.0 | Py ≥ 3.11, numpy ≥ 2, astropy ≥ 6.1.4. `Background2D` (masks, coverage masks), `SourceFinder`, `SourceCatalog`. v3 renamed many attributes (`n_pixels`, `x_centroid`, …); old names deprecated. |
| `drizzle` | 3.0.0 | C core; `drizzle.resample.Drizzle`, `drizzle.utils.calc_pixmap`; kernels `square, gaussian, point, turbo, lanczos2, lanczos3`; weight maps. `calc_pixmap` requires `wcs.pixel_shape`. |
| `reproject` | 0.21.0 | `reproject_and_coadd` combine ∈ {mean, sum, first, last, min, max} — **no median/clip**; `find_optimal_celestial_wcs`; `match_background`. |
| `astroquery` | 0.4.11 | Gaia archive / VizieR access. |
| `tweakwcs` | 0.9.2 | Linear-only alignment; FITS WCS with SIP supported. Optional cross-check only (pulls `gwcs`, `spherical_geometry`, `stsci.stimage`). |
| `ccdproc` | 2.5.1 | Works with numpy 2.4 / astropy 8 (existing 63 tests pass). |

Benchmarks (synthetic 2048 × 4608 chip, TAN-SIP, 4-core container):
`sep` bkg+extract 0.8 s; `sep.winpos` 0.06 s; photutils `Background2D` 1.0 s;
photutils FFT-convolve + `SourceFinder` w/ deblend 15.7 s; `reproject_interp`
18.5 s (11.9 s with `roundtrip_coords=False`, 5.4 s with `parallel=4`);
`reproject_adaptive` 34 s; drizzle `calc_pixmap` 8.0 s + `add_image` 1.7 s.
**WCS evaluation dominates resampling cost** → compute pixel maps on a coarse
grid and interpolate (§6.4).

---

## 4. Target architecture

```
LBCgo/
  lbcproc.py            (existing; fix relative import; add weight/mask output)
  masks.py              bad-pixel / saturation / vignetting masks, weight maps
                        (done; top level because lbcproc produces them)
  lbcregister.py        (existing astromatic back-end, improved in Phase 0;
                         go_register becomes a dispatcher)
  register/             (new, in-process back-end)
    __init__.py
    config.py           dataclasses with all tunables + defaults
    detect.py           sep-based detection for alignment
    sky.py              science-sky models (default + extended-target mode)
    refcat.py           Gaia DR3 query, cache, epoch propagation, quality cuts
    match.py            offset search + KD-tree matching
    distortion.py       load/apply static distortion model → per-chip TAN-SIP
    astrometry.py       per-exposure joint linear solve, header writing, QA
    calibrate.py        one-off fitting of the static distortion model
    photscale.py        per-exposure relative flux scaling
    coadd.py            output grid, pixmaps, drizzle resampling, clipped combine
    catalog.py          final catalogs (photutils SourceCatalog)
    qa.py               metrics, residual maps, comparison vs astromatic
  conf/distortion/
    lbcb_<version>.fits one HDU per chip (TAN-SIP header) + metadata table
    lbcr_<version>.fits
```

Public API (backwards compatible):

```python
go_register(filter_directories, lbc_chips=True,
            backend='native',            # or 'astromatic'
            do_detect=True, do_astrometry=True, do_coadd=True,
            astrometric_catalog='GAIA-DR3',
            extended_target=None,        # None | True | dict(ra, dec, a, b, pa)
            combine='clipped_mean',      # or 'median', 'mean'
            refcat_file=None,            # offline reference catalog
            verbose=True)
```
Keep `do_sextractor/do_scamp/do_swarp` as deprecated aliases mapping onto
`do_detect/do_astrometry/do_coadd`. `lbcgo()` passes `backend` through.

Intermediate products per chip file `<base>_<chip>.fits`:
`<base>_<chip>.mask.fits` (uint8 bitmask), `<base>_<chip>.weight.fits`,
`<base>_<chip>.srccat.fits` (alignment catalog, astropy Table),
updated WCS written into the chip header (original TAN WCS preserved as
alternate WCS `O`), `<filterdir>/astrometry_qa.ecsv`.

---

## 5. Phase 0 — astromatic path at its best + validation harness

Purpose: a fair benchmark and an immediate improvement for users.

### 5.1 Housekeeping
- [x] `lbcproc.py:20` → `from .lbcregister import *`; confirm `import LBCgo`
      works without `PYTHONPATH` hacks; run tests.
- [x] `pyproject.toml`: `requires-python = ">=3.11"`, `numpy>=2`,
      `astropy>=6.1.4`, `ccdproc>=2.5`; add `scipy`. (New deps added in later
      phases.) Update README/`docs/installation.rst`.
- [ ] Accept both astromatic executable names: Ubuntu/Debian packages
      install `source-extractor` (2.28.0) and `SWarp` (2.41.5), not `sex`
      and `swarp` (verified with apt on Ubuntu 24.04, SCAMP 2.10.0 is
      `scamp`). LBCgo checks/calls only `sex`, `scamp`, `swarp`
      (`lbcproc.check_external_dependencies`, `lbcregister.go_sextractor`,
      `go_swarp`), so an apt install is reported as missing even though the
      README suggests `apt-get`. Check the Homebrew/conda-forge names too.

### 5.2 Weight and mask maps (shared by both back-ends)
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
- [ ] Tune `badpix_threshold`/`vignette_threshold` on real LBCB/LBCR flats
      (defaults 0.2 / 0.5 are untested on real data).
- [x] Gain/read-noise source for the weights: `LBCgo/detector.py`. Lookup
      order: per-chip row of `conf/lbc_detector.ecsv` (channel, chip, MJD
      validity range) → `GAIN`/`RDNOISE` header keywords → nominal defaults.
      Weight headers record `GAINSRC` (`table`/`header`/`default`).
      Photon-transfer measurement `measure_gain_rdnoise_files(flat1, flat2,
      bias1, bias2)` handles unequal flat levels.
- [ ] Measure gain/read noise per chip for LBCB and LBCR from real bias and
      flat pairs (several epochs; run locally) and populate
      `conf/lbc_detector.ecsv` with validity ranges.

### 5.3 SExtractor improvements
- [ ] Pass `-WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE <weight>`; `-FLAG_IMAGE`
      from the mask.
- [ ] Alignment run (default): `-BACK_SIZE 32 -BACK_FILTERSIZE 3`,
      `-DEBLEND_MINCONT 1e-4`, `-DETECT_THRESH 5`, drop the odd
      `-ANALYSIS_THRESH 8`. Make all of these function arguments.
- [ ] Extended-target mode: provide the ellipse-masked, high-passed image
      (from §6.2) as `-BACK_TYPE MANUAL -BACK_VALUE 0` input, or a
      `CHECKIMAGE`-free variant using a precomputed background file.

### 5.4 SCAMP: one joint run per filter directory
SCAMP's focal-plane modes need one catalog per **exposure** with one
extension per chip.
- [ ] `merge_ldac(chip_cats) -> exposure_cat`: concatenate the
      (`LDAC_IMHEAD`, `LDAC_OBJECTS`) HDU pairs of the 4 chip catalogs of an
      exposure, in chip order (astropy `fits`).
- [ ] Run SCAMP **once** on all exposure catalogs of the filter directory:
      iteration 1 `MOSAIC_TYPE LOOSE`, later `FIX_FOCALPLANE`;
      `STABILITY_TYPE INSTRUMENT`; `ASTRINSTRU_KEY FILTER` (drop CFHT's
      `QRUNID`); restore the `-MOSAIC_TYPE` flag; delete the no-op
      `replace()`; keep `DISTORT_DEGREES 3`.
- [ ] `split_head(exposure_head) -> chip heads`: split the multi-section
      `.head` (sections separated by `END`) into `<base>_<chip>.head`.
- [ ] Parse the SCAMP XML (`AstromSigma_Internal`, `AstromSigma_Reference`,
      `XY_Contrast`, `AstromNDets_Reference`) per exposure/chip into
      `astrometry_qa.ecsv`; flag fits exceeding thresholds (implements the
      "auto-identify bad astrometric fits" ToDo).
- [ ] Check whether the installed SCAMP version applies Gaia proper motions
      to the observation epoch; record the version in QA output.

### 5.5 SWarp improvements
- [ ] `-SUBTRACT_BACK N` (sky handled by `register/sky.py`, §6.2, applied
      to chip images before SWarp) — or, if run standalone, `BACK_SIZE ≥ 1024`.
- [ ] `-WEIGHT_TYPE MAP_WEIGHT` with the §5.2 weights; `-FSCALE_KEYWORD
      FLXSCALE` (from SCAMP `.head`); `-COMBINE_TYPE CLIPPED` (Gruen et al.
      2014) as default, `MEDIAN` optional.

### 5.6 Validation harness (`register/qa.py`, used by every later phase)
Datasets (PI to provide paths; see §10):
| ID | Content | Purpose |
|----|---------|---------|
| V1 | Many (≥ 20–50) moderately rich, galaxy-free LBCB exposures; selection criteria in §6.3.2 "Calibration data" | distortion calibration + astrometry accuracy |
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

Acceptance for Phase 0: improved astromatic path runs end-to-end on V1–V5;
QA table produced; numbers recorded in `docs/planning/baseline_results.md`.

---

## 6. Phases 1–4 — native back-end

### 6.1 Phase 1a — detection for alignment (`register/detect.py`)
Input: chip image (float32), mask, weight; seeing from `LBCFWHM`
(fallback 1.0″ when missing **or ≤ 0**: the header uses −3600 as a
"not measured" sentinel) → FWHM_px = LBCFWHM / 0.224.
1. Detection background: `sep.Background(data, mask=mask, bw=32, bh=32,
   fw=3, fh=3)`. In extended mode additionally subtract a median-filtered
   image (box ≈ 5 × FWHM_px, rounded to odd; `scipy.ndimage.median_filter`
   on a 2×-binned image for speed, then upsample) so that smooth galaxy light
   is flattened before detection.
2. `sep.extract(data_sub, thresh=5, err=bkg.rms(), mask=mask, minarea=5,
   filter_kernel=<Gaussian, FWHM_px>, deblend_cont=1e-4)`.
3. Centroids: `sep.winpos(data_sub, x, y, sig)` with
   `sig = 2 × flux_radius(0.5) / 2.355` (SExtractor XWIN convention);
   positional errors from `winpos` flags/SNR.
4. Point-source selection: `flag == 0`, no masked pixel within 2 FWHM,
   peak < 0.8 × saturation (in flattened units, using the mask), SNR > 20,
   FWHM within ±30 % of the clipped mode, ellipticity < 0.3.
5. Output table: x, y (0-based, documented), errors, flux, fwhm, flags, snr.

Tests: synthetic chip with injected Gaussian stars on a smooth Sérsic-like
"galaxy" + sky gradient; require recovery of > 95 % of injected stars with
SNR > 20 inside and outside the galaxy, centroid error < 0.02 px at SNR 100.

### 6.2 Phase 1b — science sky (`register/sky.py`)
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

References: Watkins et al. 2024 (masking + parametric modelling vs dithered
stacking); Borlaff et al. 2019 (over-subtraction of extended outskirts);
Trujillo & Fliri 2016; Akhlaghi & Ichikawa 2015.

### 6.3 Phase 2 — astrometry

#### 6.3.1 Reference catalog (`register/refcat.py`)
- Query Gaia DR3 (`astroquery.gaia` ADQL; fallback VizieR `I/355/gaiadr3`)
  for a cone covering the dithered mosaic (radius ≈ 0.3°). Columns:
  `source_id, ra, dec, ra_error, dec_error, pmra, pmdec, parallax,
  phot_g_mean_mag, ruwe, astrometric_params_solved`.
- Cache to `<target>/refcat_gaiadr3_<hash>.fits` keyed on center/radius;
  `refcat_file=` allows fully offline operation.
- Propagate from J2016.0 to `MJD_OBS` with
  `SkyCoord(...).apply_space_motion(new_obstime=...)` for 5/6-parameter
  sources; 2-parameter sources kept with inflated errors (or dropped).
- Quality cuts: `ruwe < 1.4`, `astrometric_params_solved ∈ {31, 95}`,
  magnitude range chosen to avoid saturation (derived per exposure from the
  matched instrumental magnitudes; default G ∈ [16, 21]).

#### 6.3.2 Static distortion model (`register/distortion.py`, `calibrate.py`)
Representation: for each chip a FITS header with `RA---TAN-SIP`, a CRPIX
that places the common reference point (the raw headers' focal-plane
reference, e.g. chip 2 `CRPIX = (985, 2924)` after trimming), unit-scaled CD
(chip rotation/scale relative to the focal plane), and SIP `A/B` (order 3–4)
with inverse `AP/BP`. One file per channel and model version, with metadata:
channel, filters, valid MJD range, fit rms, number of exposures, software
version.

Calibration data (V1/V2; selection criteria). The distortion is shared by
all exposures, and each exposure adds only ~6 parameters of its own, so
the number of exposures and how well they cover the focal plane matter
more than the star density of any single field:
- **Unsaturated Gaia stars on every chip, including the corners**
  (distortion is largest at the field edge). Gaia is not the limit:
  median proper-motion errors are ~0.5 mas/yr at G = 20 (Lindegren et al.
  2021), i.e. a few mas when propagated back several years, well below
  the ~15 mas reached by Bellini & Bedin (2010); use G ≈ 16–20.5 with
  proper motions applied (§6.3.1). The bright limit is set by LBC
  saturation (exposure time and filter dependent; not known a priori):
  measure the usable G range per exposure from the saturation masks
  (§5.2) and matched Gaia magnitudes.
- **Intermediate Galactic latitude (|b| ≈ 10–30°).** High latitude gives
  too few Gaia stars per chip; the plane and cluster cores give blending
  (biased centroids) and Gaia crowding problems. Cluster outskirts are
  good (Bellini & Bedin used M67). Avoid large galaxies, nebulosity and
  very bright stars (ghosts, bleeds, halos). Rank candidate exposures by
  an ADQL count of Gaia DR3 sources with 16 < G < 20.5, RUWE < 1.4 in each
  exposure's footprint rather than by a rule of thumb.
- **Many exposures, not much shorter than ~30 s.** Atmospheric turbulence
  produces 10–30 mas astrometric errors in 30 s exposures, coherent over
  5–10′ (Bernstein et al. 2017, DECam) — comparable to the target and
  correlated across much of the LBC field. It is random from exposure to
  exposure and averages down only with numbers: aim for ≥ 20–50 exposures
  per channel (and per filter group, if filters turn out to differ).
  Shorter exposures help with saturation but add turbulence noise.
- **Spread in position angle and airmass.** Optical distortion is fixed to
  the detector; refraction and turbulence are fixed to the sky. A range of
  rotator angles separates them. Prefer low airmass; check residuals
  against airmass and Gaia BP−RP colour (differential chromatic
  refraction).
- **Dithers:** with Gaia as the absolute reference, small science dithers
  suffice. Large dithers (≳ one chip width) are a useful extra: the same
  stars fall on different chips, an internal check of chip placement
  independent of Gaia.
- **Filters, channels, epochs:** each channel separately (simultaneous
  LBCB+LBCR pointings are convenient); the most-used filters (refractive
  corrector → possible wavelength dependence; see step 5); several epochs,
  bracketing any known hardware interventions (§11 item 3).
- **Archive sources:** photometric standard fields at moderate latitude
  (Landolt/Stetson; short exposures, many filters and epochs); science
  programs on low-latitude targets without large galaxies; cluster
  outskirts such as M67 (also allows the step-6 comparison).
- **Hold out a test subset** of exposures (different nights) for the
  polynomial-order choice (step 4) and residual maps.

Calibration fit (`calibrate.py`, run by the PI on V1/V2):
1. For each calibration exposure: detect (§6.1), initial WCS from header,
   coarse offset search (§6.3.3), match to Gaia.
2. Unknowns: per-exposure 6-parameter linear transform (tangent point,
   scale, rotation, skew) + global per-chip polynomial distortion (shared by
   all exposures). Gaia positions are projected to the tangent plane of each
   exposure.
3. Solve by alternating least squares (fix distortion → fit linear terms per
   exposure; fix linear terms → fit distortion on all exposures), with
   sigma-clipping, until parameter changes < 1 mas; or a single sparse
   `scipy.optimize.least_squares` with `jac_sparsity`. Degeneracies: fix the
   distortion's constant and linear terms on a reference chip (they are
   absorbed by the per-exposure linear transform).
4. Choose polynomial order by residual rms vs order on held-out exposures
   (expect 3 or 4).
5. Stability tests: fit separately per filter and per observing season;
   compare displacement maps. Merge filters/epochs whose maps differ by
   < 10 mas rms; otherwise ship separate models with validity ranges.
6. LBCB sanity check: compare the displacement field with Bellini & Bedin
   (2010) for B/V.
7. Write the per-chip TAN-SIP headers; verify forward/inverse SIP
   round-trip < 0.01 px over each chip.

#### 6.3.3 Per-exposure solve (`register/astrometry.py`, `match.py`)
For each exposure (all available chips together):
1. Initial WCS per chip = static model + header pointing/rotation
   (`CRVAL` from header, rotation from header `PC`/`PA_PNT`). Use the static
   model's scale and chip placement, not the header's nominal 0.224″/px or
   `DETSEC` (§3.3).
2. Coarse offset: 2-D histogram (or cross-correlation) of all pairwise
   (detected − Gaia) tangent-plane offsets within ±60″; take the peak.
   Rotation is known to < 0.1° from the header; if the peak is weak, search
   rotation ±0.5° in steps.
3. Match with `scipy.spatial.cKDTree`, radius 2″ → 1″ → 0.5″ across
   iterations; reject ambiguous matches.
4. Fit the 6-parameter linear correction jointly over all chips (weighted
   least squares, 3σ clipping, 3 iterations). Optional per-chip shift terms
   (off by default; enable if the QA residual maps show chip offsets).
5. Internal refinement (default on): match all good sources (not just Gaia)
   between exposures; re-fit per-exposure linear terms minimising internal
   residuals with the Gaia solution as a prior. This is what achieves the
   ≤ 0.05″ internal target where Gaia density is low.
6. Write the final TAN-SIP WCS into each chip header (preserve the original
   as alternate WCS `O`); append QA row: N matches per chip, rms vs Gaia,
   internal rms, offset/rotation found.
7. Failure handling: fewer than 15 Gaia matches over the exposure → flag
   and fall back to internal alignment onto the best exposure; fewer than 5
   on a chip with others fine → chip inherits the exposure solution (that is
   the point of the joint model).

Tests (no network): synthetic catalogs generated from a known distortion
+ known per-exposure linear terms + noise; recover distortion to < 2 mas
and linear terms to < 1 mas; matching robust to 30″ pointing offsets and
30 % spurious detections. Mock `refcat` in tests.

### 6.4 Phase 3 — coadd (`register/coadd.py`, `photscale.py`)
1. Output grid: TAN, north up, pixel scale 0.224″ (configurable), bounds
   from all chip footprints (`reproject.mosaicking.find_optimal_celestial_wcs`
   with `resolution` fixed, or own bounding-box code).
2. Flux scaling: per exposure, relative zero point from the median
   magnitude difference of matched stars vs a reference exposure (highest
   transparency); store `FLXSCALE` in headers. Absolute calibration is out of
   scope (later: Gaia XP synthetic photometry, separate ToDo).
3. Pixel maps: evaluate the input→output mapping (chip WCS → sky → output
   pixel) on a coarse grid (every 32 px) and interpolate bicubically; require
   max interpolation error < 0.01 px (test). Build `drizzle` pixmap arrays
   from this.
4. Resample each exposure (sky-subtracted, scaled) with
   `drizzle.resample.Drizzle` (`kernel='square'`, `pixfrac=1.0` default;
   `lanczos3` optional) using the §5.2 weight maps; write one resampled
   science + weight layer per exposure to disk (`np.memmap` or zarr).
   Expected size: ~7000 × 7000 float32 ≈ 200 MB per layer.
5. Combine in row blocks: weighted mean with iterative σ-clipping
   (3σ, ≤ 3 iterations, minimum 3 inputs else plain weighted mean);
   `combine='median'` → weighted median option.
6. Outputs (names unchanged from SWarp path): `<object>.<filter>.mos.fits`,
   `.mos.weight.fits`, plus `.mos.nexp.fits` (number of contributing
   exposures) and `.mos.exptime.fits`. Copy keywords as `go_swarp` does now
   (OBJECT, …, TIME-OBS); write exposure-time-weighted AIRMASS; record the
   output flux unit (ADU scaled to the reference exposure, matching current
   behaviour).
7. Cosmic rays: handled by clipping when ≥ 3 exposures; for fewer, run
   `ccdproc.cosmicray_lacosmic` per chip (optional flag).

Performance target: full V3 dataset (12 exposures × 4 chips) coadded in
< 10 min on 4 cores, peak memory < 8 GB.

### 6.5 Phase 4 — science catalogs (`register/catalog.py`)
- Detection on the coadd with the weight map (`photutils.segmentation`;
  `SourceFinder` with deblending), using the science sky (not the
  detection background).
- `SourceCatalog` columns: positions (pixel + sky), windowed centroids,
  Kron and fixed-aperture fluxes with errors from the weight map, FWHM,
  ellipticity, flags.
- Optional later: SourceXtractor++ back-end for model-fitting photometry.

### 6.6 Phase 5 — switch default & clean-up
- [ ] Native back-end becomes default only after it meets or beats the
      Phase 0 numbers on V1–V5 (§7).
- [ ] Update README, `docs/pipeline.rst`, `docs/installation.rst`
      (astromatic tools become optional); add this plan's outcome to docs.
- [ ] Resolve ToDo items covered: "Auto-identify bad astrometric fits",
      "Alignment in presence of extended sources", "Background estimation in
      presence of extended sources", "Consider whether MONTAGE is better…"
      (answered: not needed), "Update to use SourceExtractor++" (D7).

---

## 7. Acceptance criteria (native vs improved astromatic, on V1–V5)

| Quantity | Requirement |
|----------|-------------|
| Internal registration rms | ≤ 0.05″ (goal ≤ 0.03″), and ≤ astromatic value |
| Absolute rms vs Gaia | < 0.1″, and ≤ astromatic value + 5 mas |
| Systematic residual map | no coherent pattern > 10 mas over the focal plane |
| Coadd star flux vs SWarp | agree to < 0.5 % (median), no trend with position |
| Coadd FWHM | ≤ 1.03 × median input FWHM |
| NGC 891 (V4, extended mode) | no negative moat; sky in masked annulus consistent with 0 within 1σ of the pixel-noise-limited mean; ≥ 90 % of Gaia stars projected onto the galaxy body detected and matched |
| Missing chip (V5) | runs without error; QA flags absent chip |
| Install | `pip install lbcgo` sufficient for the native path (no compilers, no conda) |

---

## 8. Dependencies after migration

Required: `numpy>=2`, `scipy>=1.13`, `astropy>=6.1.4`, `ccdproc>=2.5`,
`sep>=1.4`, `photutils>=3.0`, `drizzle>=3.0`, `astroquery>=0.4.11`.
Optional: `reproject>=0.21` (grid helper / cross-check), `tweakwcs`
(cross-check), astromatic binaries (`backend='astromatic'`).
All required packages publish binary wheels for Linux/macOS/Windows
(verified for Linux x86_64, CPython 3.11).

Licensing note: `sep` is LGPLv3+ (dynamic use from an MIT package is fine;
note it in the README).

---

## 9. Testing strategy

- Unit tests: synthetic data only, no network (mock `refcat`), fast
  (< 30 s total). Follow existing `tests/conftest.py` style.
- Real-data regression tests: marked `@pytest.mark.realdata`, skipped unless
  `LBCGO_TESTDATA` points at the V1–V5 directories; they run the QA harness
  and assert §7 thresholds.
- Keep the existing astromatic tests (they mock `Popen`); add tests for
  `merge_ldac`, `split_head`, XML QA parsing.
- photutils 3 uses new attribute names (`n_pixels`, `x_centroid`, …):
  use the new names from the start; set `photutils.future_column_names = True`.

---

## 10. Execution environment

- **Calibration (§6.3.2) and real-data validation (§5.6, §7) should run
  locally** (Claude Code on the PI's machine, or the PI running scripts):
  - data volume: a raw LBC MEF is ~85 MB (4 × 2304 × 4608 × 16 bit);
    calibrated float32 chips ~150 MB per exposure (4 × 2048 × 4608 × 4 bytes); V1–V5 plus intermediates
    plausibly tens of GB;
  - the cloud session used to write this plan had ~30 GB free disk, 4 cores,
    15 GB RAM, and its network policy **blocked** the Gaia archive
    (`gea.esac.esa.int`), VizieR (`vizier.cds.unistra.fr`) and the LBT
    archives; SCAMP's `REF_SERVER vizier.unistra.fr` would also be blocked.
- Cloud sessions are fine for code written against synthetic tests
  (Phases 0–4 unit tests). If cloud is used for real data, the environment's
  network policy must allow the Gaia/VizieR hosts and the data must be
  fetched from a reachable location.
- Astromatic binaries for Phase 0 benchmarking: install locally
  (`conda install -c conda-forge astromatic-source-extractor
  astromatic-scamp astromatic-swarp` or system packages).

---

## 11. Open items for the PI

1. Paths/IDs of the V1–V6 datasets (§5.6).
2. ~~LBCR chip layout~~: resolved for chips 1–2 (§3.3): same CRPIX scheme and
   spacing, reference point offset (+43, −11) px. Still to see: chips 3–4 of
   each camera (orientation of chip 4).
3. Any known LBC hardware changes (detector/corrector swaps) that should
   bound distortion-model validity ranges.
4. Preferred coadd flux unit (current: ADU scaled to reference exposure;
   alternative: ADU/s or e⁻/s).
5. Bias and flat pairs (per channel, several epochs) for the gain/read-noise
   table (§5.2).

---

## 12. References

- Bellini, A. & Bedin, L. R. 2010, A&A, 517, A34 — LBC-Blue geometric
  distortion correction, ~15 mas. https://www.aanda.org/10.1051/0004-6361/200913783
- Giallongo, E. et al. 2008, A&A, 482, 349 — LBC-Blue performance; distortion
  ≤ 1.75 %. https://arxiv.org/abs/0801.1474
- Bertin, E. & Arnouts, S. 1996, A&AS, 117, 393 — SExtractor.
- Bertin, E. 2006, ASP Conf. Ser., 351, 112 — SCAMP.
- Bertin, E. et al. 2002, ASP Conf. Ser., 281, 228 — SWarp.
- Gruen, D., Seitz, S. & Bernstein, G. M. 2014, PASP, 126, 158 — clipped-mean
  stacking in SWarp.
- Bertin, E. et al. 2020, ASP Conf. Ser., 527, 461 — SourceXtractor++.
  https://aspbooks.org/custom/publications/paper/527-0461.html
- Barbary, K. 2016, JOSS, 1(6), 58 — SEP.
- Bradley, L. et al., photutils (Zenodo; cite the version used).
- Fruchter, A. S. & Hook, R. N. 2002, PASP, 114, 144 — Drizzle.
- Shupe, D. L. et al. 2005, ASP Conf. Ser., 347, 491 — SIP convention.
- Ginsburg, A. et al. 2019, AJ, 157, 98 — astroquery.
- Gaia Collaboration, Vallenari, A. et al. 2023, A&A, 674, A1 — Gaia DR3.
- Watkins, A. E. et al. 2024, MNRAS, 528, 4289 — sky subtraction strategies
  in the LSB regime. https://arxiv.org/abs/2401.12297
- Borlaff, A. et al. 2019, A&A, 621, A133 — missing light of the HUDF
  (over-subtraction of extended outskirts).
- Trujillo, I. & Fliri, J. 2016, ApJ, 823, 123 — LSB imaging/sky treatment.
- Akhlaghi, M. & Ichikawa, T. 2015, ApJS, 220, 1 — NoiseChisel.
- Lindegren, L. et al. 2021, A&A, 649, A2 — Gaia EDR3 astrometric solution
  (uncertainties vs magnitude). https://www.aanda.org/10.1051/0004-6361/202039709
- Bernstein, G. M. et al. 2017, PASP, 129, 074503 — DECam astrometric
  calibration; atmospheric turbulence residuals. https://arxiv.org/abs/1703.01679
