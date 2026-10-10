# LBCgo registration & coaddition: migration plan

Status (2026-10-09): Phase 0 code done on synthetic tests; work order
revised by the PI (§1.1, D8): the astromatic path is made robust and gains
an extended-target mode first, the static distortion model and native
back-end come last. Written 2026-10-04 as the agreed plan.
Audience: a later implementation session (human or Claude). Read this whole
document before writing code; section 2 records decisions that are not to be
re-litigated without the PI (J. C. Howk).

This document is planning only. Companion documents:
- `docs/planning/astromatic_path_plan.md`: the Phase 0 record and the
  current work on the astromatic path (robustness, validation harness,
  extended-target mode), with checkboxes.
- `docs/detector_gain_rdnoise.md`: the per-chip gain and read-noise work
  (method, results by epoch, open measurements).
- `docs/planning/claude_handoff.md`: working notes for Claude sessions.

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

### 1.1 Order of work (PI, 2026-10-09)

| Stage | Content | Where | Status |
|---|---|---|---|
| Phase 0 | Astromatic path at its best on synthetic tests (§5) | astromatic plan §2 | done |
| A | Astromatic path robust to the known failure modes (bad header cards, missing chips, saturated/unmatched flats, epoch keywords, flux scale, bad-fit handling, …) | astromatic plan §3 | open |
| B | Validation harness (`register/qa.py`, `register/refcat.py`) and baseline numbers on V1–V6 | astromatic plan §4 | open |
| C | Extended-target mode on the astromatic path (`register/sky.py` feeding SExtractor and SWarp) | astromatic plan §5 | open |
| D | Static distortion model (§6.3.2), then the native back-end (§6) and the switch of default (§6.6) | this plan §6; to be carved out into its own plan | deferred |

The distortion model is last because it is the most time-consuming item
(calibration data selection, 20–50 exposures per channel, local runs).
The native back-end's astrometry depends on it; the shared pieces built in
B and C (`refcat.py`, `qa.py`, `sky.py`) are reused there unchanged.

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
| D8 | Order of work (2026-10-09): the astromatic path stays the working path until the native back-end exists. First make it robust to the known failure modes and able to combine extended-target fields (§1.1 stages A–C); the static distortion model and the native back-end come last. D1 stands as the end state. |

## 3. Current state (verified 2026-10-04, commit `4b9efbc`; §3.1–3.2 updated 2026-10-09)

### 3.1 Pipeline flow
`lbcgo()` → `go_overscan` → `go_flatfield` → `make_targetdirectories` →
`go_extractchips` (writes `<base>_<chip>.fits`, single-extension, one per chip,
with `.mask.fits`/`.weight.fits` sidecars; moves the MEF `*_flat.fits` into
`data/`) → `go_register(fltr_dirs, …)`.

`go_register` runs `go_sextractor` per chip file, then (default since Phase
0) one joint SCAMP run per filter directory on per-exposure merged catalogs
(`go_scamp_joint`; `scamp_joint=False` restores the old per-chip
`go_scamp`), then `go_swarp` on all chip files of the filter directory.

### 3.2 Defects / limitations found in the astromatic usage (2026-10-04)
All fixed in Phase 0 (astromatic plan §2) except where the Status column
says otherwise.

| Location | Issue | Status |
|----------|-------|--------|
| `lbcregister.py:537` | SCAMP is invoked **once per chip catalog**, so every chip is solved independently with a 3rd-order polynomial (`scamp.lbc.conf:66`, `DISTORT_DEGREES 3`). No focal-plane or cross-exposure constraint. Under-constrained on star-poor chips → likely source of the "bad astrometric fits" ToDo item. | fixed: joint SCAMP run |
| `lbcregister.py:248` | `-MOSAIC_TYPE` flag commented out; the per-iteration `mosaic_type` values are dead code. `scamp.lbc.conf:48` → `UNCHANGED`. | fixed |
| `lbcregister.py:251` | `cmd_flags.replace('INSTRUMENT','EXPOSURE')` discards its result (and the string never contains `INSTRUMENT`): `astrometric_method` is a no-op. | fixed (removed) |
| `lbcregister.py:119` | Command-line `-DETECT_THRESH 5.0 -ANALYSIS_THRESH 8.0` overrides the config's 1.5/1.5. | fixed |
| `sextractor.lbc.conf:23,54` | `DEBLEND_MINCONT 0.005`, `BACK_SIZE 64` (≈14″). Stars on a bright galaxy hold < 0.5 % of the galaxy segment's flux → merged into the galaxy segment, not deblended; galaxy structure inflates the RMS map → raised thresholds. | fixed (alignment defaults) |
| `swarp.lbc.conf:77,84` | `SUBTRACT_BACK Y`, `BACK_SIZE 128` (≈29″), not overridden in `go_swarp`. **SWarp itself removes extended galaxy light in the coadd**, independent of SExtractor. | partial: `subtract_back`, `back_size=1024`; extended mode is Stage C |
| `swarp.lbc.conf:13,24`; `lbcregister.py:380` | `WEIGHT_TYPE NONE`, `COMBINE_TYPE MEDIAN`, `FSCALE_KEYWORD NONE`: no weights (vignetted, noisy field corners get full weight), no transparency/exposure scaling (SCAMP's `FLXSCALE` ignored). | fixed: weights, `FLXSCALE`, clipped mean |
| `go_scamp` | XML diagnostics are read (`astrometric_dispersion`) but never used. | partial: QA table written; `bad` flag not acted on (Stage A9) |
| `lbcproc.py:20` | `from lbcregister import *` (absolute, not relative) → `import LBCgo` fails when installed; tests only pass with `LBCgo/LBCgo` on `PYTHONPATH`. | fixed |
| pipeline | No weight/mask images are produced anywhere (bad columns, saturation, vignetting). | fixed (`masks.py`) |

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
- **LBCR vs LBCB layout** (chips 1–2 compared from headers; chip 4 is
  rotated on the sky, see "Detector arrangement" below):
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
  - Header scale 0.224″/px is nominal for both. Giallongo et al. (2008,
    §3.1) measure LBCB at **0.2275″ ± 0.0001 at the centre** and a **median
    of 0.2254″ ± 0.0001**, with filter-to-filter variation "affecting the
    fourth decimal digit" (scale decreases outward: pincushion). The header
    is therefore 0.6–1.6 % small: of order 20–45 px (4–10″) at ~2900 px from
    the reference point, before the distortion itself. Start the
    per-exposure solve from the fitted static model, not the bare header
    WCS (§6.3.3).
  - **Chip gaps (resolved):** Giallongo et al. (2008, §2.1): "the gaps
    between the vertical chips are 1 mm, which corresponds to 74 pixels or
    16.7 arcsec"; between the vertical chips and the horizontal chip 4,
    1.03 mm (76 px, 17.2″). This matches the 74 px implied by CRPIX; the
    49 px implied by `DETSEC` is not physical. Do not use `DETSEC`
    geometrically.
  - **Optical centre:** Giallongo et al. place the geometrical field centre
    at pixel (1024, 2919) of chip 2 (coordinate convention not stated;
    compare header `CRPIX` (1035, 2924) raw or (985, 2924) after trimming).
    The distortion fit should solve for the distortion centre, starting
    from this value (§6.3.2).
- **Detector arrangement and optics** (Speziali et al. 2008, SPIE 7014,
  70144T; read in full):
  - Both cameras use the same arrangement: four E2V 42-90 CCDs, three side
    by side and **chip 4 rotated 90° on the sky** above them (paper Figs. 3,
    4), ~0.23″/px, plus two small technical chips for guiding/active optics.
  - **The rotation is on the sky only.** In the raw readout frame all four
    chips have the same layout (2304 × 4608, overscan along x), as the PI's
    raw NGC 891 display shows, so `go_overscan`'s `overscan_axis=1` is
    correct for every chip (verification item A14 in the astromatic plan). Chip 4's rotation
    lives in its WCS (CD matrix); the distortion model (§6.3.2) must take
    its orientation from the header, not assume it matches chips 1–3.
  - The two correctors were built with "the same focal plane scale and even
    the geometrical distortions ... forced to be the same" (§3.1). But the
    red corrector is BK7 (blue: fused silica) and has a 10 % larger field of
    view, chosen "to remove the small vignetting that affected the LBCB".
    → Vignetting masks matter mainly for LBCB (astromatic plan §2.2, A11); the LBCB distortion
    model is a good starting guess for LBCR but is fitted separately
    (§6.3.2).
  - Detectors are coplanar to ±13.5 µm (one pixel) without shimming.
  - The paper notes that the distortion "is visible on the raw frames as a
    background enhancement in the central part of the images, that must not
    be confused with genuine flat field effect": the **pixel solid angle
    varies across the field**. See the surface-brightness requirement in
    §6.4.
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
  seen → nominal values, not per-chip measurements. Measured values (from
  Giallongo et al. 2008, Table 1, and LBCgo's own photon-transfer
  measurements of 2010-03, 2014-06/07 and 2025-05) are in `conf/lbc_detector.ecsv`,
  read by `LBCgo/detector.py`; header values are the fallback. The record of
  that work, including an unexplained difference of the LBCB gains
  between 2010 and 2014–2025, is `docs/detector_gain_rdnoise.md`.
  Impact on registration: in the sky-limited regime a gain error rescales
  all exposures of a chip alike (coadd weights barely change); it matters
  where read noise is not negligible (LBCB U band) and for absolute flux
  errors.
- Typical observing pattern (from OB `j1419.ob`): `NDIT = 3` dither positions,
  offsets (0,0), (−40,−80), (−20,+60)″; **one exposure per filter per dither
  position** → ~3 exposures per filter per OB (repeated OBs add more).
  Dithers (20–90″) fill chip gaps but are ≪ a large galaxy.
- LBC-Blue optical distortion (Giallongo et al. 2008, §2.1, §3.1, Fig. 4):
  pincushion, "always below 1.75 % even at the edge of the field"; their
  distortion map shows 1 %, 1.5 % and 2 % contours. Their astrometric
  solution (AstromC) was close to the optical-design prediction;
  "second-order corrections vary from frame to frame because of different
  elevation, filter or position angle. These variations are however very
  small." A static average correction gives ~15 mas relative precision in
  B, V (Bellini & Bedin 2010). A pure-TAN header is off by several arcsec
  near the field edge.
- Other LBCB facts relevant here (Giallongo et al. 2008):
  - Unvignetted field 27′ diameter; 5 % light loss at the edge of the
    corrected field. The chips span 23.6′ × 25.3′, so chip corners (up to
    ~17′ from the centre) lie outside it: vignetting there is expected.
  - Flats from twilight + night sky; their flat-field illumination profile
    (Fig. 5) is "corrected for pixel scale variation across the field",
    i.e. they treat the pixel-area effect separately, as §6.4 requires.
  - Ghosts: none measurable with Bessel U, B, V or custom G, R. With the
    interference U-LBC filter, a bright star's primary ghost holds
    2.8 ± 0.7 % of its flux: a ring 75 px across plus a diffuse 200 px
    component shifted radially outward. A sky ghost adds ~0.15 % near the
    field centre. The header filter name `SDT_Uspec` (e.g. the NGC 891
    data) **is** the U-LBC filter (PI, 2026-10-04).
  - Electronic cross-talk between chips/channels: coefficients ~3 × 10⁻⁵;
    the LBC team's pipeline corrects it. LBCgo does not.
  - Linearity residual < 1 % over the full 16-bit range; full well
    > 150,000 e⁻ before blooming, above the ADC limit (65535 ADU ≈ 130,000
    e⁻ at ~2 e⁻/ADU). Saturation is therefore the ADC limit; the 0.9 ×
    `SATURATE` mask threshold (astromatic plan §2.2) is conservative.
  - Bias: they fit pre-scan and over-scan line by line; LBCgo fits a
    4th-order polynomial to the over-scan only.

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
  register/             (new, in-process back-end; refcat.py and qa.py are
                         built in stage B, sky.py in stage C, for the
                         astromatic path; the rest in stage D)
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
calibration/            (top level, not shipped in the package; §9.1)
  README.md             provenance convention
  NEW_PRODUCT_TEMPLATE.md  README skeleton for a new calibration product
  <product>_<version>/  README.md, inputs.ecsv, run.py per product
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

Moved to `docs/planning/astromatic_path_plan.md` on 2026-10-09:
- record of the Phase 0 changes (housekeeping, weight/mask maps,
  SExtractor, joint SCAMP, SWarp): astromatic plan §2 (all done);
- remaining work as Stage A (robustness to known failure modes): §3;
- the validation harness (formerly §5.6: datasets V1–V6, metrics) as
  Stage B: §4;
- the gain/read-noise items formerly in §5.2: `docs/detector_gain_rdnoise.md`.

---

## 6. Phases 1–4 — native back-end (stage D: after the astromatic stages A–C)

The static distortion model (§6.3.2) is the first and largest part of
stage D (D8); the native astrometry (§6.3.3) needs it. To be carved out
into its own planning document when the work starts.

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
Built in stage C for the astromatic path and reused here unchanged. The
design (default and extended-target modes, ghosts, references) is in
`docs/planning/astromatic_path_plan.md` §5.


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
version. Chip 4 is rotated 90° on the sky relative to chips 1–3 (§3.3):
take each chip's starting CD matrix from its raw header rather than
assuming a common orientation. Fit the distortion centre rather than fixing
it at the header `CRPIX`; start from Giallongo et al.'s optical centre,
pixel (1024, 2919) of chip 2 (§3.3). The two correctors were designed with the
same scale and distortion (Speziali et al. 2008), so the fitted LBCB model
is a good starting point for LBCR, but LBCR is fitted separately (different
glass, larger field, different reference point).

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
  (astromatic plan §2.2) and matched Gaia magnitudes.
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
  (Landolt/Stetson; short exposures, many filters and epochs; LBCB
  commissioning used SA98 and SA113); science programs on low-latitude
  targets without large galaxies; cluster outskirts such as M67 (also
  allows the step-6 comparison). LBCB commissioning observed NGC 7789,
  NGC 2419 and M67 specifically to study the PSF across the field and the
  astrometric distortion (Giallongo et al. 2008, §3); if those data are
  in your archive they are natural V1 candidates.
- **Hold out a test subset** of exposures (different nights) for the
  polynomial-order choice (step 4) and residual maps.

Calibration inputs (form of the data):
- **Raw** multi-extension `lbcb.*` / `lbcr.*` frames as delivered by the
  archive (~85 MB each). `calibrate.py` reads them directly and does
  overscan subtraction and trimming itself, reusing `go_overscan`'s
  per-chip code (factor it out rather than duplicating it), so positions are
  in the same trimmed chip coordinates the pipeline uses. It also writes the
  saturation mask: saturated stars give biased centroids.
- **No bias frames:** the overscan removes the bias level; a separate bias
  frame does not change centroids.
- **Flat optional**, used only for the bad-pixel mask (`masks.flat_mask`)
  and for more uniform detection depth in the vignetted corners. Any flat
  from a nearby epoch will do; flat fielding does not move centroids
  measurably (the pixel-to-pixel term is ~0.01 px for a bright star,
  estimated, below the 10–30 mas per-exposure turbulence), and the
  pixel-area effect changes fluxes, not positions.
- Header keywords used: `MJD_OBS` (proper-motion epoch), `CRVAL`/`CRPIX`/
  `CD` (starting WCS), `PA_PNT`, `ROTANGLE`, `AIRMASS`, `FILTER`,
  `INSTRUME`/`DETECTOR`. `go_overscan` strips `ROTANGLE`/`PARANGLE` from
  the chip headers but keeps them in the primary header.
- Outputs: the model file in `LBCgo/conf/distortion/` (with a `PROVENANCE`
  keyword naming its `calibration/` directory) and that directory's
  README, `inputs.ecsv` and `run.py` (§9.1). Matched star catalogs
  (exposure, chip, x, y, Gaia `source_id`, MJD) go to a Zenodo deposit or
  release asset so the fit can be redone without the raw frames.

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
   compare displacement maps after removing a linear term (a pure scale
   change between filters is absorbed by the per-exposure linear fit).
   Merge filters/epochs whose maps differ by < 10 mas rms; otherwise ship
   separate models with validity ranges. Expectation: Giallongo et al.
   report filter-to-filter variation "of the order of 0.01 %,
   corresponding to about 1 pixel at the edge of the FoV" (the two numbers
   do not quite agree: 0.01 % is ~0.3–0.5 px at the edge); if it is a
   change of scale it is absorbed, but if it changes the shape of the
   distortion it is ≫ 10 mas and per-filter models are needed.
6. LBCB sanity check: compare the displacement field with Bellini & Bedin
   (2010) for B/V. Also compare the LBCB and LBCR displacement fields:
   they were designed to be the same, so large differences beyond chip
   placement point to a fitting problem or a real optical difference worth
   understanding.
7. Produce a relative pixel-area map per chip (determinant of the
   distortion Jacobian, normalized to the reference point) and ship it with
   the model. §6.4 needs it for photometry on single, unresampled frames
   (not for resampling).
8. Write the per-chip TAN-SIP headers; verify forward/inverse SIP
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
   Optional per-exposure 2nd-order terms (off by default): Giallongo et al.
   (2008) found second-order corrections varying with elevation, filter and
   position angle ("very small"). Enable only if residual maps from V1/V2
   show coherent per-exposure quadratic patterns above the turbulence
   level, and regularize them towards zero.
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
**Surface brightness (requirement).** The distortion changes the pixel solid
angle across the field (pincushion: pixels cover less sky towards the
edge; Speziali et al. 2008 note the resulting central "background
enhancement" in raw frames). A sky flat divides this out, so after
`go_flatfield` each chip is in **surface-brightness** units: the sky is flat,
but the summed counts of a star on a single frame are biased by the local
pixel area (percent level near the edge, given the ≤ 1.75 % distortion
quoted by Giallongo et al. 2008; the exact area variation comes from the
fitted model). Therefore:
- Resample chips as surface brightness, i.e. without a position-dependent
  pixel-area rescaling. `drizzle` 3.0 with its defaults (`iscale=1`,
  `in_units='cps'`) already does this; verified 2026-10-04: a uniform input
  stays uniform under a strong SIP distortion and under a 2× change of
  pixel scale, and the output sum of a point source scales with the local
  Jacobian (27.92 vs 25 × J = 27.92 at a corner). So drizzle the
  sky-flattened chips directly. **Do not** multiply them by a pixel-area
  map first: that would apply the area correction twice. Do not set
  `iscale`/`pixel_scale_ratio` to anything but a constant.
- Required test: a synthetic star field with known fluxes and a known
  distortion, "sky-flattened" (divided by its relative pixel area), must
  come out of the coadd with fluxes independent of position to < 0.2 %.
- Any photometry on single, unresampled frames (step 2, QA) must
  **multiply** the summed flux by the relative pixel area a at the source
  position (§6.3.2 step 7): the flattened frame holds raw/a, so a star's
  sum is S/a. Or measure on the resampled images instead.

1. Output grid: TAN, north up, pixel scale 0.224″ (configurable), bounds
   from all chip footprints (`reproject.mosaicking.find_optimal_celestial_wcs`
   with `resolution` fixed, or own bounding-box code).
2. Flux scaling: per exposure, relative zero point from the median
   magnitude difference of matched stars vs a reference exposure (highest
   transparency), with single-frame fluxes corrected by the pixel-area map
   (see "Surface brightness" above; a dithered star sits at different
   field positions in different exposures); store `FLXSCALE` in headers. Absolute calibration is out of
   scope (later: Gaia XP synthetic photometry, separate ToDo).
3. Pixel maps: evaluate the input→output mapping (chip WCS → sky → output
   pixel) on a coarse grid (every 32 px) and interpolate bicubically; require
   max interpolation error < 0.01 px (test). Build `drizzle` pixmap arrays
   from this.
4. Resample each exposure (sky-subtracted, scaled; surface-brightness
   units, see above) with `drizzle.resample.Drizzle` (`kernel='square'`, `pixfrac=1.0` default;
   `lanczos3` optional) using the weight maps (astromatic plan §2.2); write one resampled
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
      stage B baseline numbers on V1–V5 (§7; `docs/planning/baseline_results.md`).
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

### 9.1 Calibration provenance

Convention set out in `calibration/README.md` (PI decision 2026-10-04):
- Code that derives a calibration lives in the package, with tests
  (`LBCgo/detector.py` now; `LBCgo/register/calibrate.py` later).
- Products the pipeline reads live in `LBCgo/conf/` and point back to
  their provenance (`source` column; `PROVENANCE` header keyword).
- How each product was made lives in a top-level `calibration/<product>/`
  directory: README (date, who, LBCgo commit, product checksums, results),
  `inputs.ecsv` (archive filenames / `OBS_ID`s, MJD, channel, filter,
  role) and `run.py` (orchestration only; data paths from an environment
  variable). New versions get new directories.
- Raw frames stay out of git. Intermediate catalogs small enough (≲ 5 MB
  compressed) can be committed; otherwise Zenodo or a release asset.
- `.gitignore` ignores `temp*` and `_*`: do not name files in
  `calibration/` "template…" or with a leading underscore.

---

## 10. Execution environment

- **Calibration (§6.3.2) and real-data validation (astromatic plan §4, §7) should run
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

1. Paths/IDs of the V1–V6 datasets (astromatic plan §4). Needed now for
   stage B (PI agreed 2026-10-09 to run the baseline before native code).
2. ~~LBCR chip layout~~: resolved for chips 1–2 (§3.3): same CRPIX scheme and
   spacing, reference point offset (+43, −11) px. Chip 4 is rotated 90° on
   the sky with the same readout layout (Speziali et al. 2008; PI). Left: the
   low-priority header check (astromatic plan A14). We will need to examine the headers from all four chips to verify their known locations with respect to the central chip (chip #2 in extension 2). 
3. There are no hardware changes that should affect the result over time (other than perhaps some natural drift that could change things over time).
   Caveat (2026-10-10): the LBCB flat-based gains of chips 2–4 are 8–12 %
   higher in 2010 than in 2014 and 2025 (which agree to ~1 %), with no
   known hardware or controller change. It is either a change between
   2010-03 and 2014-06 or a biased 2010 measurement; the PI is measuring
   epochs in that window (`docs/detector_gain_rdnoise.md` §3, §5). This concerns the electronics;
   it says nothing yet about the optics/distortion.
4. Preferred coadd flux unit is ADU/s.
5. Bias and flat pairs (per channel, several epochs) for the gain/read-noise
   table will be collected; in progress (`docs/detector_gain_rdnoise.md` §5).

---

## 12. References

- Bellini, A. & Bedin, L. R. 2010, A&A, 517, A34 — LBC-Blue geometric
  distortion correction, ~15 mas. https://www.aanda.org/10.1051/0004-6361/200913783
- Giallongo, E. et al. 2008, A&A, 482, 349 — LBC-Blue performance: distortion
  ≤ 1.75 %, pixel scale 0.2275″ centre / 0.2254″ median, chip gaps 74/76 px,
  per-chip gain and read noise (Table 1), ghosts, cross-talk, linearity.
  https://www.aanda.org/10.1051/0004-6361:20078402 (read in full)
- Speziali, R. et al. 2008, Proc. SPIE 7014, 70144T — LBC description and
  performance, both cameras: detector arrangement, correctors designed with
  equal scale and distortion, LBCR vignetting, read noise, pixel-area
  background effect. doi:10.1117/12.790132 (read in full)
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
