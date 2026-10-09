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
    correct for every chip (verification item in §5.1). Chip 4's rotation
    lives in its WCS (CD matrix); the distortion model (§6.3.2) must take
    its orientation from the header, not assume it matches chips 1–3.
  - The two correctors were built with "the same focal plane scale and even
    the geometrical distortions ... forced to be the same" (§3.1). But the
    red corrector is BK7 (blue: fused silica) and has a 10 % larger field of
    view, chosen "to remove the small vignetting that affected the LBCB".
    → Vignetting masks matter mainly for LBCB (§5.2); the LBCB distortion
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
  seen → nominal values, not per-chip measurements.
  **Measured LBCB values** (Giallongo et al. 2008, Table 1; variance method
  on flat sequences, 2006 commissioning):
  | | chip 1 | chip 2 | chip 3 | chip 4 |
  |---|---|---|---|---|
  | gain (e⁻/ADU) | 1.96 | 2.09 | 2.06 | 1.98 |
  | read noise (e⁻) | 11.4 | 11.6 | 11.6 | 11.2 |
  So the header gain is 11–16 % low for LBCB, while the header read noise
  (12 e⁻) is close. The paper also gives 11 e⁻ at 500 kpix/s/ch for the
  controller. **LBCR:** no measured per-chip values in either paper;
  Speziali et al. (2008) give "< 10 e⁻ @500 Kpix/s/ch" for the red
  controller. (A per-chip table attributed to 2010 commissioning on LBTO/
  Arizona pages, seen only in search snippets, gives the same LBCB gains and
  LBCR gains of 2.08–2.14 e⁻/ADU; its LBCB read-noise values, 4.8–5.2 ADU ≈
  10 e⁻, do not match Table 1, so treat that table as unverified.)
  Impact: in the sky-limited regime a gain error rescales all exposures of a
  chip alike (coadd weights barely change); it matters where read noise is
  not negligible (LBCB U band: at a 150-ADU sky the header and Table 1
  values give variances ~20–30 % apart, depending on chip) and for absolute
  flux errors.
  Handled by `LBCgo/detector.py`: per-chip table
  `conf/lbc_detector.ecsv` overrides headers; header values are the
  fallback (§5.2). PI decision (2026-10-04): seed the table with the LBCB
  Table 1 values (branch `claude/seed-lbcb-detector-table`); LBCR stays on
  header values until measured.
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
    `SATURATE` mask threshold (§5.2) is conservative.
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

- [ ] (Low priority; PI expects it to hold.) Confirm from raw headers of
      both cameras that chips 3 and 4 have the same readout layout as chips
      1–2 (`NAXIS1/2 = 2304/4608`, `TRIMSEC [51:2098,…]`, `BIASSEC
      [2099:2304,…]`), i.e. that chip 4's 90° rotation is on the sky only
      and `go_overscan`'s `overscan_axis=1` is right for all chips. Record
      chip 4's CD matrix for §6.3.2.

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
      (defaults 0.2 / 0.5 are untested on real data). Tune per camera: the
      LBCR corrector's larger field was designed to remove the vignetting
      seen in LBCB (§3.3), so expect the threshold to matter mainly for
      LBCB.
- [x] Gain/read-noise source for the weights: `LBCgo/detector.py`. Lookup
      order: per-chip row of `conf/lbc_detector.ecsv` (channel, chip, MJD
      validity range) → `GAIN`/`RDNOISE` header keywords → nominal defaults.
      Weight headers record `GAINSRC` (`table`/`header`/`default`).
      Photon-transfer measurement `measure_gain_rdnoise_files(flat1, flat2,
      bias1, bias2)` handles unequal flat levels, and measures the variances
      in 50-px blocks (`cell`): the first real run (2025-05 flats, branch
      `202505_calibration`) gave LBCB gains 6–12 % below Giallongo et al.
      Table 1, consistent with illumination differences between the two
      twilight flats of a pair (different exposure time, time, rotator
      angle), which a whole-region variance turns into a low gain (−6 % for
      a 0.5 % peak-to-peak mismatch over 1000 px in simulation; < 0.5 % with
      blocks). Pick consecutive flats of one sequence at the same rotator
      angle; avoid z/Y-band flats (fringing). `lookup_detector_params` warns
      when two matching rows share `mjd_start` (`detector_table_conflicts`
      lists them): give new rows a finite `mjd_start`.
- [x] Level dependence of the photon-transfer gain. The second real run
      (branch `202505_calibration`, commit `83d8100`, several flat pairs per
      chip at different levels) shows the apparent gain rising linearly with
      flat level: +1.8–2.1 % per 10,000 ADU on LBCB and +2.9–5.1 % on LBCR,
      with residuals of 0.1–0.4 %. On LBCB that is several times what the
      < 1 % linearity residuals of Giallongo et al. (2008) allow, and is the
      signature of the brighter-fatter effect: charge pushed into
      neighbouring pixels lowers the per-pixel variance, raises the
      apparent gain, and leaves positive nearest-neighbour covariances in
      the flat difference (Antilogus et al. 2014, JINST 9, C03048; Astier
      et al. 2019, A&A 629, A36). (An earlier version of this item called
      the LBCR CCDs "thick, deep-depletion"; none of the references read
      for this plan says so. Unverified.) A median over sets therefore depends on which
      levels were observed. Adopted method (`detector.summarize_gain_rdnoise`):
      per channel/chip, fit gain = g0 + slope × level and adopt the
      zero-level intercept g0 as the conversion gain (median if fewer than
      3 sets or no spread in level); read noise = g0 × median read noise in
      ADU (the ADU value does not depend on the gain). `measure_ptc` also
      reports `rho_x`, `rho_y` (lag-1 correlation of the flat difference,
      in 50-px blocks) and `gain_nn`, the gain with those covariances added
      back to the variance. If the brighter-fatter effect explains the
      trend, `rho` grows with level and `gain_nn` is nearly flat (verified
      on a simulation in `tests/test_detector.py`: true gain 1.75, apparent
      1.79–2.00, intercept 1.749, `gain_nn` 1.75–1.77). `gain_nn` ignores
      longer-range covariances, so it is a test, not the adopted value.
      For the weight maps the intercept is the right value: sky levels are
      low, and the per-pixel gain sets the per-pixel variance. (It is not
      always the gain for fluxes summed over pixels: see `gain_flux`
      below.)
- [x] Re-run of `202505_calibration` with the new `run.py` (commit
      `eede421`, five sets per chip; assessed 2026-10-06): nearest-neighbour
      covariances explain 75–103 % of the LBCB slope (`gain_nn` slope
      −0.05 to +0.51 %/10k ADU), i.e. brighter-fatter, but only 11–30 % of
      the LBCR slope (`gain_nn` keeps +2.5 to +4.0 %/10k ADU), although
      LBCR's ρ_y grows no faster than LBCB's. Other findings:
      (a) a lag-1 anti-correlation along x (ρ_x down to −0.019 on LBCR
      chip 2, roughly level-independent, so proportional to shot noise and
      apparently electronic; the same chip has the high read noise and
      `bias_rho_x` ≈ −0.05), which pushes `gain_nn` 1–8 % above the
      intercept; (b) the intercept moves by up to 1.1 % between a linear
      and a quadratic fit, so the systematic uncertainty (1–3 %) exceeds
      `gain0_err` (0.1–0.4 %); (c) LBCB chips 2 and 3 read noise rises
      ~5 % through the bias sequence (01:28–01:41 UT; cause unknown —
      `biascheck` is the PROPID of every LBC bias, not a special
      start-of-night set (PI, 2026-10-06); compare biases from later in
      the night); (d) LBCB read
      noise in ADU agrees with Giallongo Table 1 (RN/gain) for chip 2 and
      within 4–8 % for chips 3–4, while the gains are 16–18 % lower: the
      two measurements differ in electron scale rather than in ADC
      conversion; compare with the `GAIN` keywords of 2025 headers.
- [x] Covariance sum to separate brighter-fatter from non-linearity:
      `measure_ptc(max_lag=3)` sums the correlation coefficients of the
      flat difference over all lags with |dx|, |dy| ≤ 3 (`rho_sum`;
      `_cell_covariance_sum`, removing a plane per 50-px block and
      correcting the −3(1 + S)/n bias of each lag) and reports
      `gain_sum` = gain with var × (1 + `rho_sum`). Charge conservation
      makes `gain_sum` free of the brighter-fatter effect (within the lag
      range); non-linearity creates no covariances and survives.
      `gain_fit.ecsv` gains `rho_sum_slope_per_10k`, `gain_sum_median`,
      `gain_sum_slope_pct_per_10k` ± err (err: the larger of the residual
      and propagated errors). Simulations (`tests/test_detector.py`): with
      half the charge sharing at 2 px, the gain slope 5.6 %/10k ADU leaves
      2.4 in `gain_nn` and 0.1 ± 0.5 in `gain_sum`; with a sublinear
      response the gain slope 4.9 stays 4.9 ± 0.6 in `gain_sum`. Noise:
      ~1 % per set for a 1000 × 1000 box (48 lags), so run this test with
      the whole chip (`run.py --box 0`, ~0.3 %).
- [x] Whole-chip re-run of `202505_calibration` (`--box 0`, commit
      `a2d7108`; assessed 2026-10-06). Adding back the covariances summed
      to 3 px (`gain_sum`) removes the level dependence on every chip:
      `gain_sum` slopes −0.78 to +0.05 %/10k ADU (all sets) and −0.65 to
      +0.67 (close pairs only), against per-pixel gain slopes of +1.9–2.0
      (LBCB) and +3.2–5.0 (LBCR). So the LBCR trend is the brighter-fatter
      effect with covariances beyond lag 1 (`rho_sum` grows 0.035–0.054
      per 10k ADU on LBCR, 0.020–0.026 on LBCB), not non-linearity:
      LBCR non-linearity ≲ 0.3 % at 10k ADU and ≲ 0.7 % at 22k ADU
      (from |`gain_sum` slope| ≲ 1 %/10k ADU, apparent gain ∝ 1 + 3βN).
      No separate LBCR linearity test is needed for the gain. The
      intercept g0 changes by ≤ 0.13 % when the widely spaced pairs are
      dropped.
- [x] Two gains (`detector.OPTIONAL_TABLE_COLUMNS`). Correlations present
      at zero signal make the per-pixel gain and the flux gain differ:
      median `gain_sum` / g0 = −4.6 % on LBCR chip 1 (positive serial
      correlation at low level, falling with level: CTI-like), +2.8 % on
      LBCR chip 2 (serial anti-correlation, ρ_x ≈ −0.020: electronic),
      −1.2 to +1.2 % elsewhere. Both differences are reproduced by
      1/(1 + S) with S the summed correlation extrapolated to zero signal.
      For any linear readout kernel with weights summing to H (CTI: H = 1;
      undershoot: H < 1), the mean of an aperture sum scales as H and its
      variance as H², so mean/variance of aperture sums — `gain_sum` once
      `max_lag` covers the kernel — is the electrons per ADU of a flux.
      Hence: `gain` (per pixel, zero-level intercept) for per-pixel
      variance and weight maps; `gain_flux` (median `gain_sum`) for
      Poisson errors of source fluxes and flux→electron conversion.
      `gain_flux` is an optional column of `conf/lbc_detector.ecsv` (NaN =
      unknown; `read_detector_table` adds it to older tables;
      `lookup_gain_flux` falls back to `gain` and says so). Nothing in the
      pipeline uses it yet. Correlated read noise: `bias_rho_sum` 0.03–0.24
      (LBCR chip 2 negative), so read noise in an aperture is up to ~11 %
      above the independent-pixel value; per-pixel read noise is
      unaffected.
- [x] `gain_sum` depends on the time between the two flats: χ²/dof of
      `gain_sum` about a line 1.1–12.5; averaged over chips, pairs 41–60 s
      apart read +0.3 to +0.6 % (sets 7–8: −0.1, −0.3 %) and pairs
      143–212 s apart −0.5 to −0.8 % (lowest: set 5, which also mixes
      pa0/pa180). Presumably the twilight changes between exposures and
      leaves small-scale structure that the 48-lag sum weights heavily;
      g0 is insensitive. `summarize_gain_rdnoise(flux_max_dt=60)` uses only
      pairs ≤ 60 s apart for `gain_flux` (NaN if none), with an error
      from the scatter of those pairs; `run.py` records `flat_dt` per set
      and warns when the two flats come from different OBs (`lbcobnam`).
- [x] 2025-05-27 rows merged into `conf/lbc_detector.ecsv` (PI decision
      2026-10-07): 8 rows (LBCB + LBCR, `gain`, `rdnoise`, `gain_flux`)
      from `calibration/gain_rdnoise_lbc_202505/`, `mjd_start` = 60822,
      open-ended. They supersede the 2006 LBCB rows from MJD 60822; the
      2006 rows still apply to earlier dates and to an unknown date, and
      LBCR before 60822 uses header values. Expected uncertainty of
      `gain`: ~1 % (linear vs quadratic fit), not the 0.1–0.4 % of
      `gain0_err`.
- [x] 2014-06/07 rows merged into `conf/lbc_detector.ecsv` (PI decision
      2026-10-09): 8 rows from `calibration/gain_rdnoise_lbc_201406/`,
      `mjd_start` = 56830 (2014-06-22 UT, earliest flats used), open-ended
      in the product, superseded by the 2025 rows at 60822. LBCB agrees with
      2025 to ~1 %; LBCR uses three flat pairs only (provisional).
- [ ] Measure further epochs (the PI plans several soon): move
      `mjd_start` earlier if older data agree, add date-limited rows if
      they do not. Before date-limiting the seeded LBCB rows, compare
      with the `GAIN` keywords of the 2025 headers and the LBT/LBC
      team's current values: both gains are 0.82–0.92 × the
      2006 values while read noise in ADU agrees for chip 2, i.e. the
      electron scales differ.
- [ ] More data: flat pairs at 1–4k ADU (shorter extrapolation to zero
      level), consecutive pairs ≤ 60 s apart at the same rotator angle
      (for `gain_flux`), and biases from later in the night (read-noise
      drift on LBCB chips 2–3).
- [ ] Measure gain/read noise per chip for LBCB and LBCR from real bias and
      flat pairs (several epochs; run locally) and populate
      `conf/lbc_detector.ecsv` with validity ranges. Check the LBCB results
      against Giallongo et al. (2008) Table 1 (§3.3), which seeds the
      LBCB rows for now (open date range, source = the paper; PI decision
      2026-10-04, PR jchowk/LBCgo#6; provenance in
      `calibration/gain_rdnoise_lbcb_giallongo2008/`). Replace or
      date-limit those rows once measured values exist, and record each
      measurement run in its own `calibration/` directory (§9.1).
- [ ] (Optional, matters for low-surface-brightness work.) Electronic
      cross-talk ~3 × 10⁻⁵ (Giallongo et al. 2008): a saturated star
      imprints ~2 ADU ghosts in the other chips/channels. Either correct it
      (needs the coefficient matrix) or mask those positions in
      extended-target mode.

### 5.3 SExtractor improvements
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
      `BACK_VALUE 0`. **Pending:** the producer of that image (§6.2).
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

### 5.4 SCAMP: one joint run per filter directory
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
      high-PM stars sit at a median (−48, −35) mas. **Still to check:** that
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
  (as done in the validation script). `go_swarp` path handling was fixed in §5.5.

### 5.5 SWarp improvements
Implemented in `go_swarp` (`go_register(swarp_args=dict(...))` forwards
overrides); the packaged `swarp.lbc.conf` defaults were changed to match.
- [x] Background: `go_swarp(subtract_back=True, back_size=1024)` by default
      (standalone use, since `register/sky.py` §6.2 does not exist yet);
      `subtract_back=False` gives `-SUBTRACT_BACK N` for use once the sky is
      removed upstream or for extended targets.
- [x] `-WEIGHT_TYPE MAP_WEIGHT` with the §5.2 `<base>.weight.fits` sidecars
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
  frame). Real-data comparison belongs to the §5.6 harness.
- **Not verified:** that real SCAMP `FLXSCALE` values are sensible for LBC
  data (needs PI datasets).

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
- Ghosts (LBCB, Giallongo et al. 2008): with the U-LBC interference filter
  (header `FILTER = 'SDT_Uspec'`, so this applies to the NGC 891 U data),
  mask each bright star's ghost (ring 75 px + diffuse 200 px component,
  shifted radially outward, 2.8 % of the star's flux) before fitting the
  sky; the ~0.15 % sky ghost near the field centre is part of the sky model
  or flat, not a source. Not needed for Bessel U, B, V or G, R.

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

1. Paths/IDs of the V1–V6 datasets (§5.6). This will come later.
2. ~~LBCR chip layout~~: resolved for chips 1–2 (§3.3): same CRPIX scheme and
   spacing, reference point offset (+43, −11) px. Chip 4 is rotated 90° on
   the sky with the same readout layout (Speziali et al. 2008; PI). Left: the
   low-priority header check in §5.1. We will need to examine the headers from all four chips to verify their known locations with respect to the central chip (chip #2 in extension 2). 
3. There are no hardware changes that should affect the result over time (other than perhaps some natural drift that could change things over time).
4. Preferred coadd flux unit is ADU/s.
5. Bias and flat pairs (per channel, several epochs) for the gain/read-noise
   table (§5.2) will be collected. This is a later priority.

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
