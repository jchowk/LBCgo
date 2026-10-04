- [ ] Work out how to handle only saturated flats. (Right now this just kills the process, even if there are remaining filters w/o saturated flats.)

- [ ] Identify high background images?

- [ ] FIXPIX
- [ ] Auto-identify bad astrometric fits.

- [ ] DOCUMENTATION: co-pointing files (with one chip) will cause issues. Code assumes all files have same number of chips when processing. (ERROR: IndexError: list index out of range).

- [ ] Alignment in presence of extended sources.
- [ ] Background estimation in presence of extended sources.

- [ ] BESTstars / Gaia XRP spectra for photometric calibration.
- [ ] Update to use SourceExtractor++
- [ ] Consider whether MONTAGE is better than SCAMP+SWARP

- [ ] `go_scrub` is a no-op: `rm data/` lacks `-r`, and `rm astro*pdf` via `Popen` has no shell glob expansion. It also works relative to the cwd; should take `image_directory`. Fixing it means `lbcgo(clean=True)` will actually delete `data/`.

- [ ] `make_flatfield` normalizes on a hard-coded box (`data[500:1500, 2000:2500]`) sized for full chips; a smaller ROI/readout gives an empty or off-chip region. Derive the box from the chip shape.

- [ ] Consider Claude's notes on potential concerns:
## Design limitations documented in tests (not bugs, but worth knowing)
* `make_flatfield` has a hardcoded normalization section data[500:1500, 2000:2500] designed for full 2048×4612 LBC chips — returns empty on any smaller array. The 8 make_flatfield tests still pass (they check structure, not values); the 3 go_flatfield tests use a synthetic master flat to work around this.