- [ ] DOCUMENTATION
 
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

- [ ] Consider Claude's notes on potential concerns:
## Design limitations documented in tests (not bugs, but worth knowing)
* `make_flatfield` has a hardcoded normalization section data[500:1500, 2000:2500] designed for full 2048×4612 LBC chips — returns empty on any smaller array. The 8 make_flatfield tests still pass (they check structure, not values); the 3 go_flatfield tests use a synthetic master flat to work around this.
* `go_extractchips` has the same mv-via-Popen cwd bug as make_targetdirectories had before we fixed it (line 940). Worth fixing the same way (shutil.move with absolute path).

* `LBCgo/whatidid.astrom.py` is an orphan script with a dot in its name sitting inside the package directory — not a valid Python module name and invisible to normal imports. Should probably be moved to docs/ or references/.

* `LBCgo/examples/` contains five scripts (pypit.py, pyhrs.py, modsProc.py, reduceccd_example.py, wht_basicredux.py) that are borrowed from other reduction pipelines and have broken relative imports. They are excluded from the docs but their presence in the package directory is a hygiene issue.

* `LBCgo/__init__.py` is empty — the package exposes no public __version__ or __all__, so LBCgo.__version__ (referenced in installation.rst) doesn't actually work yet. It should either be added there or the install verification example should use import LBCgo; print(LBCgo.lbcproc.__file__) instead.
