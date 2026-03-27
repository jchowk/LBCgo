Pipeline Overview
=================

LBCgo reduces LBC multi-extension FITS files through a sequence of stages,
controlled by the top-level function :func:`~LBCgo.lbcproc.lbcgo`. Each stage
can also be called independently via the lower-level functions in
:mod:`LBCgo.lbcproc` and :mod:`LBCgo.lbcregister`.

Stage 1: Overscan subtraction and trimming
------------------------------------------

:func:`~LBCgo.lbcproc.go_overscan`

A fourth-order polynomial is fit to the overscan strip (defined by the
``BIASSEC`` header keyword) and subtracted from each row. The science region
is then trimmed to ``TRIMSEC``. Output files carry the ``_over.fits`` suffix.

Several LBC header keywords that confuse downstream astrometric tools are
removed at this stage: ``BIASSEC``, ``TRIMSEC``, ``DATASEC``, ``ROTANGLE``,
``PARANGLE``, and all duplicate astrometric keywords ending in ``A``
(introduced by the LBC readout electronics).

Stage 2: Master bias (optional)
--------------------------------

:func:`~LBCgo.lbcproc.make_bias` / :func:`~LBCgo.lbcproc.go_bias`

Bias frames (``IMAGETYP=zero``) are median-combined with 3σ sigma clipping to
produce ``zero.fits``. Bias subtraction is disabled by default
(``bias_proc=False``) because the LBC overscan provides a reliable per-image
bias level. Enable with ``bias_proc=True`` if your data require it.

Stage 3: Master flat field
---------------------------

:func:`~LBCgo.lbcproc.make_flatfield` / :func:`~LBCgo.lbcproc.go_flatfield`

Sky flat frames (``IMAGETYP=flat``, ``OBJECT=SkyFlat``) are overscan-subtracted
and median-combined per filter. Frames with a mean pixel value above 55,000 ADU
are rejected as saturated before combination. The combined flat is normalised
and applied to science frames. Output flat-fielded files carry the ``_flat.fits``
suffix.

Stage 4: Target directory organisation
---------------------------------------

:func:`~LBCgo.lbcproc.make_targetdirectories`

Flat-fielded science files are moved into subdirectories organised by
``OBJECT`` header keyword and filter, e.g. ``NGC891/g-SLOAN/``. Object names
with spaces are compressed (``NGC 891`` → ``NGC891``).

Stage 5: Chip extraction
------------------------

:func:`~LBCgo.lbcproc.go_extractchips`

Each multi-extension FITS file is split into per-chip single-extension files
named ``<base>_1.fits`` through ``<base>_4.fits``. These are the inputs to the
astrometric pipeline. Chips can be subsetted via ``lbc_chips``.

Stage 6: Astrometric registration (optional)
---------------------------------------------

:func:`~LBCgo.lbcregister.go_register`

Controlled by ``do_astrometry=True`` (the default). Requires SExtractor, SCAMP,
and SWarp to be installed.

**SExtractor** (:func:`~LBCgo.lbcregister.go_sextractor`): Detects sources and
writes a FITS_LDAC catalog (``*.cat``) for each chip image.

**SCAMP** (:func:`~LBCgo.lbcregister.go_scamp`): Matches the catalog against
an astrometric reference (default: GAIA-DR3) and writes a ``.head`` file
containing the WCS solution. SCAMP is run iteratively (default 3 iterations)
with progressively tighter tolerance parameters to improve the astrometric fit.

**SWarp** (:func:`~LBCgo.lbcregister.go_swarp`): Resamples each chip image
onto a common grid and co-adds them into a final mosaic. Output files are named
``<object>.<filter>.mos.fits`` with a companion weight map.

Known limitations
-----------------

- **V-BESSEL shared between cameras**: The V-BESSEL filter is available in both
  LBC-Blue and LBC-Red. If both cameras observed in V-BESSEL on the same night,
  run the pipeline separately with ``lbcr=False`` and ``lbcb=False`` to prevent
  inappropriate co-addition.

- **Unmatched flat fields**: If a flat field is present for a filter but no science
  frame was taken in that filter (or vice versa), the behaviour can be unexpected.
  Inspect your calibration file inventory before running.

- **Partial-readout test images**: Flat or bias frames taken as partial-readout
  test exposures will cause the pipeline to fail without a helpful error message.
  Remove these from ``raw/`` before running.

- **Astrometric fit quality**: No automated quality assessment of the astrometric
  solution is currently implemented. Inspect the SCAMP ``.xml`` diagnostic files
  or the ``astro*.pdf`` plots to identify chips with poor solutions before running
  SWarp.

- **Extended objects**: Background estimation and astrometric fitting can be
  unreliable for fields dominated by large extended sources (e.g., nearby galaxies).

Credits
-------

LBCgo is built on code originally developed by David Sands and subsequently
extended in scripts made available by Ben Weiner
(`bjweiner/LBC-reduction <https://github.com/bjweiner/LBC-reduction>`_).
It makes extensive use of the astropy-affiliated package
`ccdproc <https://ccdproc.readthedocs.io>`_ and the
`Astromatic <http://astromatic.iap.fr>`_ tools by Emmanuel Bertin.
