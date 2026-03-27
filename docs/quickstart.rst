Quickstart
==========

Directory layout
----------------

LBCgo expects raw LBC FITS files in a ``raw/`` subdirectory of your working
directory. Files should follow the LBT naming convention
``lbc[b|r].YYYYMMDD.NNNNNN.fits``. Calibration frames (biases, sky flats) and
science frames can all live together in ``raw/``; the pipeline sorts them by
``IMAGETYP`` header keyword.

Expected layout before running::

    my_night/
    └── raw/
        ├── lbcb.20230101.000001.fits   # LBC-Blue science frame
        ├── lbcb.20230101.000002.fits
        ├── lbcr.20230101.000010.fits   # LBC-Red science frame
        └── ...

Standard run
------------

For a typical night where all calibration data are present, reduce everything
in one call::

    from LBCgo.lbcproc import lbcgo
    lbcgo()

This processes both LBC-Blue and LBC-Red data for all filters found in
``./raw/``, applies overscan correction, master flat fields, extracts
individual chips, and performs astrometric registration against GAIA-DR3.

Output files are organised into subdirectories by target name and filter,
e.g. ``NGC891/g-SLOAN/``.

Common options
--------------

**Process a single filter only**::

    lbcgo(filter_names=['g-SLOAN'])

**Skip astrometry** (useful for a quick first look, or when the external
tools are not installed)::

    lbcgo(do_astrometry=False)

**LBC-Blue data only**::

    lbcgo(lbcb=True, lbcr=False)

**Missing or disabled chips** (e.g., chip 3 was offline)::

    lbcgo(lbc_chips=[1, 2, 4])

**Apply bias frames** (not done by default because LBC bias levels are
stable and well-characterised by the overscan)::

    lbcgo(bias_proc=True)

Multi-filter example (both cameras)
------------------------------------

To reduce SDSS r-band from LBC-Red and g-band from LBC-Blue separately::

    from LBCgo.lbcproc import lbcgo

    # Red camera: r-band
    lbcgo(lbcr=True, lbcb=False, filter_names=['r-SLOAN'])

    # Blue camera: g and U
    lbcgo(lbcr=False, lbcb=True, filter_names=['g-SLOAN', 'SDT_Uspec'])

Running astrometry separately
------------------------------

The astrometric step can be deferred or re-run independently using
:func:`~LBCgo.lbcregister.go_register`::

    from glob import glob
    from LBCgo.lbcregister import go_register

    # Run full astrometric pipeline on previously reduced data
    fltr_dirs = glob('NGC891/*/')
    go_register(fltr_dirs)

Run the steps individually if you need to inspect intermediate results or
re-run just SCAMP and SWarp::

    # Step 1: Source extraction only
    go_register(fltr_dirs, do_sextractor=True, do_scamp=False, do_swarp=False)

    # Step 2: Astrometric calibration (more iterations for precision)
    go_register(fltr_dirs, do_sextractor=False, do_scamp=True, do_swarp=False,
                scamp_iterations=5)

    # Step 3: Image co-addition (after removing any frames with bad solutions)
    go_register(fltr_dirs, do_sextractor=False, do_scamp=False, do_swarp=True)

Output products
---------------

After a full reduction run, the working directory contains:

- ``<target>/<filter>/`` — flat-fielded individual chip images (``*_1.fits`` … ``*_4.fits``)
- ``<target>/<filter>/<object>.<filter>.mos.fits`` — astrometrically registered co-add (if ``do_astrometry=True``)
- ``<target>/<filter>/<object>.<filter>.mos.weight.fits`` — corresponding weight map
- ``data/`` — intermediate files (overscan-corrected, flat-fielded MEF images)
