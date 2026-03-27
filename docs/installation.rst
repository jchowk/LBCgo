Installation
============

Python package
--------------

LBCgo is available on PyPI::

    pip install lbcgo

For development installation from source::

    git clone https://github.com/jchowk/LBCgo.git
    cd LBCgo
    pip install -e .

Conda environment
-----------------

A pre-defined conda environment file is included in the repository::

    conda env create -f lbcgo_environment.yaml
    conda activate lbcgo

Python dependencies (``numpy``, ``astropy``, ``ccdproc``) are installed
automatically via either method above.

External tools (Astromatic suite)
----------------------------------

Astrometric processing requires three external tools from the
`Astromatic <http://astromatic.iap.fr>`_ software suite:

- **SExtractor** — source detection and cataloguing
- **SCAMP** — astrometric calibration against reference catalogs
- **SWarp** — image resampling and co-addition

These are only needed if ``do_astrometry=True`` (the default) is used in
:func:`~LBCgo.lbcproc.lbcgo`. If you are only doing basic reduction, they
can be omitted.

**macOS (Homebrew)**::

    brew install sextractor scamp swarp

**Ubuntu/Debian**::

    apt-get install sextractor scamp swarp

**conda-forge**::

    conda install -c conda-forge astromatic-scamp astromatic-swarp sextractor

Verifying the installation
--------------------------

After installation, check that Python dependencies are present::

    python -c "import LBCgo; print(LBCgo.__version__)"

If using astrometric processing, verify the external tools are on your PATH::

    which sex
    which scamp
    which swarp
