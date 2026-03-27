LBCgo
=====

LBCgo is a Python pipeline for reducing data from the Large Binocular Telescope's
Large Binocular Camera (LBC). It processes multi-extension FITS files from LBC-Blue
and LBC-Red through overscan removal, bias subtraction, flat fielding, chip
extraction, and optional astrometric registration via the Astromatic suite
(SExtractor, SCAMP, SWarp).

.. warning::

   LBCgo is under active development. Use with care and attention to the
   :doc:`known limitations <pipeline>`.

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Reference

   pipeline
   autoapi/index
