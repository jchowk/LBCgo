# LBCgo: LBC data reduction pipeline

README GOES HERE...

## Dependencies:

* `astropy`
* `CCDProc`
* ``

The external C++ codes `Scamp`, `Swarp`, and `Sextractor` (provided by Emmanuel Bertin and collaborators; http://astromatic.iap.fr).

## Credit:

This pipeline is built on code initially developed by David Sands, and eventually incorporated into scripts made available by Ben Weiner
(https://github.com/bjweiner/LBC-reduction).

`LBCgo` was designed to simplify the process of LBC reduction, removing the need for IDL or IRAF in favor of Python. This package continues to require Scamp, Swarp and Sextractor provided by Emmanuel Bertin (http://astromatic.iap.fr). It makes extensive use of the `astropy`-affiliated package `CCDProc`.
