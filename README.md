# LBCgo: LBC data reduction pipeline

README GOES HERE...

## Dependencies:

* `astropy`
* `CCDProc`
* ``

The external C++ codes `Scamp`, `Swarp`, and `Sextractor` (provided by Emmanuel Bertin and collaborators; http://astromatic.iap.fr).

## Credit:

This pipeline is built on code initially developed by David Sands, and eventually incorporated into scripts made available by Ben Weiner
(http://mingus.as.arizona.edu/~bjw/software/). Subsequent development by Neil Crighton (https://github.com/nhmc/LBC_redux) started the push toward Python; several of the routines here are built off of that work.

`LBCgo` was designed to simplify the process of LBC reduction, removing the need for IDL or IRAF in favor of Python. This package continues to require Scamp, Swarp and Sextractor provided by Emmanuel Bertin (http://astromatic.iap.fr).
