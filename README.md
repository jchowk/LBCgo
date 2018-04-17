# LBCgo: LBC data reduction pipeline

README GOES HERE...

## Dependencies:

* `astropy`
* `CCDProc`
* ``

The external C++ codes `Scamp`, `Swarp`, and `Sextractor` (provided by Emmanuel Bertin and collaborators; http://astromatic.iap.fr).

## Running LBCgo:

For "standard" situations, the `LBCgo` can be run in one step from the python command line. In this case, all of the data in the `raw/` directory are taken on the same night and have appropriate calibrations. In this case, running `LBCgo` from the command line is as simple as:
```
ipython> lbcgo()
```
Before doing this, copy the parameter files from `LBCgo/LBCgo/conf/` into the current working directory (an eventual fix won't require this step).

Alternatively, it can be useful to process each filter separately or even to avoid doing the astrometric steps until a later time. In this case, one may do:
```
ipython> lbcgo(filter_names=['I-BESSEL'], do_astrometry=False)
```

The astrometric portion of the reduction can be done later using, for example reducing the I-BESSEL data for the target PG1338+101:
```
ipython> fltr_dirs=glob('PG1338+101/I-BESSEL/')
ipython> go_register(fltr_dirs, do_sextractor=True,
            do_scamp=True, do_swarp=True)
```

## Some things that might go wrong:

Testing has revealed some occasional issues with the astrometric solution for the individual chips. This can be difficult to diagnose. The registration step using `SWARP` can warn you of some obvious cases, and these can subsequently be removed before rerunning the `SWARP` step by doing, e.g.:
```
ipython> go_register(fltr_dirs, do_sextractor=False,
            do_scamp=False, do_swarp=True)
```

There are several issues related to missing or inappropriate files that the current code does not deal with gracefully. The most common is missing flat fields or missing configuration files (found in `LBCgo/LBCgo/conf/`).


## Credit:

This pipeline is built on code initially developed by David Sands, and eventually incorporated into scripts made available by Ben Weiner
(https://github.com/bjweiner/LBC-reduction).

`LBCgo` was designed to simplify the process of LBC reduction, removing the need for IDL or IRAF in favor of Python. This package continues to require `SCAMP`, `SWARP`, and `SExtractor` provided by Emmanuel Bertin (http://astromatic.iap.fr). It makes extensive use of the `astropy`-affiliated package `CCDProc`.


## Known bugs / limitations:

* As of yet no tests are performed to separate LBCB / LBCR images taken with the V-BESSEL filter (which exists in both imagers). Care must be taken to avoid having both in the same directory.

* If flat field images are present, but no image is taken in that flat, an unfortunate behavior results (existing flat fields are divided by the unmatched flats).

* Flat field images taken as "test" images, including only a partial read-out of a single CCD, will cause the code to bail without a helpful error message.
