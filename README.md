# LBCgo: LBC data reduction pipeline

WARNING: This code is currently under continued development. While the basic functionality exists, it should be used with some care and attention.

## Dependencies:

Python dependencies:
* `astropy`
* `CCDProc`
* `numpy`


## External Dependencies:
* `SExtractor` - Source extraction software
* `SCAMP` - Astrometric calibration software  
* `SWarp` - Image resampling and co-addition software

These can be installed via:
- **macOS**: `brew install sextractor scamp swarp`
- **Ubuntu/Debian**: `apt-get install sextractor scamp swarp`
- **From source**: Available at http://astromatic.iap.fr

These can also be installed using `conda-forge`:
```
conda install astromatic-scamp
conda install astromatic-swarp
conda install sextractor
```

## Installing:

`LBCgo` is available through `pip`. The easiest path to installation is:

```
pip install lbcgo
```

## Running LBCgo:

For "standard" situations, the `LBCgo` can be run in one step from the python command line. In this case, all of the data in the `raw/` directory are taken on the same night and have appropriate calibrations. In this case, running `LBCgo` from the command line is as simple as:
```python
# Modern Python import style:
from LBCgo import lbcgo

# Run the main reduction pipeline
lbcgo()
```

Before doing this, copy the parameter files from `LBCgo/LBCgo/conf/` into the current working directory (an eventual fix won't require this step).

Alternatively, it can be useful to process each filter separately or even to avoid doing the astrometric steps until a later time. In this case, one may do:
```
# Run the main reduction pipeline with just the I-band filter:
lbcgo(filter_names=['I-BESSEL'], do_astrometry=False)
```

The astrometric portion of the reduction can be done later using, for example reducing the I-BESSEL data for the target PG1338+101:
```
from lbcgo.lbcregister import go_register
fltr_dirs=glob('PG1338+101/I-BESSEL/')
go_register(fltr_dirs, do_sextractor=True,
            do_scamp=True, do_SWarp=True)
```

#### Missing chips:

`LBCgo` can be used if the images were taken when one of the LBC CCDs was off-line. The approach to doing this is to explicitly specify the chips to include in the data reduction steps:
```
lbcgo(lbc_chips=[1,2,4])
```
This is useful, as there were several months in 2011 when LBCCHIP3 was inoperable.

## Some things that might go wrong:

Testing has revealed some occasional issues with the astrometric solution for the individual chips. This can be difficult to diagnose. The registration step using `SWarp` can warn you of some obvious cases, and these can subsequently be removed before rerunning the `SWarp` step by doing, e.g.:
```
go_register(fltr_dirs, do_sextractor=False,
            do_scamp=False, do_SWarp=True)
```

There are several issues related to missing or inappropriate files that the current code does not deal with gracefully. The most common is missing flat fields or missing configuration files (found in `LBCgo/LBCgo/conf/`).


## Credit:

This pipeline is built on code initially developed by David Sands, and eventually incorporated into scripts made available by Ben Weiner
(https://github.com/bjweiner/LBC-reduction).

`LBCgo` was designed to simplify the process of LBC reduction, removing the need for IDL or IRAF in favor of Python. This package continues to require `SCAMP`, `SWarp`, and `SExtractor` provided by Emmanuel Bertin (http://astromatic.iap.fr). It makes extensive use of the `astropy`-affiliated package `CCDProc`.


## Known bugs / limitations:

* As of yet no tests are performed to separate LBCB / LBCR images taken with the V-BESSEL filter (which exists in both imagers). Care must be taken to avoid having both in the same directory.

* If flat field images are present, but no image is taken in that flat, an unfortunate behavior results (existing flat fields are divided by the unmatched flats).

* Flat field images taken as "test" images, including only a partial read-out of a single CCD, will cause the code to bail without a helpful error message.
