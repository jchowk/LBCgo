!cp /Users/howk/python/LBCgo/LBCgo/conf/* ./
%load_ext autoreload
%autoreload 2

#After adding the LBCgo directories to PYTHONPATH:
from lbcproc import *
from lbcregister import *

from astropy.io import fits
from ccdproc import CCDData, ImageFileCollection


lbc_chips=[1,2,3,4]
fltr_dirs=glob('I-BESSEL/')
go_register(fltr_dirs, lbc_chips=lbc_chips,
            do_sextractor=True, do_scamp=False, do_swarp=False)

go_register(fltr_dirs, lbc_chips=lbc_chips,scamp_iterations=3,
            do_sextractor=False, do_scamp=True, do_swarp=False)




## Found an image with *2* high S/N objects in the SCAMP solution.
!cp I-BESSEL/lbcr.20100318.091047_3.fits Test/

lbc_chips=[1,2,3,4]
fltr_dirs=glob('Test/')
go_register(fltr_dirs, lbc_chips=lbc_chips,
            do_sextractor=True, do_scamp=False, do_swarp=False)

go_register(fltr_dirs, lbc_chips=lbc_chips,scamp_iterations=4,
            astroref_catalog='GSC-2.3',
            do_sextractor=False, do_scamp=True, do_swarp=False)

go_register(fltr_dirs, lbc_chips=lbc_chips,scamp_iterations=5,
            do_sextractor=False, do_scamp=True, do_swarp=False)


go_register(fltr_dirs, lbc_chips=lbc_chips,scamp_iterations=3,
            do_sextractor=False, do_scamp=False, do_swarp=True)






#Try to do this with multiple exposures using the INSTRUMENT stability type.

scamp \
I-BESSEL/lbcr.20100318.091047_1.cat \
I-BESSEL/lbcr.20100318.092032_1.cat \
I-BESSEL/lbcr.20100318.093008_1.cat \
I-BESSEL/lbcr.20100318.095531_1.cat \
I-BESSEL/lbcr.20100319.094300_1.cat \
I-BESSEL/lbcr.20100319.095307_1.cat \
I-BESSEL/lbcr.20100318.091047_2.cat \
I-BESSEL/lbcr.20100318.092032_2.cat \
I-BESSEL/lbcr.20100318.093008_2.cat \
I-BESSEL/lbcr.20100318.095531_2.cat \
I-BESSEL/lbcr.20100319.094300_2.cat \
I-BESSEL/lbcr.20100319.095307_2.cat \
I-BESSEL/lbcr.20100318.091047_3.cat \
I-BESSEL/lbcr.20100318.092032_3.cat \
I-BESSEL/lbcr.20100318.093008_3.cat \
I-BESSEL/lbcr.20100318.095531_3.cat \
I-BESSEL/lbcr.20100319.094300_3.cat \
I-BESSEL/lbcr.20100319.095307_3.cat \
I-BESSEL/lbcr.20100318.091047_4.cat \
I-BESSEL/lbcr.20100318.092032_4.cat \
I-BESSEL/lbcr.20100318.093008_4.cat \
I-BESSEL/lbcr.20100318.095531_4.cat \
I-BESSEL/lbcr.20100319.094300_4.cat \
I-BESSEL/lbcr.20100319.095307_4.cat \
-c scamp.lbc.conf  -STABILITY_TYPE INSTRUMENT \
-PIXSCALE_MAXERR 1.2 -POSANGLE_MAXERR 3.0 -POSITION_MAXERR 1.0 \
-DISTORT_DEGREES 3 \
-MOSAIC_TYPE LOOSE \
-ASTREF_CATALOG GAIA-DR1 -AHEADER_SUFFIX .ahead \
-CROSSID_RADIUS 10.0

scamp \
I-BESSEL/lbcr.20100318.091047_1.cat \
I-BESSEL/lbcr.20100318.092032_1.cat \
I-BESSEL/lbcr.20100318.093008_1.cat \
I-BESSEL/lbcr.20100318.095531_1.cat \
I-BESSEL/lbcr.20100319.094300_1.cat \
I-BESSEL/lbcr.20100319.095307_1.cat \
I-BESSEL/lbcr.20100318.091047_2.cat \
I-BESSEL/lbcr.20100318.092032_2.cat \
I-BESSEL/lbcr.20100318.093008_2.cat \
I-BESSEL/lbcr.20100318.095531_2.cat \
I-BESSEL/lbcr.20100319.094300_2.cat \
I-BESSEL/lbcr.20100319.095307_2.cat \
I-BESSEL/lbcr.20100318.091047_3.cat \
I-BESSEL/lbcr.20100318.092032_3.cat \
I-BESSEL/lbcr.20100318.093008_3.cat \
I-BESSEL/lbcr.20100318.095531_3.cat \
I-BESSEL/lbcr.20100319.094300_3.cat \
I-BESSEL/lbcr.20100319.095307_3.cat \
I-BESSEL/lbcr.20100318.091047_4.cat \
I-BESSEL/lbcr.20100318.092032_4.cat \
I-BESSEL/lbcr.20100318.093008_4.cat \
I-BESSEL/lbcr.20100318.095531_4.cat \
I-BESSEL/lbcr.20100319.094300_4.cat \
I-BESSEL/lbcr.20100319.095307_4.cat \
-c scamp.lbc.conf  -STABILITY_TYPE INSTRUMENT \
-PIXSCALE_MAXERR 1.1 -POSANGLE_MAXERR 1.0 -POSITION_MAXERR 0.5 \
-DISTORT_DEGREES 3 \
-MOSAIC_TYPE FIX_FOCALPLANE \
-ASTREF_CATALOG GAIA-DR1 -AHEADER_SUFFIX .head \
-CROSSID_RADIUS 5.0


######################################################################

scamp \
I-BESSEL/lbcr.20100318.091047_2.cat \
I-BESSEL/lbcr.20100318.092032_2.cat \
I-BESSEL/lbcr.20100318.093008_2.cat \
I-BESSEL/lbcr.20100318.095531_2.cat \
I-BESSEL/lbcr.20100319.094300_2.cat \
I-BESSEL/lbcr.20100319.095307_2.cat \
-c scamp.lbc.conf  -STABILITY_TYPE INSTRUMENT \
-PIXSCALE_MAXERR 1.1 -POSANGLE_MAXERR 1.0 -POSITION_MAXERR 1.0 \
-DISTORT_DEGREES 3 \
-MOSAIC_TYPE LOOSE \
-ASTREF_CATALOG GAIA-DR1 -AHEADER_SUFFIX .ahead \
-CROSSID_RADIUS 10.0 

scamp \
I-BESSEL/lbcr.20100318.091047_2.cat \
I-BESSEL/lbcr.20100318.092032_2.cat \
I-BESSEL/lbcr.20100318.093008_2.cat \
I-BESSEL/lbcr.20100318.095531_2.cat \
I-BESSEL/lbcr.20100319.094300_2.cat \
I-BESSEL/lbcr.20100319.095307_2.cat \
-c scamp.lbc.conf  -STABILITY_TYPE INSTRUMENT \
-PIXSCALE_MAXERR 1.1 -POSANGLE_MAXERR 1.0 -POSITION_MAXERR 0.5 \
-DISTORT_DEGREES 3 \
-MOSAIC_TYPE FIX_FOCALPLANE \
-ASTREF_CATALOG GAIA-DR1 -AHEADER_SUFFIX .head \
-CROSSID_RADIUS 5.0 

scamp \
I-BESSEL/lbcr.20100318.091047_2.cat \
I-BESSEL/lbcr.20100318.092032_2.cat \
I-BESSEL/lbcr.20100318.093008_2.cat \
I-BESSEL/lbcr.20100318.095531_2.cat \
I-BESSEL/lbcr.20100319.094300_2.cat \
I-BESSEL/lbcr.20100319.095307_2.cat \
-c scamp.lbc.conf  -STABILITY_TYPE INSTRUMENT \
-PIXSCALE_MAXERR 1.025 -POSANGLE_MAXERR 1.0 -POSITION_MAXERR 0.5 \
-DISTORT_DEGREES 3 \
-MOSAIC_TYPE FIX_FOCALPLANE \
-ASTREF_CATALOG GAIA-DR1 -AHEADER_SUFFIX .head \
-CROSSID_RADIUS 1.0 -MATCH N

##### 3
scamp \
I-BESSEL/lbcr.20100318.091047_3.cat \
I-BESSEL/lbcr.20100318.092032_3.cat \
I-BESSEL/lbcr.20100318.093008_3.cat \
I-BESSEL/lbcr.20100318.095531_3.cat \
I-BESSEL/lbcr.20100319.094300_3.cat \
I-BESSEL/lbcr.20100319.095307_3.cat \
-c scamp.lbc.conf  -STABILITY_TYPE INSTRUMENT \
-PIXSCALE_MAXERR 1.1 -POSANGLE_MAXERR 1.0 -POSITION_MAXERR 1.0 \
-DISTORT_DEGREES 3 \
-MOSAIC_TYPE LOOSE \
-ASTREF_CATALOG GAIA-DR1 -AHEADER_SUFFIX .ahead \
-CROSSID_RADIUS 10.0 

scamp \
I-BESSEL/lbcr.20100318.091047_3.cat \
I-BESSEL/lbcr.20100318.092032_3.cat \
I-BESSEL/lbcr.20100318.093008_3.cat \
I-BESSEL/lbcr.20100318.095531_3.cat \
I-BESSEL/lbcr.20100319.094300_3.cat \
I-BESSEL/lbcr.20100319.095307_3.cat \
-c scamp.lbc.conf  -STABILITY_TYPE INSTRUMENT \
-PIXSCALE_MAXERR 1.1 -POSANGLE_MAXERR 1.0 -POSITION_MAXERR 0.5 \
-DISTORT_DEGREES 3 \
-MOSAIC_TYPE FIX_FOCALPLANE \
-ASTREF_CATALOG GAIA-DR1 -AHEADER_SUFFIX .head \
-CROSSID_RADIUS 5.0 -MATCH N

scamp \
I-BESSEL/lbcr.20100318.091047_3.cat \
I-BESSEL/lbcr.20100318.092032_3.cat \
I-BESSEL/lbcr.20100318.093008_3.cat \
I-BESSEL/lbcr.20100318.095531_3.cat \
I-BESSEL/lbcr.20100319.094300_3.cat \
I-BESSEL/lbcr.20100319.095307_3.cat \
-c scamp.lbc.conf  -STABILITY_TYPE INSTRUMENT \
-PIXSCALE_MAXERR 1.025 -POSANGLE_MAXERR 1.0 -POSITION_MAXERR 0.5 \
-DISTORT_DEGREES 3 \
-MOSAIC_TYPE FIX_FOCALPLANE \
-ASTREF_CATALOG GAIA-DR1 -AHEADER_SUFFIX .head \
-CROSSID_RADIUS 1.0 -MATCH N

