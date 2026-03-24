# Copy configuration files
# ! cp /Users/howk/python/LBCgo/LBCgo/conf/* ./

# Import LBCgo procedures
from lbcproc import *
from lbcregister import *

# Red side data: SDSS r-band:
fltr=['r-SLOAN']
fltr_dirs = glob('*/r-SLOAN/')
lbcgo(lbcb=False,filter_names=fltr)

# Red side data: SDSS i-band excluding Chip 4:
fltr = ['i-SLOAN']
lbcgo(lbcb=False,filter_names=fltr,clean=False,chips=[1,2,3])

# Blue side data: U, g
fltr = ['SDT_Uspec','g-SLOAN']
lbcgo(filter_names=fltr,lbcr=False)

# Testing the alignment steps one-by-one
fltr_dirs = glob('*/i-SLOAN/')
# First Source Extractor
go_register(fltr_dirs,
            do_sextractor=True, do_scamp=False, do_swarp=False)
# Next SCAMP
go_register(fltr_dirs,
            do_sextractor=False, do_scamp=True, do_swarp=False,
            scamp_iterations=5)
# Last SWARP
go_register(fltr_dirs,
            do_sextractor=False, do_scamp=False, do_swarp=True)
