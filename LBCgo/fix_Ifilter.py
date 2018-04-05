import numpy as np
from glob import glob

from astropy.io import fits

from ccdproc import  ImageFileCollection

# Suppress some of the WCS warnings
import warnings
from astropy.utils.exceptions import AstropyWarning,AstropyUserWarning
warnings.filterwarnings('ignore', category=AstropyWarning, append=True)
warnings.filterwarnings('ignore', category=AstropyUserWarning, append=True)


def fix_Ifilter(raw_directory):
        keywds = ['object', 'filter', 'exptime', 'imagetyp',
                'propid', 'lbcobnam',
                'airmass', 'HA', 'objra', 'objdec']

        lbc_file_base = 'lbcr.*.*.fits*'
        ##### Create an ImageFileCollection object to hold the raw data list.
        if np.int(ccdproc.__version__[0]) == 2:
            ic0 = ImageFileCollection(raw_directory, keywords=keywds,
                                      glob_include=lbc_file_base)
            num_images = np.size(ic0.summary['object'])
            # Exit if (for some reason) no images are found in the raw_directory.
            if num_images == 0:
                print('WARNING: No images found.')
                return None
        else:
            raw_lbc_files = glob(raw_directory+lbc_file_base)
            ic0 = ImageFileCollection(raw_directory, keywords=keywds,
                                      filenames=raw_lbc_files)
            # Exit if (for some reason) no images are found in the raw_directory.
            if num_images == 0:
                print('WARNING: No images found.')
                return None

        filter_names = ic0.values('filter',unique=True)

        if 'I-BESSEL(?)' in filter_names:
            # Pull out the "I-BESSEL(?)" images:
            icX = ImageFileCollection(raw_directory,
                            keywords=keywds,
                            filenames=(ic0.files_filtered(filter='I-BESSEL(?)')).tolist())

            for fl in icX.files:
                fits.setval(fl,'FILTER',value='I-BESSEL')
                print("Changed FILTER keyword: 'I-BESSEL(?)' --> 'I-BESSEL' in {0}".format(fl))
        else:
            print("No files contain FILTER = 'I-BESSEL(?)'.")
