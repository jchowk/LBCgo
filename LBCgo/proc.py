
# Holder: information useful for CR rejection
# lacos_im = OrderedDict()
# lacos_im['gain']=1.75
# lacos_im['readn']=12.0
# lacos_im['sigclip']=4.5
# lacos_im['sigfrac']=0.3
# lacos_im['objlim']=1.0
# lacos_im['niter']=2

# Import the utils code here. It has to be done outside of a function in py3.


def proclbc(raw_directory='./raw/', image_directory='./',
            filter_names=None, bias_proc=False,
            verbose=True, clean=True):
    """Process a directory of LBC data.

    By default, all raw data are to be in the ./raw/ directory before processing (raw_directory='./raw/').
    By default, all new images are written in the CWD (image_directory='./').




    """
    import numpy as np
    from glob import glob

    from .utils import make_bias, make_bpm, make_flatfield
    from .utils import do_overscan, do_bias, do_flatfield

    import ccdproc
    from ccdproc import ImageFileCollection, CCDData

    from astropy.io import fits
    import astropy.units as u
    import astropy.constants as c


    # It will go something like this...

    # Make sure the input directories have trailing slashes:
    if raw_directory[-1] != '/':
        raw_directory += '/'
    if image_directory[-1] != '/':
        image_directory += '/'

    ########## Collect basic image information
    # Find the raw files we'll use
    lbc_file_base = 'lbc?.????????.??????.fits'
    raw_lbc_files = glob(raw_directory+lbc_file_base)

    # What information do we want from them?
    keywds = ['object', 'filter', 'exptime', 'imagetyp', 'propid', 'lbcobnam',
              'airmass', 'HA', 'objra', 'objdec']

    ##### Create an ImageFileCollection object to hold the raw data list.
    if np.int(ccdproc.__version__[0]) == 2:
        ic0 = ImageFileCollection(image_directory, keywords=keywds,
                                  glob_include=lbc_file_base)
    else:
        ic0 = ImageFileCollection(image_directory, keywords=keywds,
                                  filenames=raw_lbc_files)

    ######### Create the master bias frame (if requested)
    # if bias_proc == True:
    # TODO implement bias creation, testing
    if bias_proc == True:
        make_bias(ic0,
                  image_directory=image_directory,raw_directory=raw_directory)

    ###### Per filter:
    # Do this on a per filter basis. Right now I'm just laying
    # out where to go.

    # The filters to go through are all those in the IC file unless otherwise specified.
    if filter_names == None:
        filter_names = np.unique(ic0.summary['filter'])

    # Loop through each of the filters
    for filter in filter_names:
        # Find the images with the filter of interest.
        if verbose == True:
            print('Processing {0} files.'.format(filter))

        # List of images in the current filter
        ic1 = ImageFileCollection(raw_directory,
                        keywords=keywds,
                        filenames=(ic0.files_filtered(filter=filter)).tolist())

        # Make master flat fields. Could be done for all filters at once, but keeping it here for now.
        make_flatfield(ic1,verbose=verbose,
                       raw_directory=raw_directory,image_directory=image_directory)

        # Remove overscan, trim object files.
        overfiles = do_overscan(ic1,verbose=verbose,
                                raw_directory=raw_directory,image_directory=image_directory)

        # Apply bias frames
        if bias_proc == True:
            # zero_files = do_bias(verbose=verbose)
            print('')

        # Image collection of object frames overscanned & bias subtracted
        ic2 = ImageFileCollection(image_directory,
                    keywords=keywds,filenames=overfiles)

        # Apply the flatfields
        flatfiles = do_flatfield(ic2, image_directory=image_directory, verbose=verbose)

        # Image collection of object frames overscanned & bias subtracted + flattened
        ic3 = ImageFileCollection(image_directory,
                    keywords=keywds,filenames=flatfiles)

        # Create directories for extracting individual chips.
        tgt_dirs, fltr_dirs = mk_targetdirectories(ic3,
                                        image_directory = image_directory, verbose = verbose)

        # TODO: Extract individual chips to object / filter directories
        # Extract individual chips
        extract_chips(ic3, image_directory = image_directory, verbose = verbose)

        # TODO: do_sextractor
        # TODO: do_scampswarp



    # Let's do some clean-up.
    if clean:
        # We will eventually ...
        #remove _over, _zero, _flat files.
        print('')
    else:
        # We will eventually ...
        # Move _over, _zero, _flat files.
        print('')