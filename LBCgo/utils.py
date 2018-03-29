import numpy as np
import os
import shlex
from subprocess import call, Popen
from glob import glob

import pdb

import astropy.units as u
from astropy.modeling import models

from astropy.io import fits

# Suppress some of the WCS warnings
import warnings
from astropy.utils.exceptions import AstropyWarning,AstropyUserWarning
warnings.filterwarnings('ignore', category=AstropyWarning, append=True)
warnings.filterwarnings('ignore', category=AstropyUserWarning, append=True)

import ccdproc
from ccdproc import  ImageFileCollection,CCDData

# Holder: information useful for CR rejection
# lacos_im = OrderedDict()
# lacos_im['gain']=1.75
# lacos_im['readn']=12.0
# lacos_im['sigclip']=4.5
# lacos_im['sigfrac']=0.3
# lacos_im['objlim']=1.0
# lacos_im['niter']=2

# TODO FIXPIX!
# TODO Cosmic ray cleaning

# TODO Create option for cleaning intermediate steps (_over, _zero, _flat)
# TODO Do gain correction, uncertainty system?
# TODO More general MEF flat fielding?
# TODO Logging system

# TODO Saturation correction?
# TODO image weights
# TODO clean up


def do_overscan(image_collection, objects_only=True,
                image_directory='./',raw_directory='./raw/',
                verbose=True, return_files = True):
    """Remove overscan and trim images."""
    ##
    ## Fit, subtract overscan
    ##
    ## Notes:
    ##   image_collection -- A ccdproc-style image collection.
    ##   image_directory  -- Where to store the output images.
    ##   raw_directory    -- Where to find the raw data
    ##   return_files     -- Return a list of subtracted files (default: True)

    # Hardwire the num. of LBC chips
    lbc_chips = [1,2,3,4]
    num_lbc_chips = np.size(lbc_chips)

    # By default we're only doing this to the object files.
    # Flats and biases are corrected when producing the master calibration images.
    if objects_only == True:
        over_files = image_collection.files_filtered(imagetyp='object')
    else:
        over_files = image_collection.files
    over_files_out=[]

    # Loop through the files
    for filename in over_files:
        # Create the output filename.
        output_filename = filename.split('.fits')[0]+'_over.fits'

        # Set up the output HDU list
        # Capture the 0th header
        base_header = fits.getheader(raw_directory+filename)
        master_hdu = fits.PrimaryHDU(header=base_header)
        # Start output HDU list:
        output_hdu = fits.HDUList([master_hdu])

        # Loop through the chips
        for chip in lbc_chips:
            # Create the CCDData version of this chip
            ccd = CCDData.read(raw_directory+filename, chip,
                               unit=u.adu)

            # Fit, subtract overscan
            poly_model = models.Polynomial1D(4)
            ccd = ccdproc.subtract_overscan(ccd,
                                        overscan_axis=1,
                                        model = poly_model,
                                        fits_section=ccd.header['BIASSEC'])
            # Trim the image
            ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])

            # Convert to 32 bit
            ccd.data=ccd.data.astype('float32')

            # Remove unneeded header keywords. Makes this consistent
            #   with IRAF treatment.
            temp_header = ccd.header
            del temp_header['trimsec']
            del temp_header['biassec']
            del temp_header['datasec']

            # Append the current chip into the hdu:
            output_hdu.append(ccd.to_hdu()[0])

        # Write the data
        output_hdu.writeto(output_filename, overwrite=True)

        # Keep track of what files we've written.
        over_files_out.append(output_filename)

        # Report:
        if verbose:
            print("Created {0}".format(output_filename))

    # Return the corrected filenames if requested (True is default).
    if return_files == True:
        return over_files_out

def make_bias(image_collection, bias_filename=None,
              image_directory='./', raw_directory='./raw/',
              verbose=True):
    """Make a master bias image for a collection of images."""
    ##
    ## Create master bias image
    ##
    ## Notes:
    ##   image_collection -- A ccdproc-style image collection.
    ##   image_directory  -- Where to store the output images.
    ##   raw_directory    -- Where to find the raw data


    # Report:
    if verbose:
        print('----- MAKE_BIAS: making master zero file. -----')


    # Pull out the bias images from our collection
    zero_files = image_collection.files_filtered(imagetyp='zero')
    # Create the output file name
    if bias_filename == None:
        zero_output_name = 'zero.fits'
    else:
        zero_output_name = bias_filename

    # Hardwire the number of CCDs in LBC
    lbc_chips = [1,2,3,4]
    num_lbc_chips = np.size(lbc_chips)

    # Create 2D list
    # (from https://stackoverflow.com/questions/7745562/appending-to-2d-lists-in-python)
    zero_list = [[] for i in range(num_lbc_chips)]
    num_zero_images = 0

    # Loop through the files
    for filename in zero_files:
        # Update zero counter
        num_zero_images += 1

        # This will serve to hold all of the master flats before writing
        if filename == zero_files[0]:
            # master_hdu = fits.open(raw_directory+filename)

            # Set up the output HDU list
            # Capture the 0th header
            base_header = fits.getheader(raw_directory + filename)
            master_hdu = fits.PrimaryHDU(header=base_header)
            # Start output HDU list:
            output_hdu = fits.HDUList([master_hdu])

        # Loop through the chips
        for chip in lbc_chips:
            # Create the CCDData version of this chip
            ccd = CCDData.read(raw_directory+filename, chip, unit=u.adu)

            # Fit, subtract overscan
            poly_model = models.Polynomial1D(4)
            ccd = ccdproc.subtract_overscan(ccd,
                                        overscan_axis=1,
                                        model = poly_model,
                                        fits_section=ccd.header['BIASSEC'])

            # Trim the image
            ccd = ccdproc.trim_image(ccd,
                                     fits_section=ccd.header['TRIMSEC'])

            # Convert to 32 bit
            ccd.data=ccd.data.astype('float32')

            # Add the current image to the list
            zero_list[chip-1].append(ccd)

    # After looping files, chips:
    # Create master bias from subtracted images
    for chip in lbc_chips:
        master_zero = ccdproc.combine(zero_list[chip-1],
                                      #dtype=np.float32,
                                      method='median',sigma_clip=True,
                                      sigma_clip_high_thresh=3.,
                                      sigma_clip_low_thresh=3.)
        # Convert to 32 bit
        master_zero.data=master_zero.data.astype('float32')
        # Add to output
        output_hdu.append(master_zero.to_hdu()[0])

        # Remove unneeded header keywords. Makes this consistent
        #   with IRAF treatment.
        # import IPython; IPython.embed()
        temp_header = output_hdu[chip].header
        del temp_header['trimsec']
        del temp_header['biassec']
        del temp_header['datasec']

        # Note how many files are combined
        temp_header['ncombine'] = num_zero_images

        # Fill chip-level master header
        output_hdu[chip].header = temp_header


    # Write the master bias
    output_hdu.writeto(image_directory+zero_output_name, overwrite=True)
    output_hdu.close()

    # Report:
    if verbose:
        print("Created {0} master bias frame.".format(zero_output_name))

def do_bias(image_collection, bias_file=None,
                 image_directory='./',input_directory = './',bias_directory='./',
                 verbose=True, return_files=False):
    """Apply master bias to MEF data."""
    ##
    ## Apply MEF master bias to a set of object data.
    ##
    ## Notes:
    ##   image_collection -- A ccdproc-style image collection OR file list
    ##   bias_file        -- Filename for bias.
    ##   image_directory  -- Where to store the output images (default ./)
    ##   input_directory  -- Where to find the object data (default ./)
    ##   bias_directory   -- Where to find the bias image (default ./)
    ##   return_files     -- Return a list of files that were flattened (default False)


    # Hardwire number of CCDs in LBC
    lbc_chips = [1,2,3,4]
    num_lbc_chips = np.size(lbc_chips)

    # Process the input data list
    if isinstance(image_collection,ImageFileCollection):
        # image_collection is a CCDPROC-style ImageFileCollection
        image_files = image_collection.files
    elif isinstance(image_collection,list):
        # image_collection is a list of filenames.
        image_files = image_collection
    else:
        # image_collection is something else. Let's guess that this will work.
        image_files = image_collection

    # Construct the bias file name:
    if bias_file == None:
        zero_file = 'zero.flat'
    else:
        zero_file = bias_file

    # Hold a list of files that get corrected
    zero_corrected = []

    # Load bias CCDs into a list.
    # We do this here, as it's faster than repeating this for every chip of every file.
    zero_chips = []
    for chip in lbc_chips:
         zerochip = CCDData.read(bias_directory+zero_file, chip, unit=u.adu)
         zero_chips.append(zerochip)
    if verbose:
         print('Reading bias frame {0}'.format(zero_file))


    for file in image_files:
        # Set up the output HDU list
        # Capture the 0th header
        base_header = fits.getheader(file)
        master_hdu = fits.PrimaryHDU(header=base_header)
        # Start output HDU list:
        output_hdu = fits.HDUList([master_hdu])

        # Loop through the chips
        for chip in lbc_chips:
            # Create the CCDData version of this chip
            image = CCDData.read(input_directory+file, chip, unit=u.adu)
            # Apply the flat
            image_zeroed = ccdproc.subtract_bias(image,zero_chips[chip-1])

            # Convert to 32 bit
            image_zeroed.data=image_zeroed.data.astype('float32')

            # Append the bias'd data into output HDU:
            output_hdu.append(image_zeroed.to_hdu()[0])

        # Create the output file name
        #   - Get rid of the overscan or zero labels and the .fits extension.
        #   - Add the _flat tag.
        output_filename = file.replace('_over','').replace('.fits','_zero.fits')

        # Write the output flat-fielded data
        output_hdu.writeto(image_directory + output_filename, overwrite=True)
        zero_corrected.append(output_filename)

        if verbose:
            print('Bias subtracted {0} to {1}.'.format(file,output_filename))

    if return_files == True:
        return zero_corrected


def make_flatfield(image_collection, filter_name=None, simple_masks=False,
                   image_directory='./',raw_directory='./raw/',
                   verbose=True):
    """Make a flat field image for a collection of images."""
    ##
    ## Create master flat field image for a specific flat
    ##
    ## Notes:
    ##   image_collection -- A ccdproc-style image collection.
    ##   filter_name      -- Filter for flat combination. Default (None) is to
    ##                       determine from image_collection.
    ##   make_masks       -- Make image masks based on flat fields.
    ##   image_directory  -- Where to store the output images.
    ##   raw_directory    -- Where to find the raw data


    # TODO Update make_flatfield to allow for _zero suffix if we've done bias subtraction.
    # TODO Update to step through all the filters in a list.
    # TODO Implement bias subtraction in flatfields.

    # Hardwire the number of CCDs in LBC
    lbc_chips = [1,2,3,4]
    num_lbc_chips = np.size(lbc_chips)

    # If the filter_name variable is passed, use only that filter. Otherwise, use
    #   the first in the list. That is, select the filter beforehand by in image_collection
    #   or select it by passing filter_name.
    if filter_name == None:
        # Create the output file name
        filter_name = image_collection.summary['filter'][0]

    # Report:
    if verbose:
        print('----- MAKE_FLATFIELD: making flat for {0}. -----'.format(filter_name) )

    # Create output file name:
    flat_output_name = 'flat.' + filter_name + '.fits'
    mask_output_name = 'mask.' + filter_name + '.fits'

    # Pull out the flat field images from our collection
    flt_files = image_collection.files_filtered(imagetyp='flat',filter=filter_name)

    # Create 2D list
    # (from https://stackoverflow.com/questions/7745562/appending-to-2d-lists-in-python)
    flat_list = [[] for i in range(num_lbc_chips)]
    num_flat_images = np.zeros(num_lbc_chips,dtype=np.int)

    # Loop through the files
    for filename in flt_files:
        # This will serve to hold all of the master flats before writing
        if filename == flt_files[0]:
            # Set up the output HDU list
            # Capture the 0th header
            base_header = fits.getheader(raw_directory+filename)
            master_hdu = fits.PrimaryHDU(header=base_header)

            # Start output HDU list:
            output_hdu = fits.HDUList([master_hdu])

            mask_hdu = fits.HDUList([master_hdu])
            mask_list = [[] for i in range(num_lbc_chips)]
            num_mask_images = np.zeros(num_lbc_chips)

        # Loop through the chips
        for chip in lbc_chips:
            # Create the CCDData version of this chip
            ccd = CCDData.read(raw_directory+filename, chip,
                               unit=u.adu)
            ccd_filter = ccd.header['filter']

            # Fit, subtract overscan
            poly_model = models.Polynomial1D(4)
            ccd = ccdproc.subtract_overscan(ccd, #median=True,
                                        overscan_axis=1,
                                        model = poly_model,
                                        fits_section=ccd.header['BIASSEC'])
            # Trim the image
            ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])

            # Convert to 32 bit
            ccd.data=ccd.data.astype('float32')

            # Check for saturation, matching filters.
            # **Individual chip headers only include 8 letter filter names**
            if (np.average(ccd.data) <= 62000.) and \
                    (ccd_filter == filter_name[:len(ccd_filter)]):
                # If appropriate, add the flat to our list
                if (chip == 1) & (verbose == True):
                    print('Adding {0} to {1} flatfield'.format(
                            filename,filter_name))
                # Cycle counter
                num_flat_images[chip-1] += 1
                # Add image to list
                flat_list[chip-1].append(ccd)

                # Now do the same for the masks
                if num_mask_images[chip-1] <= 1:
                    mask_list[chip-1].append(ccd)
                    num_mask_images[chip-1] += 1

    # After looping files, chips:
    # Create master flat from overscan-subtracted images
    for chip in lbc_chips:
        # Grab images for this chip alone
        chip_list = flat_list[chip-1]

        # Set up the scale factors.
        #    Normalizes all images so the rejection techn./median combine work.
        flat_scale = []
        for j in np.arange(num_flat_images[chip-1]):
            # Select section in middle of chip. No need to do whole chip.
            section = chip_list[j].data[500:1500,2000:2500]
            flat_scale.append(np.average(1./section))

        # TODO Weight the images to avoid adding too much noise. Prob weight by inverse sqrt.
        # Use ccdproc.combine to combine flat images
        master_flat = ccdproc.combine(chip_list,
                                      scale=flat_scale,
                                      method='median',sigma_clip=True,
                                      sigma_clip_high_thresh=3.,
                                      sigma_clip_low_thresh=3.)

        # Convert to 32 bit
        master_flat.data=master_flat.data.astype('float32')

        # Remove unneeded header keywords. Makes this consistent
        #   with IRAF treatment.
        del master_flat.header['trimsec']
        del master_flat.header['biassec']
        del master_flat.header['datasec']
        del master_flat.header['bzero']    # This is sometimes different than our value

        # Append the flat for this chip to the output HDU
        output_hdu.append(master_flat.to_hdu()[0])

        # Create the ratio image
        #   **Don't calculate the masks unless directed to do so! too long...**
        ratio = mask_list[chip - 1][0].divide(mask_list[chip - 1][1])

        # Calculate mask for this chip.
        if not simple_masks:
            print('Calculating base pixel mask for {0} chip {1}'.format(
                                                                filter_name,chip))

            # Calculate mask. This is _slow_ due to the median calculation.
            mask_out = ccdproc.ccdmask(ratio,findbadcolumns=True)
            # The *1 makes it an int array.
            # mask_out.data = mask_out.data * np.ones(1,dtype=np.float32)
            mask_out.data = mask_out.data * np.ones(1)
            mask_hdu.append(mask_out.to_hdu()[0])

        if simple_masks:
            print('Creating simple pixel mask for {0} chip {1}'.format(
                filter_name, chip))
            # If not calculating full mask, for now set everything to good (0).
            mask_temp = ratio.copy()
            mask_temp.data = mask_temp.data * np.zeros(1)
            mask_hdu.append(mask_temp.to_hdu()[0])

        # Same basic header for the masks
        temp_header = output_hdu[0].header
        mask_hdu[0].header = temp_header


    # TODO Include list of flatfield images combined in output header
    # Put the combined image number in the top header, as well.
    output_hdu[0].header['ncombine'] = num_flat_images[0]

    # Write the master flat
    output_hdu.writeto(image_directory+flat_output_name, overwrite=True)

    # Write the mask
    mask_hdu.writeto(image_directory+mask_output_name, overwrite=True)

    # Report:
    if verbose:
        print("\nCreated {0} flatfield frame {1}.".format(
            filter_name, flat_output_name))
        print("Created {0} mask base {1}.\n".format(
            filter_name, mask_output_name))


def do_flatfield(image_collection, flat_file=None, filter_names = None,
                 image_directory='./',input_directory = './',
                 flat_directory='./',
                 verbose=True, return_files=False):
    """Apply flat fields to MEF data."""
    ##
    ## Apply MEF flat fields to a set of object data.
    ##
    ## Notes:
    ##   image_collection -- A ccdproc-style image collection.
    ##   flat_file        -- Filename for a specific flat to be applied.
    ##   filter_names     -- List of filters for images to be flattened. If None, it will be generated.
    ##   image_directory  -- Where to store the output images (default ./)
    ##   input_directory  -- Where to find the object data (default ./)
    ##   flat_directory   -- Where to find the flat fields (default ./)
    ##   return_files     -- Return a list of files that were flattened (default False)

    import warnings
    from astropy.utils.exceptions import AstropyWarning,AstropyUserWarning
    warnings.filterwarnings('ignore', category=AstropyWarning, append=True)
    warnings.filterwarnings('ignore', category=AstropyUserWarning, append=True)

    # The pre-flatfield images will be put in a data/ directory:
    datadir = 'data/'
    if not os.path.exists(datadir):
        os.makedirs(datadir)

    # Hardwire number of CCDs in LBC
    lbc_chips = [1,2,3,4]
    num_lbc_chips = np.size(lbc_chips)

    # Identify the filters for which to do flats:
    if filter_names == None:
        filter_names = np.unique(image_collection.summary['filter'])

    # Check to see if the user passed a flat-field file and a list with multiple filters
    if (flat_file != None) & (np.size(filter_names) > 1):
        print('WARNING: a single flatfield file was passed, while multiple filters are present.')

    # Hold a list of files that get flattened.
    flattened_files = []

    # Loop over the available filters
    for filter in filter_names:
        # Work out which flatfield file we're working with
        if flat_file == None:
            flat_filename = 'flat.'+filter+'.fits'
        else:
            flat_filename = flat_file

        # For each filter, read the flat field just once, stack into a list.
        # We do this here, as it's faster than repeating this for every chip of every file.
        flatfield_chips = []
        for chip in lbc_chips:
             flatchip = CCDData.read(flat_directory+flat_filename, chip,
                                     unit=u.adu)
             flatfield_chips.append(flatchip)
        if verbose:
             print('Reading flatfield {0}'.format(flat_filename))

        # Select the files to be normalized with the current filter
        image_files = image_collection.files_filtered(filter=filter)

        for file in image_files:
            # Set up the output HDU list
            # Capture the 0th header
            base_header = fits.getheader(file)
            master_hdu = fits.PrimaryHDU(header=base_header)
            # Start output HDU list:
            output_hdu = fits.HDUList([master_hdu])

            # Loop through the chips
            for chip in lbc_chips:
                # Create the CCDData version of this chip
                image = CCDData.read(input_directory+file, chip,
                                     unit=u.adu)
                # Apply the flat
                image_normed = ccdproc.flat_correct(image,
                            flatfield_chips[chip-1],
                            min_value=0.1,
                            norm_value=np.median(flatfield_chips[chip-1]))

                # Convert to 32 bit
                image_normed.data=image_normed.data.astype('float32')

                # Append the flattened data into output HDU:
                output_hdu.append(image_normed.to_hdu()[0])

            # Create the output file name
            #   - Get rid of the overscan or zero labels and the .fits extension.
            #   - Add the _flat tag.
            output_filename = file.replace('_over','').replace('_zero','')
            output_filename = output_filename.replace('.fits','_flat.fits')

            # Write the output flat-fielded data
            output_hdu.writeto(image_directory+output_filename,overwrite=True)

            # Append the flattened image to our final list.
            flattened_files.append(output_filename)

            # Move pre-flatfield file to data directory
            cmd = 'mv '+file+' '+datadir
            mvover = Popen(shlex.split(cmd),
                  close_fds=True)
            mvover.wait()

            if verbose:
                print('Flattened {0} to {1}.'.format(file,output_filename))

    if return_files == True:
        return flattened_files


def make_bpm():
    """Create the bad pixel masks needed for...sextractor??"""

    print('')

def clean_cosmic(ccd, mbox=15, rbox=15, gbox=11, sigclip=5,
                 cleantype="medmask", cosmic_method='lacosmic'):
    # From ReduceCCD, cleanCOSMIC of R. Garcia-Benito
    #  (https://github.com/rgbIAA/reduceccd)

    ##
    ## Currently not enabled
    ##

    ctype = cosmic_method.lower().strip()
    ctypes = ['lacosmic', 'median']
    if not ctype in ctypes:
        print ('>>> Cosmic ray type "%s" NOT available [%s]' % (ctype, ' | '.join(ctypes)))
        return
    if ctype == 'lacosmic':
        ccd = ccdproc.cosmicray_lacosmic(ccd, sigclip=sigclip, cleantype=cleantype)
    elif ctype == 'median':
        ccd = ccdproc.cosmicray_median(ccd, mbox=mbox, rbox=rbox, gbox=gbox)
    if isinstance(ccd, CCDData):
        ccd.header['COSMIC'] = ctype.upper()
    return ccd


def make_targetdirectories(image_collection, image_directory='./',
                           object_names=None, verbose=True):
    """ Make directories to store data from individual targets and move the data for those targets into the directories.

    :param image_collection:
    :param image_directory:
    :param object_names:
    :param verbose:
    :return: object_directories, filter_directories
    """

    # If we specify objects, use only those objects
    if object_names == None:
        object_names = np.unique(image_collection.summary['object'])


    # Return a list of the directories
    object_directories = []
    filter_directories = []

    # Loop through the objects
    for object in object_names:

        # Name of the target directory to be created
        dirname= image_directory+object+'/'
        object_directories.append(dirname)

        # Test that the directory doesn't already exist
        if not os.path.lexists(dirname):
            os.makedirs(dirname)

            if verbose == True:
                print('Creating directory for object {0}.'.format(object))
        else:
            if verbose == True:
                print('Directory exists for object {0}.'.format(object))


        # Move the data for each object into its directory.
        object_files = image_collection.files_filtered(object=object)
        for fl in object_files:
            cmd = 'mv {0} {1}'.format(fl,dirname)
            crap = call(cmd,shell=1)


        # Create filter-specific directories and fill them. For now this just assumes that every
        # object has the same filter list. This is fine, since we're normally going through this step
        # on a filter-by-filter basis anyway (see proc.py).

        # The following should select only the filters appropriate for this object:
        keywds = ['object','filter']
        image_collectionObj = ImageFileCollection(dirname, keywords=keywds,
                                    filenames = (image_collection.files_filtered(object=object)).tolist())
        # Now select the unique filters for this objects
        filters = image_collectionObj.values('filter',unique=True)

        # Step through each filter, creating directories and moving files.
        for filter in filters:
            filter_dirname = dirname+filter+'/'

            # Test that the directory doesn't already exist
            if not os.path.lexists(filter_dirname):
                os.makedirs(filter_dirname)
                filter_directories.append(filter_dirname)

                if verbose == True:
                    print('Creating {1} directory for object {0}.'.format(object,filter))
            else:
                if verbose == True:
                    print('Directory exists for {1} for object {0}.'.format(object,filter))

            filter_files = image_collection.files_filtered(object=object, filter=filter)
            for fltfl in filter_files:
                cmd = 'mv {0} {1}'.format(dirname+fltfl, filter_dirname)
                crap = call(cmd, shell=1)

    return object_directories, filter_directories


def extract_chips(filter_directories,
                  verbose=True, return_files = False):
    """Extract individual chips from flat-fielded data in preparation for
    combination using sextractor/scamp/swarp.

    :param filter_directories:
    :param verbose:
    :param return_files:
    :return:
    """
    # If user enters just a single directory:
    if np.size(filter_directories) == 1 & ~isinstance(filter_directories,list):
        filter_directories = [filter_directories]

    # The pre-chip-extraction images will be put in a data/ directory:
    datadir = 'data/'
    if not os.path.exists(datadir):
        os.makedirs(datadir)

    # Hardwire the num. of LBC chips
    lbc_chips = [1,2,3,4]
    num_lbc_chips = np.size(lbc_chips)

    # ImageFileCollection keywords
    keywds = ['object', 'filter', 'exptime', 'objra', 'objdec']

    # Create a list of files to split apart.
    input_filenames = []
    for fltdr in filter_directories:
        fls = glob(fltdr + '*_flat.fits')
        for fl in fls:
            input_filenames.append(fl)

    chip_files = []
    # Loop through the files
    for filename in input_filenames:
        # Fill the 0th header with the right info
        base_header = fits.getheader(filename)
        master_hdu = fits.PrimaryHDU(header=base_header)

        # We will only house one extension in these files:
        master_hdu.header['NEXTEND'] = 1

        # Loop through the chips
        for chip in lbc_chips:
            # Create the output HDU list:
            output_hdu = fits.HDUList([master_hdu])

            # Create the output filename.
            filesuffix = '_{0}.fits'.format(chip)
            output_filename = filename.split('_flat.fits')[0] + filesuffix

            # This will serve to hold all of the master flats before writing
            # Fill the 0th header with the right info
            # base_header = fits.getheader(filename)
            #
            # master_hdu = fits.PrimaryHDU(header=base_header)
            # output_hdu = fits.HDUList([master_hdu])
            # Add to output HDU
            # output_hdu.append(ccd.to_hdu()[0])


            # Read the CCDData version of this chip
            ccd = CCDData.read(filename, chip,
                               unit=u.adu)
            # Create the new HDU
            output_hdu.append(ccd.to_hdu()[0])

            # Write the data
            output_hdu.writeto(output_filename, overwrite=True)

            # Keep track of what files we've written.
            chip_files.append(output_filename)

            # Report:
            if verbose:
                print("Created {0}".format(output_filename))

        # Move pre-flatfield file to data directory
        cmd = 'mv '+filename+' '+datadir
        mvflt = Popen(shlex.split(cmd),
              close_fds=True)
        mvflt.wait()

    # Return the corrected filenames if requested (True is default).
    if return_files == True:
        return chip_files
