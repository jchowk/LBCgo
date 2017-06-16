import numpy as np

from astropy.io import fits
import astropy.units as u
import astropy.constants as c

from astropy.modeling import models

import ccdproc
from ccdproc import  ImageFileCollection,CCDData

from collections import OrderedDict

# Holder: information useful for CR rejection
# lacos_im = OrderedDict()
# lacos_im['gain']=1.75
# lacos_im['readn']=12.0
# lacos_im['sigclip']=4.5
# lacos_im['sigfrac']=0.3
# lacos_im['objlim']=1.0
# lacos_im['niter']=2

# TODO Create some common block variables for information? num_lbc_chips
# TODO Create option for cleaning intermediate steps
# TODO Do gain correction
# TODO More general MEF flat fielding?
# TODO Logging system

# TODO ** make_bias
# TODO ** do_bias

# TODO cosmic ray cleaning


def do_overscan(image_collection,objects_only = True):
    ##
    ## Fit, subtract overscan
    ##
    """Remove overscan and trim images."""

    # Hardwire the num. of LBC chips
    num_lbc_chips = 4

    if objects_only == True:
        over_files = image_collection.files_filtered(imagetyp='object')
    else:
        over_files = image_collection.files

    # Loop through the files
    for j in np.arange(np.size(over_files)):
        filename = over_files[j]
        output_filename = filename.split('.fits')[0]+'_over.fits'

        # This will serve to hold all of the master flats before writing
        output_hdu = fits.open(filename)

        # Loop through the chips
        for chip in np.arange(1, num_lbc_chips + 1):
            # Create the CCDData version of this chip
            ccd = CCDData.read(filename, chip, unit=u.adu)
            ccd_filter = ccd.header['filter']

            # Fit, subtract overscan
            poly_model = models.Polynomial1D(4)
            ccd = ccdproc.subtract_overscan(ccd, #median=True,
                                        overscan_axis=1,
                                        model = poly_model,
                                        fits_section=ccd.header['BIASSEC'])
            # Trim the image
            ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])
            # Put the data back into the hdu:
            output_hdu[chip].data = ccd.data
            output_hdu[chip].header = ccd.header

        # Write the data
        # TODO Do we want to write these files or just keep them in memory to return to main?
        output_hdu.writeto(output_filename, overwrite=True)
        output_hdu.close()

        # Report:
        print("Created {0}".format(output_filename))


def do_flatfield(image_collection,flat_file=None,flat_directory='./'):
    """Apply flat fields to MEF data."""

    # Hardwire number of CCDs in LBC
    num_lbc_chips = 4

    # Identify the filters present in the image collection
    filter_names = np.unique(image_collection.summary['filter'])

    # Check to see if the user passed a flat-field file and a list with multiple filters
    if (flat_file != None) & (np.size(filter_names) > 1):
        print('WARNING: a single flatfield file was passed, while multiple filters are present.')

    # Loop over the available filters
    for filter in filter_names:
        # Work out which flatfield file we're working with
        if flat_file == None:
            flat_filename = flat_directory+'flat.'+filter+'.fits'
        else:
            flat_filename = flat_file

        # For each filter, read the flat field just once, stack into a list.
        # We do this here, as it's faster than repeating this for every chip of every file.
        flatfield_chips = []
        print('Reading flatfield {0}'.format(flat_filename))
        for chip in np.arange(1, num_lbc_chips + 1):
             flatchip = CCDData.read(flat_filename, chip, unit=u.adu)
             flatfield_chips.append(flatchip)

        # Select the files to be normalized with the current filter
        image_files = image_collection.files_filtered(filter=filter)
        num_images = np.size(image_files)

        for file in image_files:
            master_hdu = fits.open(file)
            # Loop through the chips
            for chip in np.arange(1, num_lbc_chips + 1):
                # Create the CCDData version of this chip
                image = CCDData.read(file, chip, unit=u.adu)
                # Apply the flat
                image_normed = ccdproc.flat_correct(image, flatfield_chips[chip-1])
                # Put the flattened data and header into the master HDU
                master_hdu[chip].data = image_normed.data
                master_hdu[chip].header = image_normed.header

            # Create the output file name
            output_filebase = (file.split('.fits'))[0].split('_over')[0]
            output_filename = output_filebase+'_flat.fits'
            # Write the output flat-fielded data
            master_hdu.writeto(output_filename,overwrite=True)
            print('Flattened {0} to {1}.'.format(file,output_filename))

def make_bias(image_collection):
    """Make a master bias image for a collection of images."""

    # Pull out the flat field images from our collection
    zero_files = image_collection.files_filtered(imagetyp='zero')
    # Create the output file name
    filter_name = image_collection.summary['filter'][0]
    zero_output_name = 'zero.fits'

    # Hardwire the number of CCDs in LBC
    num_lbc_chips = 4

    # Create 2D list
    # (from https://stackoverflow.com/questions/7745562/appending-to-2d-lists-in-python)
    zero_list = [[] for i in range(num_lbc_chips)]

    # Loop through the files
    for j in np.arange(np.size(zero_files)):
        filename = zero_files[j]

        # This will serve to hold all of the master flats before writing
        if j == 0:
            master_hdu = fits.open(filename)

        # Loop through the chips
        for chip in np.arange(1,num_lbc_chips+1):
            # Create the CCDData version of this chip
            ccd = CCDData.read(filename, chip, unit=u.adu)

            # Fit, subtract overscan
            poly_model = models.Polynomial1D(4)
            ccd = ccdproc.subtract_overscan(ccd, #median=True,
                                        overscan_axis=1,
                                        model = poly_model,
                                        fits_section=ccd.header['BIASSEC'])

            # Trim the image
            ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])

            # Add the current image to the list
            zero_list[chip-1].append(ccd)

    # After looping files, chips:
    # Create master bias from subtracted images
    for chip in np.arange(1,num_lbc_chips+1):
        master_zero = ccdproc.combine(zero_list[chip-1], dtype='float32',
                                      method='median',sigma_clip=True,
                                      sigma_clip_high_thresh=3.,
                                      sigma_clip_low_thresh=3.)
        master_hdu[chip].data = master_zero.data

    # Write the master flat
    master_hdu.writeto(zero_output_name, overwrite=True)
    master_hdu.close()

    # Report:
    print("Created {0} master bias frame.".format(zero_output_name))

def make_flatfield(image_collection, verbose=True):
    """Make a flat field image for a collection of images."""

    # Pull out the flat field images from our collection
    flt_files = image_collection.files_filtered(imagetyp='flat')
    # Create the output file name
    filter_name = image_collection.summary['filter'][0]
    flat_output_name = 'flat.'+filter_name+'.fits'

    # Hardwire the number of CCDs in LBC
    num_lbc_chips = 4

    # Create 2D list
    # (from https://stackoverflow.com/questions/7745562/appending-to-2d-lists-in-python)
    flat_list = [[] for i in range(num_lbc_chips)]
    num_flat_images = 0

    # Loop through the files
    for filename in flt_files:
        # This will serve to hold all of the master flats before writing
        if filename == flt_files[0]:
            master_hdu = fits.open(filename)

        # Loop through the chips
        for chip in np.arange(1,num_lbc_chips+1):
            # Create the CCDData version of this chip
            ccd = CCDData.read(filename, chip, unit=u.adu)
            ccd_filter = ccd.header['filter']

            # Fit, subtract overscan
            poly_model = models.Polynomial1D(3)
            ccd = ccdproc.subtract_overscan(ccd, #median=True,
                                        overscan_axis=1,
                                        model = poly_model,
                                        fits_section=ccd.header['BIASSEC'])
            # Trim the image
            ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])

            # Check for saturation, matching filters.
            # **Individual chip headers only include 8 letter filter names**
            if (np.average(ccd) <= 62000.) and \
                    (ccd_filter == filter_name[:len(ccd_filter)]):
                # If appropriate, add the flat to our list
                if (chip == 1) & (verbose == True):
                    print('Adding {0} to {1} flatfield'.format(
                            filename,filter_name))
                # Cycle counter
                num_flat_images += 1
                # Add image to list
                flat_list[chip-1].append(ccd)

    # After looping files, chips:
    # Create master flat from subtracted images
    for chip in np.arange(1,num_lbc_chips+1):
        # Grab images for this chip alone
        chip_list = flat_list[chip-1]

        # Set up the scale factors
        flat_scale = []
        for j in np.arange(num_flat_images/num_lbc_chips):
            # Select section in middle of chip.
            section = chip_list[j].data[500:1500,2000:2500]
            flat_scale.append(np.median(1./section))

        # Normalize scale factor to keep values close to num of ADU.
        # flat_scale = np.array(flat_scale)/np.mean(flat_scale)

        # Use ccdproc.combine to combine flat images
        master_flat = ccdproc.combine(flat_list[chip-1], dtype=np.float32,
                                      scale=flat_scale,
                                      method='median',sigma_clip=True,
                                      sigma_clip_high_thresh=3.,
                                      sigma_clip_low_thresh=3.)

        # Insert the combined image into the master HDU.
        master_hdu[chip].data = master_flat.data
        # Remove unneeded header keywords. Makes this consistent
        #   with IRAF treatment.
        master_hdu[chip].header['trimsec'] = ''
        master_hdu[chip].header['biassec'] = ''
        master_hdu[chip].header['datasec'] = ''
        # Note how many files are combined
        master_hdu[chip].header['ncombine'] = num_flat_images

    # TODO Include image list for combination
    # Put the combined image number in the top header, as well.
    master_hdu[0].header['ncombine'] = num_flat_images

    # Write the master flat
    master_hdu.writeto(flat_output_name, overwrite=True)
    master_hdu.close()

    # Report:
    if verbose == True:
        print("Created {0} flatfield frame {1}.".format(
            filter_name, flat_output_name))


def make_bpm():
    """Create the bad pixel masks needed for...sextractor??"""

    print('')



def cleanCosmic(ccd, mbox=15, rbox=15, gbox=11, sigclip=5, cleantype="medmask", cosmic_method='lacosmic'):
    # From ReduceCCD of R. Garcia-Benito
    #  (https://github.com/rgbIAA/reduceccd)

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

def proclbc():
    """Process a directory of LBC data"""
    print('')
