import numpy as np
from astropy.io import fits

import ccdproc
from ccdproc import  ImageFileCollection,CCDData

import astropy.units as u
import astropy.constants as c

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
# TODO More general MEF flat fielding
# TODO make_bias
# TODO Logging system

def do_overscan(over_files):
    ##
    ## Fit, subtract overscan
    ##
    """Remove overscan and trim images."""

    # Hardwire the num. of LBC chips
    num_lbc_chips = 4

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
            ccd = ccdproc.subtract_overscan(ccd, median=True,
                                            overscan_axis=1,
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


def do_flatfield(image_collection,flat_file=None):
    """Apply flat fields to MEF data."""

    # Hardwire number of CCDs in LBC
    num_lbc_chips = 4
    # Where are the flatfield data stored?
    flat_directory = './'

    # Loop over the filters represented (in case there are multiple)
    filter_names = np.unique(ic1.summary['filter'])

    # Check to see if the user passed a flat-field file and a list with multiple filters
    if (flat_file != None) & (np.size(filter_names) > 1):
        print('WARNING: a single flatfield file was passed, while multiple filters are present.')

    for filter in filter_names:
        # Work out which flatfield file we're working with
        if flat_file == None:
            flat_filename = flat_directory+'flat.'+filter+'.fits'
        else:
            flat_filename = flat_file

        # For each filter, read the flat field just once, stack into a list.
        # We do this here, as it's faster than repeating this for every chip of every file.
        flatfield_chips = []
        for chip in np.arange(1, num_lbc_chips + 1):
             flatchip = CCDData.read(flat_filename, chip, unit=u.adu)
             flatfield_chips.append(flatchip)

        # Select the files to be normalized with the current filter
        image_files = image_collection.files_filtered(filter=filter)
        num_images = np.size(image_files)

        for file in image_files:
            master_hdu = fits.open(file)

            # It's likely slower to reload the flatfield for each
            # Loop through the chips
            for chip in np.arange(1, num_lbc_chips + 1):
                # Create the CCDData version of this chip
                image = CCDData.read(filename, chip, unit=u.adu)
                # Apply the flat
                image = ccdproc.flat_correct(image, flatfield_chips[chip-1])
                # Put the flattened data and header into the master HDU
                master_hdu[chip].data = image.data
                master_hdu[chip].header = image.header

            # Create the output file name
            output_filebase = (file.split('.fits'))[0].split('_over')[0]
            output_filename = output_filebase+'_flat.fits'
            # Write the output flat-fielded data
            master_hdu.writeto(output_filename,overwrite=True)

def make_bias(image_collection):
    print('')

def make_flatfield(image_collection):
    """Make a flat field image for the collection of images passed.

    So far this assumes the filter is the same for all the images..."""

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

    # Loop through the files
    for j in np.arange(np.size(flt_files)):
        filename = flt_files[j]

        # This will serve to hold all of the master flats before writing
        if j == 0:
            master_hdu = fits.open(filename)

        # Loop through the chips
        for chip in np.arange(1,num_lbc_chips+1):
            # Create the CCDData version of this chip
            ccd = CCDData.read(filename, chip, unit=u.adu)
            ccd_filter = ccd.header['filter']

            # Fit, subtract overscan
            ccd = ccdproc.subtract_overscan(ccd, median=True,
                                        overscan_axis=1,
                                        fits_section=ccd.header['BIASSEC'])
            # Trim the image
            ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])

            # Check for saturation, matching filters.
            # **Individual chip headers only include 8 letter filter names**
            if (np.average(ccd) <= 55000.) and \
                    (ccd_filter == filter_name[:len(ccd_filter)]):
                # If appropriate, add the flat to our list
                flat_list[chip-1].append(ccd)

    # After looping files, chips:
    # Create master flat from subtracted images
    for chip in np.arange(1,num_lbc_chips+1):
        master_flat = ccdproc.combine(flat_list[chip-1], method='median')
        master_hdu[chip].data = master_flat.data

    # Write the master flat
    master_hdu.writeto(flat_output_name, overwrite=True)
    master_hdu.close()

    # Report:
    print("Created {0}".format(flat_output_name))

#     return master_flat

def make_bpm():
    """Create the bad pixel masks needed for...sextractor??"""

    print('')


def proc_lbc():
    """Process a directory of LBC data"""
    print('')






