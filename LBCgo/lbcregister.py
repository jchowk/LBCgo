import os
import numpy as np
from subprocess import Popen, DEVNULL
import shlex
from glob import glob
from astropy.io import fits
from ccdproc import  ImageFileCollection

def go_sextractor(inputfile,
                configfile=None,
                paramfile = None,
                convfile = 'default.conv',
                nnwfile = 'default.nnw',
                verbose=True, clean=True):
    """
    """
    # path = os.path.abspath(amodule.__file__)

    if configfile == None:
        configfile = 'sextractor.lbc.conf'

    # TODO: Use default if output.param file doesnt exist
    if paramfile == None:
        paramfile = 'sextractor.lbcoutput.param'
        file_exists = os.path.exists(paramfile)
        if not file_exists:
            print("Warning: SEXTRACTOR paramater file {0} "
                  "does not exist".format(paramfile))
            return None

    # Test for convolution file
    file_exists = os.path.exists(convfile)
    if not file_exists:
        print("Warning: SEXTRACTOR convolution file {0} "
              "does not exist".format(convfile))
        return None

    # Test for NNW file
    file_exists = os.path.exists(nnwfile)
    if not file_exists:
        print("Warning: SEXTRACTOR convolution file {0} "
              "does not exist".format(nnwfile))
        return None

    # Base filename:
    filebase = inputfile.replace('.fits','')

    # Source extractor catalog suffix
    outputsuffix = '.cat'
    outputcatalog = filebase+outputsuffix


    cmd_flags = ' -c '+ configfile + \
        ' -CATALOG_NAME '+outputcatalog + \
        ' -CATALOG_TYPE FITS_LDAC'+ \
        ' -DETECT_THRESH 5.0 -ANALYSIS_THRESH 8.0'+ \
        ' -PARAMETERS_NAME '+paramfile

    cmd = 'sex '+inputfile+cmd_flags
    # print(cmd)

    try:
        if verbose:
            print(cmd)
            sextract = Popen(shlex.split(cmd),
                             close_fds=True)
        else:
            sextract = Popen(shlex.split(cmd),
                             stdout=DEVNULL,
                             stderr=DEVNULL,
                             close_fds=True)
    except Exception as e:
        print('Whoops: source Extractor call:', (e))
        return None

    sextract.wait()


def go_scamp(inputfile,
             astrometric_catalog='GAIA-DR2',
             astrometric_method = 'exposure',
             num_iterations = 3,
             configfile=None,
             verbose=True, clean=True):

    """ Run Astromatic.net code SCAMP to calculate astrometric
     solution on an image.
    """

    # Make sure the input file is a SEXTRACTOR catalog:
    inputfile = inputfile.replace('.fits','.cat')

    # Make sure we have a configuration file:
    if configfile == None:
        # TODO: replace this with default config file for LBCgo.
        configfile = 'scamp.lbc.conf'

    # Using only a single iteration of SCAMP doesn't do well enough. Force at
    # least two iterations:
    if num_iterations < 2:
        print('WARNING: Use at least 2 SCAMP iterations. Setting num_iterations = 2...')
        num_iterations = 2

    SN_thresholds = '5,50'
    for scmpiter in np.arange(num_iterations):
        if scmpiter == 0:
            degree = '3'
            mosaic_type = 'LOOSE'
            pixscale_maxerr = '1.2'
            position_maxerr = '1'
            posangle_maxerr = '5.0'
            crossid_radius = '7.5'
            aheader_suffix = '.ahead'
        elif scmpiter == 1:
           degree = '3'
           mosaic_type = 'FIX_FOCALPLANE'
           pixscale_maxerr = '1.1'
           position_maxerr = '0.1'
           posangle_maxerr = '3.0'
           crossid_radius = '5.0'
           aheader_suffix = '.head'
        elif scmpiter == 2:
           degree = '3'
           mosaic_type = 'FIX_FOCALPLANE'
           pixscale_maxerr = '1.05'
           position_maxerr = '0.05'
           posangle_maxerr = '1.0'
           crossid_radius = '2.5'
           aheader_suffix = '.head'
        else:
           degree = '3'
           mosaic_type = 'FIX_FOCALPLANE'
           pixscale_maxerr = '1.05'
           position_maxerr = '0.025'
           posangle_maxerr = '1.0'
           crossid_radius = '1.5'
           aheader_suffix = '.head'

        cmd_flags = ' -c '+ configfile + \
            ' -PIXSCALE_MAXERR '+pixscale_maxerr+ \
            ' -POSANGLE_MAXERR '+posangle_maxerr+ \
            ' -POSITION_MAXERR '+position_maxerr+ \
            ' -DISTORT_DEGREES '+degree+ \
            ' -MOSAIC_TYPE '+mosaic_type+ \
            ' -ASTREF_CATALOG '+astrometric_catalog+ \
            ' -AHEADER_SUFFIX '+aheader_suffix+ \
            ' -CROSSID_RADIUS '+crossid_radius+\
            ' -STABILITY_TYPE INSTRUMENT'+\
            ' -SN_THRESHOLDS '+SN_thresholds

        if astrometric_method == 'exposure':
            cmd_flags.replace('INSTRUMENT','EXPOSURE')

        # Create the final command:
        cmd = 'scamp '+inputfile+cmd_flags

        try:
            if verbose:
                print('########### SCAMP iteration {0} for {1} '
                '########### \n'.format(scmpiter+1,
                inputfile.replace('.cat','')))
                # Diagnostics
                print(cmd)

                scamp = Popen(shlex.split(cmd),
                                   close_fds=True)
            else:
                scamp = Popen(shlex.split(cmd),
                                   stdout=DEVNULL,
                                   stderr=DEVNULL,
                                   close_fds=True)
        except Exception as e:
            print('Whoops: source Extractor call:', (e))
            return None

        scamp.wait()

    if clean:
        # Clean the GAIA catalog files
        catfiles = glob('GAIA*cat')
        for ctfls in catfiles:
            cmd = 'rm '+ctfls
            clean_dir = Popen(shlex.split(cmd),
                close_fds=True)
            clean_dir.wait()



def go_swarp(inputfiles, output_filename = None,
                configfile = None,
                verbose = True):
    """Do SWARP"""

    # Make sure we have a configuration file:
    if configfile == None:
        # TODO: replace this with default config file in LBCgo directories.
        configfile = 'swarp.lbc.conf'

    # Gather some information about the input files
    keywds = ['object', 'filter', 'exptime', 'imagetyp', 'propid', 'lbcobnam',
                  'airmass', 'HA', 'objra', 'objdec']
    ic_swarp = ImageFileCollection('./', keywords=keywds,
                                    filenames = inputfiles)

    # Check all are the same filter:
    # TODO: Set exception if not all are same filter rather than just bail.
    fltrs = ic_swarp.values('filter',unique=True)
    if np.size(fltrs) != 1:
        print('Warning: not all files in SWARP call use the same filter.')
        return None

    # Calculate mean airmass (weighted by exposure time)
    exp_airmass = np.array(ic_swarp.values('airmass'))
    exp_time = np.array(ic_swarp.values('exptime'))
    airmass = np.average(exp_airmass,weights=exp_time)

    # from IPython import embed ; embed()

    # For now grab the information from the first header:
    if output_filename == None:
        imhead = fits.getheader(inputfiles[0])
        # Shorten the filter names used:
        filter_text = imhead['FILTER']
        filter_text = filter_text.replace('-SLOAN','').replace('-BESSEL','').replace('SDT_Uspec','Uspec')
        # Create final output filename
        # TODO: Check if file exists.
        output_filename = (imhead['object']).\
            replace(' ','')+'.'+filter_text+'.mos.fits'

    # Rename the weight image
    weight_filename = output_filename.replace('.mos.fits','.mos.weight.fits')

    # Create the list of input files:
    inputfile_text = ''
    for fl in inputfiles: inputfile_text = inputfile_text+' '+fl

    cmd_flags = ' -c '+configfile+ \
        ' -IMAGEOUT_NAME '+ output_filename + \
        ' -WEIGHTOUT_NAME '+ weight_filename + \
        ' -WEIGHT_TYPE NONE '+ \
        ' -HEADER_SUFFIX ".head"'+ \
        ' -FSCALE_KEYWORD NONE -FSCALE_DEFAULT 1.0 '+\
        ' -CELESTIAL_TYPE EQUATORIAL -CENTER_TYPE ALL '+\
        ' -COMBINE_BUFSIZE 4096 ' +\
        ' -COPY_KEYWORDS '+\
        ' OBJECT,OBJRA,OBJDEC,OBJEPOCH,PROPID,PI_NAME,'+\
        'FILTER,SATURATE,RDNOISE,GAIN,EXPTIME,AIRMASS,TIME-OBS'

    # Create the final command:
    cmd = 'swarp ' + inputfile_text + cmd_flags
    # debug
    print(cmd)
    try:
        if verbose:
            swarp = Popen(shlex.split(cmd),
                          close_fds=True)
        else:
            swarp = Popen(shlex.split(cmd),
                          stdout=DEVNULL,
                          stderr=DEVNULL,
                          close_fds=True)
    except Exception as e:
        print('Whoops: SWARP call:', (e))
        return None

    swarp.wait()

    # Add airmass to header:
    fits.setval(output_filename,'AIRMASS',value=airmass)

# def go_imagequality(inputfile,
#                 configfile=None,
#                 paramfile = None,
#                 convfile = 'default.conv',
#                 nnwfile = 'default.nnw',
#                 verbose=True, clean=True):
#     """
#     """
#

def go_register(filter_directories,
                lbc_chips = [1,2,3,4],
                do_sextractor=True,
                do_scamp=True,
                do_swarp=True,
                astrometric_catalog='GAIA-DR2',
                scamp_iterations = 3):

    # TODO: Add the sextractor, scamp, swarp parameters for input.

    # If user enters just a single directory:
    if np.size(filter_directories) == 1 & ~isinstance(filter_directories,list):
        filter_directories = [filter_directories]

    for j in np.arange(np.size(filter_directories)):
        # Make sure the input directories have trailing slashes:
        drctry = filter_directories[j]
        if filter_directories[j][-1] != '/':
            filter_directories[j] += '/'


    # Loop through each of the filter directories:
    for fltdr in filter_directories:
        input_filenames = []

        # Only include the chips we want in the final image:
        for chp in lbc_chips:
            fls = glob(fltdr + '*_'+str(chp)+'.fits')
            for fl in fls:
                input_filenames.append(fl)

        # Loop through the files
        # go_sextractor = find sources
        # go_scamp = calculate astrometry
        for filename in input_filenames:
            # Find sources for alignment
            if do_sextractor:
                go_sextractor(filename)
            # Calculate the astrometry
            if do_scamp:
                go_scamp(filename, astrometric_catalog=astrometric_catalog,
                         num_iterations = scamp_iterations)

        # Stitch together the images
        # go_swarp = reproject and coadd images
        if do_swarp:
            go_swarp(input_filenames)
