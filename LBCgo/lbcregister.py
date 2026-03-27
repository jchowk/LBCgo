import os
import numpy as np
import shutil
from subprocess import Popen, DEVNULL
import shlex
from glob import glob
from astropy.io import fits
from ccdproc import  ImageFileCollection
import astropy.io.votable as votable
import LBCgo


# TODO: Offer an iterative treatment of SCAMP to get desired precision
# TODO: Photometric calibration / selection of filters

def go_sextractor(inputfile,
                configfile=None,
                paramfile = None,
                convfile = None,
                nnwfile = None,
                verbose=True):
    """Run SExtractor on a single chip image to produce a source catalog.

    Detects sources and writes a FITS_LDAC catalog (``<base>.cat``) alongside
    the input file. Default configuration files are read from the LBCgo
    package ``conf/`` directory.

    Parameters
    ----------
    inputfile : str
        Path to the chip FITS image to process.
    configfile : str or None, optional
        Path to a SExtractor configuration file. If None, uses the LBCgo
        default ``sextractor.lbc.conf``. Default: None
    paramfile : str or None, optional
        Path to a SExtractor output parameter file. If None, uses the LBCgo
        default ``sextractor.lbcoutput.param``. Default: None
    convfile : str or None, optional
        Path to a SExtractor convolution kernel file. If None, uses the
        LBCgo default ``default.conv``. Default: None
    nnwfile : str or None, optional
        Path to a SExtractor neural network weights file. If None, uses the
        LBCgo default ``default.nnw``. Default: None
    verbose : bool, optional
        Print the SExtractor command and progress. Default: True

    Raises
    ------
    RuntimeError
        If SExtractor (``sex``) is not found on the system PATH.

    Returns
    -------
    None
        Writes a FITS_LDAC catalog to ``<inputfile_base>.cat``.
        Returns None if a required configuration file is missing.
    """
    # Check if SExtractor is available
    if not shutil.which('sex'):
        raise RuntimeError("SExtractor (sex) is not available. Please install it.")
    
    # from IPython import embed; embed()

    if configfile == None:
        # Use default config file from LBCgo package directories
        configfile = os.path.join(LBCgo.__path__[0], 'conf', 'sextractor.lbc.conf')
        if verbose:
            print("Using default SExtractor configuration file: {0}".format(configfile))
        

    # Use default if output.param file doesnt exist
    if paramfile == None:
        paramfile = os.path.join(LBCgo.__path__[0], 'conf', 'sextractor.lbcoutput.param')
        if verbose:
            print("Using default SExtractor parameter file: {0}".format(paramfile))
        file_exists = os.path.exists(paramfile)
        if not file_exists:
            print("Warning: SEXTRACTOR paramater file {0} "
                  "does not exist".format(paramfile))
            return None

    # Set default convolution file if not provided
    if convfile == None:
        convfile = os.path.join(LBCgo.__path__[0], 'conf', 'default.conv')
        if verbose:
            print("Using default SExtractor convolution file: {0}".format(convfile))

    # Test for convolution file
    file_exists = os.path.exists(convfile)
    if not file_exists:
        print("Warning: SEXTRACTOR convolution file {0} "
              "does not exist".format(convfile))
        return None

    # Set default neural network weights file if not provided
    if nnwfile == None:
        nnwfile = os.path.join(LBCgo.__path__[0], 'conf', 'default.nnw')
        if verbose:
            print("Using default SExtractor neural network weights file: {0}".format(nnwfile))

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

    # SExtractor flags
    cmd_flags = ' -c '+ configfile + \
        ' -CATALOG_NAME '+outputcatalog + \
        ' -CATALOG_TYPE FITS_LDAC'+ \
        ' -DETECT_THRESH 5.0 -ANALYSIS_THRESH 8.0'+ \
        ' -PARAMETERS_NAME '+paramfile

    # Put together the SExtractor command
    cmd = 'sex '+inputfile+cmd_flags

    try:
        if verbose:
            print('########### SEXTRACTOR run for {0} '
                  '########### \n'.format(inputfile.replace('.cat', '')))
            print(cmd)
            sextract = Popen(shlex.split(cmd),
                             close_fds=True)
        else:
            sextract = Popen(shlex.split(cmd),
                             stdout=DEVNULL,
                             stderr=DEVNULL,
                             close_fds=True)
    except Exception as e:
        print('Oops: source extractor call:', (e))
        return None

    sextract.wait()


def go_scamp(inputfile,
             astrometric_catalog='GAIA-DR3',
             astrometric_method = 'exposure',
             num_iterations = 3,
             configfile=None,
             verbose=True):

    """Run SCAMP to compute an astrometric solution for a chip catalog.

    Matches the SExtractor FITS_LDAC catalog against an astrometric reference
    catalog and writes a ``.head`` WCS solution file. SCAMP is run iteratively
    with progressively tighter tolerance parameters to refine the solution.
    A minimum of 2 iterations is enforced.

    Parameters
    ----------
    inputfile : str
        Path to the chip FITS image or SExtractor catalog (``*.fits`` or
        ``*.cat``). The function converts ``.fits`` to ``.cat`` automatically.
    astrometric_catalog : str, optional
        Reference catalog for cross-matching. Common options: ``'GAIA-DR3'``,
        ``'GAIA-DR2'``, ``'2MASS'``, ``'USNO-B1'``. Default: ``'GAIA-DR3'``
    astrometric_method : str, optional
        Mosaic type strategy passed to SCAMP. Default: ``'exposure'``
    num_iterations : int, optional
        Number of SCAMP iterations. Fewer than 2 will be raised to 2.
        Default: 3
    configfile : str or None, optional
        Path to a SCAMP configuration file. If None, uses the LBCgo default
        ``scamp.lbc.conf``. Default: None
    verbose : bool, optional
        Print SCAMP commands and progress. Default: True

    Raises
    ------
    RuntimeError
        If SCAMP is not found on the system PATH.

    Returns
    -------
    None
        Writes a ``.head`` WCS solution file alongside the input catalog.
    """

    # Check if SCAMP is available
    if not shutil.which('scamp'):
        raise RuntimeError("SCAMP is not available. Please install it.")

    # Make sure the input file is a SEXTRACTOR catalog:
    inputfile = inputfile.replace('.fits','.cat')
    xmlfile = inputfile.replace('.cat','.xml')

    # Make sure we have a configuration file:
    if configfile == None:
        # Use default SCAMP config file from LBCgo package directories
        configfile = os.path.join(LBCgo.__path__[0], 'conf', 'scamp.lbc.conf')
        if verbose:
            print("Using default SCAMP configuration file: {0}".format(configfile))

    # Using only a single iteration of SCAMP doesn't do well enough. Force at
    # least two iterations:
    if num_iterations < 2:
        print('WARNING: Use at least 2 SCAMP iterations. Setting num_iterations = 2...')
        num_iterations = 2

    # Perform the iterations
    for scmpiter in np.arange(num_iterations):
        if scmpiter == 0:
            mosaic_type = 'LOOSE'
            pixscale_maxerr = '1.2'
            position_maxerr = '1'
            posangle_maxerr = '5.0'
            crossid_radius = '7.5'
            aheader_suffix = '.ahead'
        elif scmpiter == 1:
           mosaic_type = 'FIX_FOCALPLANE'
           pixscale_maxerr = '1.1'
           position_maxerr = '0.1'
           posangle_maxerr = '3.0'
           crossid_radius = '5.0'
           aheader_suffix = '.head'
        elif scmpiter == 2:
           mosaic_type = 'FIX_FOCALPLANE'
           pixscale_maxerr = '1.05'
           position_maxerr = '0.05'
           posangle_maxerr = '1.0'
           crossid_radius = '5.0'
           aheader_suffix = '.head'
        else:
           mosaic_type = 'FIX_FOCALPLANE'
           pixscale_maxerr = '1.05'
           position_maxerr = '0.025'
           posangle_maxerr = '1.0'
           crossid_radius = '2.5'
           aheader_suffix = '.head'

        cmd_flags = ' -c '+ configfile + \
            ' -PIXSCALE_MAXERR '+pixscale_maxerr+ \
            ' -POSANGLE_MAXERR '+posangle_maxerr+ \
            ' -POSITION_MAXERR '+position_maxerr+ \
            ' -ASTREF_CATALOG '+astrometric_catalog+ \
            ' -AHEADER_SUFFIX '+aheader_suffix+ \
            ' -CROSSID_RADIUS '+crossid_radius+\
            ' -XML_NAME '+xmlfile
        # ' -MOSAIC_TYPE '+mosaic_type+ \

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
            print('Oops: source Extractor call:', (e))
            return None

        scamp.wait()


    # Read XML file after last iteration
    scamp_diagnostic = (votable.parse(xmlfile)).get_first_table().array
    xy_dispersion = scamp_diagnostic['AstromSigma_Reference'].data
    astrometric_dispersion = np.sqrt(np.sum(xy_dispersion**2))

    # TODO: Do something with the astrometric dispersion


def go_swarp(inputfiles,
             output_filename = None,
             configfile = None,
             verbose = True):
    """Resample and co-add chip images using SWarp.

    Reads SCAMP-produced ``.head`` WCS solutions, reprojects all input chip
    images onto a common astrometric grid, and co-adds them into a single
    mosaic. The output filename is derived from the ``OBJECT`` and ``FILTER``
    headers of the first input file if not specified.

    All input files must share the same filter; a ``ValueError`` is raised
    otherwise.

    Parameters
    ----------
    inputfiles : list of str
        Paths to chip FITS images to co-add. All must have the same filter.
    output_filename : str or None, optional
        Output mosaic filename. If None, derived from the object name and
        filter as ``<object>.<filter>.mos.fits``. Default: None
    configfile : str or None, optional
        Path to a SWarp configuration file. If None, uses the LBCgo default
        ``swarp.lbc.conf``. Default: None
    verbose : bool, optional
        Print the SWarp command and progress. Default: True

    Raises
    ------
    RuntimeError
        If SWarp is not found on the system PATH.
    ValueError
        If the input files span more than one filter.

    Returns
    -------
    None
        Writes ``<output_filename>`` and ``<output_filename>.weight.fits``
        to the current directory. Mean exposure-time-weighted airmass is
        written to the output header.
    """

    # Check if SWarp is available
    if not shutil.which('swarp'):
        raise RuntimeError("SWarp is not available. Please install it.")

    # Make sure we have a configuration file:
    if configfile == None:
        # Use default config file from LBCgo package directories
        configfile = os.path.join(LBCgo.__path__[0], 'conf', 'swarp.lbc.conf')
        if verbose:
            print("Using default SWARP configuration file: {0}".format(configfile))


    # Gather some information about the input files
    keywds = ['object', 'filter', 'exptime', 'imagetyp', 'propid', 'lbcobnam',
                  'airmass', 'HA', 'objra', 'objdec']
    ic_swarp = ImageFileCollection('./', keywords=keywds,
                                    filenames = inputfiles)

    # Check all are the same filter:
    fltrs = ic_swarp.values('filter',unique=True)
    if np.size(fltrs) != 1:
        filters_str = ', '.join(str(f) for f in fltrs)
        raise ValueError(f"All input files must have the same filter for SWARP combination. "
                        f"Found {np.size(fltrs)} different filters: {filters_str}")

    # Calculate mean airmass (weighted by exposure time)
    exp_airmass = np.array(ic_swarp.values('airmass'))
    exp_time = np.array(ic_swarp.values('exptime'))
    airmass = np.average(exp_airmass,weights=exp_time)

    # For now grab the information from the first header:
    if output_filename == None:
        imhead = fits.getheader(inputfiles[0])
        # Shorten the filter names used:
        filter_text = imhead['FILTER']
        filter_text = filter_text.replace('-SLOAN','').replace('-BESSEL','').replace('SDT_Uspec','Uspec')
        # Create final output filename
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

    try:
        if verbose:
            print('########### SWARP image combination '
                  '########### \n')
            print(cmd)
            swarp = Popen(shlex.split(cmd),
                          close_fds=True)
        else:
            swarp = Popen(shlex.split(cmd),
                          stdout=DEVNULL,
                          stderr=DEVNULL,
                          close_fds=True)
    except Exception as e:
        print('Oops: SWARP call:', (e))
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
                lbc_chips = True,
                do_sextractor=True,
                do_scamp=True,
                do_swarp=True,
                astrometric_catalog='GAIA-DR3',
                scamp_iterations = 3,
                verbose=True):
    """Perform astrometric registration and image combination for LBC chip-extracted data.
    
    This function coordinates the complete astrometric processing pipeline for LBC data
    that has been processed through flat fielding and chip extraction. It sequentially
    runs source extraction (SExtractor), astrometric calibration (SCAMP), and image
    combination (SWARP) on individual CCD chip images to produce final registered
    and co-added mosaics.
    
    Processing Steps:
    1. Source extraction on individual chip images using SExtractor
    2. Astrometric solution calculation using SCAMP with iterative refinement
    3. Image resampling and combination using SWARP to create final mosaics
    
    Parameters
    ----------
    filter_directories : str or list of str
        Directory path(s) containing chip-extracted FITS files. Each directory should
        contain individual chip images (e.g., 'object_1.fits', 'object_2.fits', etc.)
        from the chip extraction step.
    lbc_chips : bool or list of int, optional
        CCD chips to process. If True, processes all 4 chips [1,2,3,4].
        Can specify subset as list (e.g., [1,3]). Default: True
    do_sextractor : bool, optional
        Run SExtractor for source detection on each chip image. Creates catalogs
        needed for astrometric calibration. Default: True
    do_scamp : bool, optional
        Run SCAMP for astrometric calibration. Calculates WCS solutions using
        reference catalog cross-matching. Default: True
    do_swarp : bool, optional
        Run SWARP for image resampling and combination. Creates final co-added
        mosaics with corrected astrometry. Default: True
    astrometric_catalog : str, optional
        Reference catalog for astrometric calibration. Common options include
        'GAIA-DR3', 'GAIA-DR2', '2MASS', 'USNO-B1', etc. Default: 'GAIA-DR3'
    scamp_iterations : int, optional
        Number of SCAMP iterations for astrometric solution refinement.
        More iterations improve precision but increase processing time.
        Minimum of 2 recommended. Default: 3
    verbose : bool, optional
        Print detailed processing information and command outputs. Default: True
        
    Returns
    -------
    None
        Function performs file operations and creates output files in the input
        directories. Final products are astrometrically-calibrated combined
        images with '.mos.fits' extension and corresponding weight maps.
        
    Notes
    -----
    - Input directories should contain chip-extracted FITS files from go_extractchips()
    - Requires external tools: SExtractor, SCAMP, and SWARP from astromatic.net
    - Creates intermediate files (.cat, .xml, .head) during processing
    - Final mosaics are named using object name and filter (e.g., 'M31.g.mos.fits')
    - Processing is done per filter directory to maintain filter separation
    - SCAMP uses iterative refinement with progressively tighter tolerances
    
    Examples
    --------
    Basic astrometric processing for all chips:
    >>> go_register(['M31/g-SLOAN/', 'M31/r-SLOAN/'])
    
    Process only specific chips with custom catalog:
    >>> go_register(['NGC4321/V/'], lbc_chips=[1,2], 
    ...              astrometric_catalog='2MASS')
    
    Source extraction and astrometry only (no final combination):
    >>> go_register(['target/filter/'], do_swarp=False)
    
    High-precision astrometry with more iterations:
    >>> go_register(['science/'], scamp_iterations=5)
    """

    # TODO: Add the sextractor, scamp, swarp parameters for input.

    ###### Define which chips to extract if default is chosen:
    if lbc_chips == True:
        lbc_chips = [1,2,3,4]

    # If user enters just a single directory:
    if np.size(filter_directories) == 1 and not isinstance(filter_directories,list):
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
                go_scamp(filename,
                         astrometric_catalog=astrometric_catalog,
                         num_iterations=scamp_iterations,
                         verbose=verbose)

        # Stitch together the images
        # go_swarp = reproject and coadd images
        if do_swarp:
            go_swarp(input_filenames,
                     verbose=verbose)
