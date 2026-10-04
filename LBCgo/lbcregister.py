import os
import numpy as np
import shutil
from subprocess import Popen, DEVNULL
import shlex
from glob import glob
from astropy.io import fits
from ccdproc import  ImageFileCollection
import astropy.io.votable as votable
from astropy.table import Table
import LBCgo
import re
import tempfile
from . import masks as lbcmasks



# Executable names differ between packagings: Debian/Ubuntu ship SExtractor
# as ``source-extractor`` and SWarp as ``SWarp``; Homebrew/conda use ``sex``
# and ``swarp``. The first name found on PATH wins.
ASTROMATIC_NAMES = {'sex': ('sex', 'source-extractor'),
                    'scamp': ('scamp',),
                    'swarp': ('swarp', 'SWarp')}


def find_astromatic_tool(tool):
    """Return the executable name for an astromatic tool, or None if absent.

    Parameters
    ----------
    tool : {'sex', 'scamp', 'swarp'}
        Generic tool name; the alternatives in ``ASTROMATIC_NAMES`` are tried
        in order.
    """
    for name in ASTROMATIC_NAMES[tool]:
        if shutil.which(name):
            return name
    return None


# TODO: Offer an iterative treatment of SCAMP to get desired precision
# TODO: Photometric calibration / selection of filters

def go_sextractor(inputfile,
                configfile=None,
                paramfile = None,
                convfile = None,
                nnwfile = None,
                verbose=True,
                detect_thresh=5.0,
                analysis_thresh=None,
                back_size=32,
                back_filtersize=3,
                deblend_mincont=1e-4,
                use_weight=True,
                use_flags=True,
                weight_file=None,
                flag_file=None,
                subtracted_image=None):
    """Run SExtractor on a single chip image to produce a source catalog.

    Detects sources and writes a FITS_LDAC catalog (``<base>.cat``) alongside
    the input file. Default configuration files are read from the LBCgo
    package ``conf/`` directory. The defaults are tuned for alignment
    catalogs (small background mesh, aggressive deblending so stars on a
    bright galaxy are not merged into its segment).

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
    detect_thresh : float, optional
        ``DETECT_THRESH`` in sigma. Default: 5.0
    analysis_thresh : float or None, optional
        ``ANALYSIS_THRESH`` in sigma. None sets it equal to ``detect_thresh``
        (the config-file value of 1.5 would otherwise apply). Default: None
    back_size : int, optional
        ``BACK_SIZE`` background mesh in pixels (32 px = 7 arcsec). Default: 32
    back_filtersize : int, optional
        ``BACK_FILTERSIZE`` median filter size in meshes. Default: 3
    deblend_mincont : float, optional
        ``DEBLEND_MINCONT``. Stars on a bright galaxy hold well under 0.5 %
        of the galaxy segment's flux, so the SExtractor default (0.005) does
        not separate them. Default: 1e-4
    use_weight : bool, optional
        Pass the chip's ``<base>.weight.fits`` sidecar (written by
        ``go_flatfield``) as a ``MAP_WEIGHT`` image. Silently skipped when no
        weight file exists. Default: True
    use_flags : bool, optional
        Pass the chip's ``<base>.mask.fits`` sidecar as ``FLAG_IMAGE``;
        flags are then available as ``IMAFLAGS_ISO`` if requested in the
        parameter file. Skipped when no mask file exists. Default: True
    weight_file, flag_file : str or None, optional
        Explicit weight / flag images overriding the sidecar names.
    subtracted_image : str or None, optional
        Extended-target mode. A background-subtracted (e.g. ellipse-masked,
        high-passed) version of the chip image with the same pixel grid and
        WCS. SExtractor then runs on it with ``BACK_TYPE MANUAL``,
        ``BACK_VALUE 0`` so that galaxy structure is not re-estimated as
        background. The catalog is still named after ``inputfile``.
        Default: None

    Raises
    ------
    RuntimeError
        If SExtractor (``sex`` or ``source-extractor``) is not found on the
        system PATH.

    Returns
    -------
    None
        Writes a FITS_LDAC catalog to ``<inputfile_base>.cat``.
        Returns None if a required configuration file is missing.
    """
    # Check if SExtractor is available (name depends on the packaging)
    sex_exe = find_astromatic_tool('sex')
    if sex_exe is None:
        raise RuntimeError("SExtractor (sex or source-extractor) is not "
                           "available. Please install it.")

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

    if analysis_thresh is None:
        analysis_thresh = detect_thresh

    # Weight map (inverse variance; 0 = bad pixel) and flag image (mask)
    if weight_file is None and use_weight:
        weight_file = lbcmasks.sidecar_name(inputfile, 'weight')
    if weight_file is not None and not os.path.exists(weight_file):
        if use_weight and verbose:
            print("No weight map found for {0}; running unweighted.".format(inputfile))
        weight_file = None

    if flag_file is None and use_flags:
        flag_file = lbcmasks.sidecar_name(inputfile, 'mask')
    if flag_file is not None and not os.path.exists(flag_file):
        flag_file = None

    # Extended-target mode: background already removed from the input
    detect_image = inputfile if subtracted_image is None else subtracted_image

    # SExtractor's option parser splits on whitespace and does not honour
    # quotes, so a path containing a space (e.g. a Dropbox folder) silently
    # breaks the run. Stage symlinks in a space-free temporary directory.
    paths = {'image': detect_image, 'config': configfile, 'param': paramfile,
             'conv': convfile, 'nnw': nnwfile, 'weight': weight_file,
             'flag': flag_file, 'catalog': outputcatalog}
    needs_links = any(v is not None and ' ' in str(v) for v in paths.values())
    stage = None
    if needs_links or flag_file is not None:
        stage = tempfile.TemporaryDirectory(prefix='lbcgo_sex_')
    if needs_links:
        staged = {}
        for role, path in paths.items():
            if path is None:
                continue
            name = os.path.basename(path).replace(' ', '_')
            staged[role] = os.path.join(stage.name, name)
            if role != 'catalog':
                os.symlink(os.path.abspath(path), staged[role])
        # Side-car header (.head) written by SCAMP, if any
        head = os.path.splitext(detect_image)[0] + '.head'
        if os.path.exists(head):
            os.symlink(os.path.abspath(head),
                       os.path.splitext(staged['image'])[0] + '.head')
        paths = {**paths, **staged}

    # With a flag image, add the per-source mask-flag columns to a staged copy
    # of the parameter file (SExtractor fails if they are requested without a
    # FLAG_IMAGE, so they cannot live in the default file).
    if flag_file is not None:
        with open(paths['param']) as fh:
            active = [ln.split('#')[0].strip() for ln in fh]
        extra = [c for c in ('IMAFLAGS_ISO(1)', 'NIMAFLAGS_ISO(1)')
                 if c not in active]
        if extra:
            augmented = os.path.join(stage.name, 'flags.param')
            with open(paths['param']) as fh, open(augmented, 'w') as out:
                out.write(fh.read().rstrip('\n') + '\n' + '\n'.join(extra) + '\n')
            paths['param'] = augmented

    # SExtractor flags
    cmd_flags = ['-c', paths['config'],
                 '-CATALOG_NAME', paths['catalog'],
                 '-CATALOG_TYPE', 'FITS_LDAC',
                 '-PARAMETERS_NAME', paths['param'],
                 '-FILTER_NAME', paths['conv'],
                 '-STARNNW_NAME', paths['nnw'],
                 '-DETECT_THRESH', str(detect_thresh),
                 '-ANALYSIS_THRESH', str(analysis_thresh),
                 '-BACK_SIZE', str(back_size),
                 '-BACK_FILTERSIZE', str(back_filtersize),
                 '-DEBLEND_MINCONT', str(deblend_mincont)]
    if weight_file is not None:
        cmd_flags += ['-WEIGHT_TYPE', 'MAP_WEIGHT',
                      '-WEIGHT_IMAGE', paths['weight']]
    if flag_file is not None:
        cmd_flags += ['-FLAG_IMAGE', paths['flag']]
    if subtracted_image is not None:
        cmd_flags += ['-BACK_TYPE', 'MANUAL', '-BACK_VALUE', '0']

    # Put together the SExtractor command
    cmd = [sex_exe, paths['image']] + cmd_flags

    try:
        if verbose:
            print('########### SEXTRACTOR run for {0} '
                  '########### \n'.format(inputfile.replace('.cat', '')))
            print(shlex.join(cmd))
            sextract = Popen(cmd,
                             close_fds=True)
        else:
            sextract = Popen(cmd,
                             stdout=DEVNULL,
                             stderr=DEVNULL,
                             close_fds=True)
        sextract.wait()
        # Bring the catalog out of the staging directory
        if stage is not None and os.path.exists(paths['catalog']):
            shutil.move(paths['catalog'], outputcatalog)
    except Exception as e:
        print('Oops: source extractor call:', (e))
        return None
    finally:
        if stage is not None:
            stage.cleanup()


def scamp_iteration_settings(iteration):
    """Return the SCAMP tolerances for a 0-based iteration number.

    Iteration 0 is loose (``MOSAIC_TYPE LOOSE``, reads ``.ahead`` files);
    later iterations fix the focal plane, read the previous ``.head`` as
    input astrometry, and tighten the match tolerances.
    """
    if iteration == 0:
        return dict(mosaic_type='LOOSE', pixscale_maxerr='1.2',
                    position_maxerr='1', posangle_maxerr='5.0',
                    crossid_radius='7.5', aheader_suffix='.ahead')
    if iteration == 1:
        return dict(mosaic_type='FIX_FOCALPLANE', pixscale_maxerr='1.1',
                    position_maxerr='0.1', posangle_maxerr='3.0',
                    crossid_radius='5.0', aheader_suffix='.head')
    if iteration == 2:
        return dict(mosaic_type='FIX_FOCALPLANE', pixscale_maxerr='1.05',
                    position_maxerr='0.05', posangle_maxerr='1.0',
                    crossid_radius='5.0', aheader_suffix='.head')
    return dict(mosaic_type='FIX_FOCALPLANE', pixscale_maxerr='1.05',
                position_maxerr='0.025', posangle_maxerr='1.0',
                crossid_radius='2.5', aheader_suffix='.head')


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
        Deprecated and ignored (it never had an effect). The mosaic type is
        set per iteration; for focal-plane solutions use :func:`go_scamp_joint`.
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
    if find_astromatic_tool('scamp') is None:
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
        it = scamp_iteration_settings(scmpiter)

        cmd_flags = ' -c '+ configfile + \
            ' -MOSAIC_TYPE '+it['mosaic_type']+ \
            ' -PIXSCALE_MAXERR '+it['pixscale_maxerr']+ \
            ' -POSANGLE_MAXERR '+it['posangle_maxerr']+ \
            ' -POSITION_MAXERR '+it['position_maxerr']+ \
            ' -ASTREF_CATALOG '+astrometric_catalog+ \
            ' -AHEADER_SUFFIX '+it['aheader_suffix']+ \
            ' -CROSSID_RADIUS '+it['crossid_radius']+\
            ' -XML_NAME '+xmlfile

        # Create the final command:
        cmd = find_astromatic_tool('scamp')+' '+inputfile+cmd_flags

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


CHIP_FILE_RE = re.compile(r'^(?P<base>.*)_(?P<chip>\d+)\.fits$')


def group_chips_by_exposure(chip_files):
    """Group ``<base>_<chip>.fits`` names by exposure.

    Returns
    -------
    dict
        ``{base: [chip files sorted by chip number]}`` in order of first
        appearance; ``base`` includes the directory.
    """
    groups = {}
    for fl in chip_files:
        m = CHIP_FILE_RE.match(fl)
        if m is None:
            raise ValueError("Not a chip file (<base>_<chip>.fits): "
                             "{0}".format(fl))
        groups.setdefault(m.group('base'), []).append((int(m.group('chip')), fl))
    return {b: [f for _, f in sorted(v)] for b, v in groups.items()}


def merge_ldac(chip_catalogs, output_catalog):
    """Concatenate the LDAC extensions of several catalogs into one.

    SCAMP's focal-plane modes need one catalog per exposure with one
    (``LDAC_IMHEAD``, ``LDAC_OBJECTS``) HDU pair per chip. The pairs are
    appended in the order the catalogs are given (chip order); the primary
    HDU comes from the first catalog. ``FITSFILE``/``FITSEXT``/``FITSNEXT``
    in each ``LDAC_IMHEAD`` are renumbered to describe the merged file.

    Parameters
    ----------
    chip_catalogs : list of str
        SExtractor FITS_LDAC catalogs, one per chip.
    output_catalog : str
        Merged catalog to write (overwritten).

    Returns
    -------
    str
        ``output_catalog``
    """
    # Work on the raw FITS bytes: an astropy round trip rewrites the
    # space-padded 80-character header cards in LDAC_IMHEAD as NUL-padded,
    # which makes SCAMP fault (SIGSEGV/SIGBUS).
    blocks, n_imhead = [], 0
    for icat, cat in enumerate(chip_catalogs):
        with open(cat, 'rb') as fh:
            raw = fh.read()
        with fits.open(cat, memmap=False) as hdul:
            # The primary comes from the first catalog only; the rest are
            # (IMHEAD, OBJECTS) pairs
            for i in range(0 if icat == 0 else 1, len(hdul)):
                info = hdul.fileinfo(i)
                block = bytearray(raw[info['hdrLoc']:info['datLoc'] + info['datSpan']])
                is_imhead = hdul[i].name == 'LDAC_IMHEAD'
                n_imhead += is_imhead
                blocks.append((is_imhead, block, info['datLoc'] - info['hdrLoc'],
                               hdul[i].header['NAXIS1'] if is_imhead else 0))

    # SCAMP locates each extension through FITSEXT/FITSNEXT (and names the
    # field via FITSFILE) in its LDAC_IMHEAD. Chip catalogs all say 1 of 1;
    # renumber as SExtractor does for a native multi-extension image.
    exposure_name = os.path.basename(output_catalog).replace('.cat', '.fits')
    ext = 0
    with open(output_catalog, 'wb') as out:
        for is_imhead, block, dstart, nbytes in blocks:
            if is_imhead:
                ext += 1
                new_values = {'FITSFILE': exposure_name, 'FITSEXT': ext,
                              'FITSNEXT': n_imhead}
                for off in range(dstart, dstart + nbytes, 80):
                    key = block[off:off + 8].decode('ascii').strip()
                    if key in new_values:
                        old = fits.Card.fromstring(block[off:off + 80].decode('ascii'))
                        block[off:off + 80] = fits.Card(
                            key, new_values[key], old.comment).image.encode('ascii')
            out.write(block)
    return output_catalog


def split_head(exposure_head, chip_heads):
    """Split a multi-section SCAMP ``.head`` file into per-chip files.

    SCAMP writes one section per catalog extension, each terminated by an
    ``END`` card, in extension (chip) order.

    Parameters
    ----------
    exposure_head : str
        The ``.head`` file SCAMP wrote for a merged exposure catalog.
    chip_heads : list of str
        Output filenames, one per chip, in the same order as the merged
        catalog's extensions.

    Raises
    ------
    ValueError
        If the number of sections differs from ``len(chip_heads)``.
    """
    with open(exposure_head) as fh:
        lines = fh.read().splitlines()
    sections, current = [], []
    for ln in lines:
        current.append(ln)
        if ln.strip() == 'END':
            sections.append(current)
            current = []
    if len(sections) != len(chip_heads):
        raise ValueError("{0} has {1} header sections but {2} chips were "
                         "expected".format(exposure_head, len(sections),
                                           len(chip_heads)))
    for name, sect in zip(chip_heads, sections):
        with open(name, 'w') as out:
            out.write('\n'.join(sect) + '\n')


def read_head_sections(head_file):
    """Return one ``astropy.io.fits.Header`` per ``END``-terminated section."""
    with open(head_file) as fh:
        lines = fh.read().splitlines()
    headers, current = [], []
    for ln in lines:
        current.append(ln)
        if ln.strip() == 'END':
            headers.append(fits.Header.fromstring('\n'.join(current), sep='\n'))
            current = []
    return headers


def scamp_qa_table(xmlfile, groups, max_ref_rms=0.2, min_xy_contrast=2.0,
                   astref_catalog='GAIA-DR3', astref_epoch=None):
    """Build the astrometry QA table from a joint SCAMP run.

    Per-chip numbers come from the ``.head`` sections (``ASTIRMS``/``ASTRRMS``
    are in degrees, converted to arcsec; ``FLXSCALE``), per-exposure numbers
    from the SCAMP XML ``Fields`` table (``XY_Contrast``, reference-match
    count and rms).

    Parameters
    ----------
    xmlfile : str
        SCAMP XML written by the last iteration.
    groups : dict
        ``{exposure base: [chip files]}`` as from
        :func:`group_chips_by_exposure`; ``<base>_exp.head`` must exist.
    max_ref_rms : float, optional
        A chip is flagged when either axis' reference rms (arcsec) exceeds
        this, or is exactly 0 (no reference matches).
    min_xy_contrast : float, optional
        An exposure is flagged when SCAMP's ``XY_Contrast`` is below this.
    astref_catalog, astref_epoch : str, optional
        Recorded in the table metadata with the SCAMP version.

    Returns
    -------
    astropy.table.Table
        One row per exposure chip; ``bad`` is True where any check failed and
        ``reason`` lists the failed checks.
    """
    fields = votable.parse(xmlfile).get_first_table().to_table()
    by_cat = {os.path.basename(str(r['Catalog_Name'])): r for r in fields}

    rows, scamp_version = [], None
    for base, chips in groups.items():
        heads = read_head_sections(base + '_exp.head')
        if len(heads) != len(chips):
            raise ValueError("{0}_exp.head has {1} sections for {2} chips"
                             .format(base, len(heads), len(chips)))
        field = by_cat.get(os.path.basename(base) + '_exp.cat')
        if scamp_version is None:
            for h in heads[:1]:
                m = re.search(r'SCAMP version (\S+)', ' '.join(map(str, h['HISTORY'])))
                scamp_version = m.group(1) if m else None
        contrast = float(field['XY_Contrast']) if field is not None else np.nan
        ref_ndeg = int(field['NDeg_Reference']) if field is not None else 0
        for ext, (chipfile, h) in enumerate(zip(chips, heads), start=1):
            m = CHIP_FILE_RE.match(chipfile)
            irms = [3600.*float(h.get('ASTIRMS{0}'.format(i), 0.)) for i in (1, 2)]
            rrms = [3600.*float(h.get('ASTRRMS{0}'.format(i), 0.)) for i in (1, 2)]
            reasons = []
            if max(rrms) == 0.:
                reasons.append('no_ref_match')
            elif max(rrms) > max_ref_rms:
                reasons.append('ref_rms')
            if not contrast >= min_xy_contrast:
                reasons.append('low_contrast')
            rows.append((os.path.basename(base), int(m.group('chip')), ext,
                         irms[0], irms[1], rrms[0], rrms[1],
                         float(h.get('FLXSCALE', np.nan)), ref_ndeg, contrast,
                         bool(reasons), ','.join(reasons)))
    qa = Table(rows=rows, names=['exposure', 'chip', 'ext',
                                 'int_rms_x', 'int_rms_y',
                                 'ref_rms_x', 'ref_rms_y', 'flxscale',
                                 'n_ref_dof', 'xy_contrast', 'bad', 'reason'],
               dtype=[str, int, int, float, float, float, float, float, int,
                      float, bool, str])
    for col in ('int_rms_x', 'int_rms_y', 'ref_rms_x', 'ref_rms_y'):
        qa[col].unit = 'arcsec'
    qa.meta.update(scamp_version=scamp_version, astref_catalog=astref_catalog,
                   astref_epoch=astref_epoch, max_ref_rms=max_ref_rms,
                   min_xy_contrast=min_xy_contrast)
    return qa


def go_scamp_joint(chip_files,
                   astrometric_catalog='GAIA-DR3',
                   num_iterations=3,
                   configfile=None,
                   astref_epoch='FIELDS_AVERAGE',
                   qa_file='astrometry_qa.ecsv',
                   max_ref_rms=0.2,
                   min_xy_contrast=2.0,
                   verbose=True):
    """Solve the astrometry of all exposures in a filter directory jointly.

    Merges the four chip catalogs (``<base>_<chip>.cat``) of each exposure
    into ``<base>_exp.cat``, runs SCAMP once over all exposure catalogs
    (``STABILITY_TYPE INSTRUMENT``, ``ASTRINSTRU_KEY FILTER``: one
    distortion/focal-plane solution per filter directory), then splits the
    resulting ``.head`` into ``<base>_<chip>.head`` for SWarp and writes a QA
    table.

    Parameters
    ----------
    chip_files : list of str
        Chip images ``<base>_<chip>.fits`` (their ``.cat`` must exist) from
        one filter directory.
    astrometric_catalog : str, optional
        SCAMP ``ASTREF_CATALOG``. Default: ``'GAIA-DR3'``
    num_iterations : int, optional
        Iterations with tightening tolerances (minimum 2). Default: 3
    configfile : str or None, optional
        SCAMP configuration; None uses ``scamp.lbc.conf``.
    astref_epoch : {'FIELDS_AVERAGE', 'ORIGINAL', 'MANUAL'}, optional
        SCAMP ``ASTREFEPOCH_TYPE``: epoch to which reference-catalog proper
        motions are propagated. Default: ``'FIELDS_AVERAGE'``
    qa_file : str or None, optional
        QA table (ECSV) written into the filter directory; None to skip.
    max_ref_rms : float, optional
        Flag chips whose reference-catalog astrometric rms (arcsec, either
        axis) exceeds this. Default: 0.2 (not yet tuned on real data)
    min_xy_contrast : float, optional
        Flag exposures whose SCAMP ``XY_Contrast`` falls below this.
        Default: 2.0 (not yet tuned on real data)
    verbose : bool, optional

    Returns
    -------
    astropy.table.Table
        The QA table (one row per exposure chip), or None if SCAMP failed.
    """
    scamp_exe = find_astromatic_tool('scamp')
    if scamp_exe is None:
        raise RuntimeError("SCAMP is not available. Please install it.")
    if num_iterations < 2:
        print('WARNING: Use at least 2 SCAMP iterations. Setting num_iterations = 2...')
        num_iterations = 2

    if configfile is None:
        configfile = os.path.join(LBCgo.__path__[0], 'conf', 'scamp.lbc.conf')
    groups = group_chips_by_exposure(chip_files)
    workdir = os.path.dirname(os.path.abspath(next(iter(groups.values()))[0]))

    # One merged catalog per exposure
    exp_cats = {}
    for base, chips in groups.items():
        chip_cats = [c.replace('.fits', '.cat') for c in chips]
        missing = [c for c in chip_cats if not os.path.exists(c)]
        if missing:
            raise FileNotFoundError("Missing SExtractor catalogs: "
                                    "{0}".format(missing))
        exp_cats[base] = merge_ldac(chip_cats, base + '_exp.cat')
    if len({len(c) for c in groups.values()}) > 1:
        print("WARNING: exposures have different numbers of chips; SCAMP's "
              "FIX_FOCALPLANE assumes a common layout.")

    # SCAMP's option parser splits on spaces: run in the catalog directory
    # with relative names, and use a space-free copy of the config if needed.
    stage = None
    if ' ' in configfile:
        stage = tempfile.TemporaryDirectory(prefix='lbcgo_scamp_')
        staged = os.path.join(stage.name, os.path.basename(configfile))
        shutil.copy(configfile, staged)
        configfile = staged
    xmlname = 'scamp.xml'
    rel_cats = [os.path.relpath(os.path.abspath(c), workdir)
                for c in exp_cats.values()]
    try:
        for scmpiter in range(num_iterations):
            it = scamp_iteration_settings(scmpiter)
            cmd = ([scamp_exe] + rel_cats +
                   ['-c', configfile,
                    '-MOSAIC_TYPE', it['mosaic_type'],
                    '-PIXSCALE_MAXERR', it['pixscale_maxerr'],
                    '-POSANGLE_MAXERR', it['posangle_maxerr'],
                    '-POSITION_MAXERR', it['position_maxerr'],
                    '-CROSSID_RADIUS', it['crossid_radius'],
                    '-AHEADER_SUFFIX', it['aheader_suffix'],
                    '-ASTREF_CATALOG', astrometric_catalog,
                    '-ASTREFEPOCH_TYPE', astref_epoch,
                    '-STABILITY_TYPE', 'INSTRUMENT',
                    '-ASTRINSTRU_KEY', 'FILTER',
                    '-MERGEDOUTCAT_TYPE', 'NONE',
                    '-XML_NAME', xmlname])
            if verbose:
                print('########### SCAMP joint iteration {0} for {1} '
                      'exposures ###########'.format(scmpiter + 1, len(rel_cats)))
                print(shlex.join(cmd))
                proc = Popen(cmd, cwd=workdir, close_fds=True)
            else:
                proc = Popen(cmd, cwd=workdir, stdout=DEVNULL,
                             stderr=DEVNULL, close_fds=True)
            if proc.wait() != 0:
                raise RuntimeError("SCAMP exited with status {0} (iteration "
                                   "{1})".format(proc.returncode, scmpiter + 1))
    except RuntimeError:
        raise
    except Exception as e:
        print('Oops: SCAMP call:', (e))
        return None
    finally:
        if stage is not None:
            stage.cleanup()

    # Per-chip heads for SWarp
    for base, chips in groups.items():
        split_head(base + '_exp.head', [c.replace('.fits', '.head') for c in chips])

    qa = scamp_qa_table(os.path.join(workdir, xmlname), groups,
                        max_ref_rms=max_ref_rms,
                        min_xy_contrast=min_xy_contrast,
                        astref_catalog=astrometric_catalog,
                        astref_epoch=astref_epoch)
    nbad = int(np.sum(qa['bad']))
    if nbad:
        print("WARNING: {0} chip astrometric solution(s) flagged; see "
              "{1}".format(nbad, qa_file))
    if qa_file is not None:
        qa.write(os.path.join(workdir, qa_file), overwrite=True)
    return qa


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
    swarp_exe = find_astromatic_tool('swarp')
    if swarp_exe is None:
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
    cmd = swarp_exe + ' ' + inputfile_text + cmd_flags

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
                sextractor_args = None,
                scamp_joint = True,
                scamp_args = None,
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
    sextractor_args : dict or None, optional
        Extra keyword arguments passed to :func:`go_sextractor` (e.g.
        ``dict(detect_thresh=3, back_size=64, use_weight=False)``).
        Default: None
    scamp_joint : bool, optional
        Solve all exposures of a filter directory with one SCAMP run
        (:func:`go_scamp_joint`; focal-plane and distortion shared across
        exposures, writes ``astrometry_qa.ecsv``). False restores the legacy
        independent per-chip solutions (:func:`go_scamp`). Default: True
    scamp_args : dict or None, optional
        Extra keyword arguments for :func:`go_scamp_joint` (e.g.
        ``dict(max_ref_rms=0.1)``). Default: None
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

        # go_sextractor = find sources (per chip)
        if do_sextractor:
            for filename in input_filenames:
                go_sextractor(filename, verbose=verbose,
                              **(sextractor_args or {}))

        # go_scamp_joint = one solution for all exposures in this directory;
        # go_scamp = legacy independent solution for each chip
        if do_scamp and scamp_joint and input_filenames:
            go_scamp_joint(input_filenames,
                           astrometric_catalog=astrometric_catalog,
                           num_iterations=scamp_iterations,
                           verbose=verbose,
                           **(scamp_args or {}))
        elif do_scamp:
            for filename in input_filenames:
                go_scamp(filename,
                         astrometric_catalog=astrometric_catalog,
                         num_iterations=scamp_iterations,
                         verbose=verbose)

        # Stitch together the images
        # go_swarp = reproject and coadd images
        if do_swarp:
            go_swarp(input_filenames,
                     verbose=verbose)
