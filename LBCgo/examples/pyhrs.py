# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the base processing for HRS data

# See

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import numpy as np
from astropy.extern import six

from astropy import units as u
import scipy.ndimage as nd

import ccdproc

from .hrsorder import HRSOrder

__all__ = ['ccd_process', 'create_masterbias', 'hrs_process', 'blue_process',
           'red_process', 'flatfield_science_order', 'flatfield_science']


def ccd_process(ccd, oscan=None, trim=None, error=False, masterbias=None,
                bad_pixel_mask=None, gain=None, rdnoise=None,
                oscan_median=True, oscan_model=None):
    """Perform basic processing on ccd data.

       The following steps can be included:
        * overscan correction
        * trimming of the image
        * create edeviation frame
        * gain correction
        * add a mask to the data
        * subtraction of master bias

       The task returns a processed `ccdproc.CCDData` object.

    Parameters
    ----------
    ccd: `ccdproc.CCDData`
        Frame to be reduced

    oscan: None, str, or, `~ccdproc.ccddata.CCDData`
        For no overscan correction, set to None.   Otherwise proivde a region
        of `ccd` from which the overscan is extracted, using the FITS
        conventions for index order and index start, or a
        slice from `ccd` that contains the overscan.

    trim: None or str
        For no trim correction, set to None.   Otherwise proivde a region
        of `ccd` from which the image should be trimmed, using the FITS
        conventions for index order and index start.

    error: boolean
        If True, create an uncertainty array for ccd

    masterbias: None, `~numpy.ndarray`,  or `~ccdproc.CCDData`
        A materbias frame to be subtracted from ccd.

    bad_pixel_mask: None or `~numpy.ndarray`
        A bad pixel mask for the data. The bad pixel mask should be in given
        such that bad pixels havea value of 1 and good pixels a value of 0.

    gain: None or `~astropy.Quantity`
        Gain value to multiple the image by to convert to electrons

    rdnoise: None or `~astropy.Quantity`
        Read noise for the observations.  The read noise should be in
        `~astropy.units.electron`


    oscan_median :  bool, optional
        If true, takes the median of each line.  Otherwise, uses the mean

    oscan_model :  `~astropy.modeling.Model`, optional
        Model to fit to the data.  If None, returns the values calculated
        by the median or the mean.

    Returns
    -------
    ccd: `ccdproc.CCDData`
        Reduded ccd

    Examples
    --------

    1. To overscan, trim, and gain correct a data set:

    >>> import numpy as np
    >>> from astropy import units as u
    >>> from hrsprocess import ccd_process
    >>> ccd = CCDData(np.ones([100, 100]), unit=u.adu)
    >>> nccd = ccd_process(ccd, oscan='[1:10,1:100]', trim='[10:100, 1,100]',
                           error=False, gain=2.0*u.electron/u.adu)


    """
    # make a copy of the object
    nccd = ccd.copy()

    # apply the overscan correction
    if isinstance(oscan, ccdproc.CCDData):
        nccd = ccdproc.subtract_overscan(nccd, overscan=oscan,
                                         median=oscan_median,
                                         model=oscan_model)
    elif isinstance(oscan, six.string_types):
        nccd = ccdproc.subtract_overscan(nccd, fits_section=oscan,
                                         median=oscan_median,
                                         model=oscan_model)
    elif oscan is None:
        pass
    else:
        raise TypeError('oscan is not None, a string, or CCDData object')

    # apply the trim correction
    if isinstance(trim, six.string_types):
        nccd = ccdproc.trim_image(nccd, fits_section=trim)
    elif trim is None:
        pass
    else:
        raise TypeError('trim is not None or a string')

    # test subtracting the master bias
    if isinstance(masterbias, ccdproc.CCDData):
        nccd = ccdproc.subtract_bias(nccd, masterbias)
    elif isinstance(masterbias, np.ndarray):
        nccd.data = nccd.data - masterbias
    elif masterbias is None:
        pass
    else:
        raise TypeError(
            'masterbias is not None, numpy.ndarray,  or a CCDData object')

    # create the error frame
    if error and gain is not None and rdnoise is not None:
        nccd = ccdproc.create_deviation(nccd, gain=gain, readnoise=rdnoise)
    elif error and (gain is None or rdnoise is None):
        raise ValueError(
            'gain and rdnoise must be specified to create error frame')

    # apply the bad pixel mask
    if isinstance(bad_pixel_mask, np.ndarray):
        nccd.mask = bad_pixel_mask
    elif bad_pixel_mask is None:
        pass
    else:
        raise TypeError('bad_pixel_mask is not None or numpy.ndarray')

    # apply the gain correction
    if isinstance(gain, u.quantity.Quantity):
        nccd = ccdproc.gain_correct(nccd, gain)
    elif gain is None:
        pass
    else:
        raise TypeError('gain is not None or astropy.Quantity')

    return nccd


def hrs_process(image_name, ampsec=[], oscansec=[], trimsec=[],
                masterbias=None, error=False, bad_pixel_mask=None, flip=False,
                rdnoise=None, oscan_median=True, oscan_model=None):
    """Processing required for HRS observations.  If the images have multiple
       amps, then this will process each part of the image and recombine them
       into for the final results

    Parameters
    ----------
    image_name: str
       Name of file to be processed

    ampsec: list
       List of ampsections.  This list should have the same length as the
       number of amps in the data set.  The sections should be given
       in the format of fits_sections (see below).

    oscansec: list
       List of overscan sections.  This list should have the same length as the
       number of amps in the data set.  The sections should be given
       in the format of fits_sections (see below).

    trimsec: list
       List of overscan sections.  This list should have the same length as the
       number of amps in the data set.  The sections should be given
       in the format of fits_sections (see below).

    error: boolean
        If True, create an uncertainty array for ccd

    masterbias: None, `~numpy.ndarray`,  or `~ccdproc.CCDData`
        A materbias frame to be subtracted from ccd.

    bad_pixel_mask: None or `~numpy.ndarray`
        A bad pixel mask for the data. The bad pixel mask should be in given
        such that bad pixels havea value of 1 and good pixels a value of 0.


    flip: boolean
        If True, the image will be flipped such that the orders run from the
        bottom of the image to the top and the dispersion runs from the left
        to the right.

    rdnoise: None or `~astropy.Quantity`
        Read noise for the observations.  The read noise should be in
        `~astropy.units.electron`

    oscan_median :  bool, optional
        If true, takes the median of each line.  Otherwise, uses the mean

    oscan_model :  `~astropy.modeling.Model`, optional
        Model to fit to the data.  If None, returns the values calculated
        by the median or the mean.

    Returns
     -------
    ccd: `~ccdproc.CCDData`
        Data processed and


    Notes
    -----

    The format of the `fits_section` string follow the rules for slices that
    are consistent with the FITS standard (v3) and IRAF usage of keywords like
    TRIMSEC and BIASSEC. Its indexes are one-based, instead of the
    python-standard zero-based, and the first index is the one that increases
    most rapidly as you move through the array in memory order, opposite the
    python ordering.

    The 'fits_section' argument is provided as a convenience for those who are
    processing files that contain TRIMSEC and BIASSEC. The preferred, more
    pythonic, way of specifying the overscan is to do it by indexing the data
    array directly with the `overscan` argument.

    """
    # read in the data
    ccd = ccdproc.CCDData.read(image_name, unit=u.adu)

    try:
        namps = ccd.header['CCDAMPS']
    except KeyError:
        namps = ccd.header['CCDNAMPS']

    # thow errors for the wrong number of amps
    if len(ampsec) != namps:
        raise ValueError('Number of ampsec does not equal number of amps')
    if len(oscansec) != namps:
        raise ValueError('Number of oscansec does not equal number of amps')
    if len(trimsec) != namps:
        raise ValueError('Number of trimsec does not equal number of amps')

    if namps == 1:
        if ccd.header['OBSTYPE'] == 'Bias':
            gain = None
        else:
            gain = float(ccd.header['gain'].split()[0]) * u.electron / u.adu

        nccd = ccd_process(ccd, oscan=oscansec[0], trim=trimsec[0],
                           error=error, masterbias=masterbias,
                           bad_pixel_mask=bad_pixel_mask, gain=gain,
                           rdnoise=rdnoise, oscan_median=oscan_median,
                           oscan_model=oscan_model)
    else:
        ccd_list = []
        xsize = 0
        ysize = 0
        gain_ave = 0

        # determine size of the image
        ys, xs = ccd.data.shape
        x1 = 0
        x2 = 0
        y1 = 0
        y2 = ys
        for i in range(namps):
            amp = ccdproc.utils.slices.slice_from_string(ampsec[i], fits_convention=True)
            ay1, ay2, _ = amp[0].indices(ys)
            ax1, ax2, _ = amp[1].indices(xs)
            trim = ccdproc.utils.slices.slice_from_string(trimsec[i], fits_convention=True)
            ty1, ty2, _ = trim[0].indices(ys)
            tx1, tx2, _ = trim[1].indices(xs)

            dx1 = max(tx1, ax1)
            x1 = min(x1 or dx1, dx1)

            # assumes bias frames are given in increasing x-order
            x2 = ax2 - (ax2 - ax1 - tx2)

        data = np.zeros((y2 - y1, x2 - x1))
        if ccd.mask is not None:
            mask = np.zeros(((y2 - y1, x2 - x1)))
        else:
            mask = None

        if ccd.uncertainty is not None:
            raise NotImplementedError(
                'Support for uncertainties not implimented yet')
        else:
            uncertainty = None

        for i in range(namps):
            gain_ave += float(ccd.header['gain'].split()[i])
        gain_ave = gain_ave / namps * u.electron / u.adu
        for i in range(namps):
            # trim the image to just the amplifier
            cc = ccdproc.trim_image(ccd, fits_section=ampsec[i])

            # determine its position in the area after correction
            amp = ccdproc.utils.slices.slice_from_string(ampsec[i], fits_convention=True)
            ay1, ay2, _ = amp[0].indices(ys)
            ax1, ax2, _ = amp[1].indices(xs)
            trim = ccdproc.utils.slices.slice_from_string(trimsec[i], fits_convention=True)
            ty1, ty2, _ = trim[0].indices(ys)
            tx1, tx2, _ = trim[1].indices(xs)
            dx1 = max(ax1 - x1, 0)
            dx2 = ax2 - x1
            dy1 = ay1
            dy2 = ay2

            # gain correct the arrray
            gain = float(ccd.header['gain'].split()[i]) * u.electron / u.adu
            ncc = ccd_process(cc, oscan=oscansec[i], trim=trimsec[i],
                              error=False, masterbias=None, gain=gain,
                              bad_pixel_mask=None, rdnoise=rdnoise,
                              oscan_median=oscan_median,
                              oscan_model=oscan_model)

            # normalize all to have the same gain
            ncc = ncc.divide(gain_ave)

            data[dy1:dy2, dx1:dx2] = ncc.data
            if ccd.mask is not None:
                mask[dy1:dy2, dx1:dx2] = ncc.mask

        # create new entryy
        nccd = ccdproc.CCDData(data, unit=ncc.unit, mask=mask,
                               uncertainty=uncertainty)
        nccd.header = ccd.header
        nccd = ccd_process(nccd, masterbias=masterbias, error=error, gain=gain_ave,
                           rdnoise=rdnoise, bad_pixel_mask=bad_pixel_mask)

    if flip:
        nccd.data = nccd.data[::-1, ::-1]
        if (nccd.mask is not None):
            nccd.mask = nccd.mask[::-1, ::-1]
        if (nccd.uncertainty is not None):
            nccd.uncertainty = nccd.uncertainty[::-1, ::-1]

    return nccd


def blue_process(infile, masterbias=None, error=False, rdnoise=None, oscan_correct=False):
    """Process a blue frame
    """
    # check to make sure it is a blue file
    ccd = ccdproc.CCDData.read(infile, unit=u.adu)
    try:
        namps = ccd.header['CCDAMPS']
    except KeyError:
        namps = ccd.header['CCDNAMPS']

    # reduce file
    try:
        blueamp = [ccd.header['AMPSEC'].strip()]
        if oscan_correct:
            bluescan = [ccd.header['BIASSEC'].strip()]
        else:
            bluescan = [None]
        bluetrim = [ccd.header['DATASEC'].strip()]
        # ugly hack for when two amps
        if namps > 1: raise Exception()
    except:
        blueamp = ['[1:1050,:]', '[1051:2100,:]']
        if oscan_correct:
            bluescan = ['[1:26,:]', '[1025:1050,:]']
        else:
            bluescan = [None, None]
        bluetrim = ['[27:1050,:]', '[1:1024,:]']

    flip = True
    ccd = hrs_process(infile, ampsec=blueamp, oscansec=bluescan,
                      trimsec=bluetrim, masterbias=masterbias, error=error,
                      rdnoise=rdnoise, flip=flip)
    # this is in place to deal with changes from one amp to two
    if namps == 1:
        ccd.data = ccd.data[:, ::-1]
        if (ccd.mask is not None):
            ccd.mask = ccd.mask[:, ::-1]
        if (ccd.uncertainty is not None):
            ccd.uncertainty = ccd.uncertainty[:, ::-1]

    return ccd


def red_process(infile, masterbias=None, error=None, rdnoise=None, oscan_correct=False):
    """Process a blue frame
    """

    ccd = ccdproc.CCDData.read(infile, unit=u.adu)
    try:
        namps = ccd.header['CCDAMPS']
    except KeyError:
        namps = ccd.header['CCDNAMPS']

    if namps == 1:
        redamp = ['[1:4122,1:4112]']
        if oscan_correct:
            redscan = ['[1:25,1:4112]']
        else:
            redscan = [None]
        redtrim = ['[27:4122,1:4112]']
    else:
        redamp = ['[1:2074,1:2056]', '[1:2074,2057:4112]', '[2075:4148,1:2056]', '[2075:4122,2057:4112]']
        redscan = namps * [None]
        redtrim = ['[27:2074,:]', '[27:2074,:]', '[1:2048,:]', '[1:2048,:]']

    ccd = hrs_process(infile, ampsec=redamp, oscansec=redscan,
                      trimsec=redtrim, masterbias=masterbias, error=error,
                      rdnoise=rdnoise, flip=False)
    return ccd


def create_masterbias(image_list):
    """Create a master bias frame from a list of images

    Parameters
    ----------
    image_list: list
        List contain the file names to be processed

    Returns
    -------
    masterbias: ccddata.CCDData
        Combine master bias from the biases supplied in image_list

    """
    # determine whether they are red or blue
    if os.path.basename(image_list[0]).startswith('H'):
        func = blue_process
    elif os.path.basename(image_list[0]).startswith('R'):
        func = red_process
    else:
        raise TypeError('These are not standard HRS frames')

    # reduce the files
    ccd_list = []
    for image_name in image_list:
        ccd = func(image_name, masterbias=None, error=False)
        ccd_list.append(ccd)

    # combine the files
    cb = ccdproc.Combiner(ccd_list)
    nccd = cb.median_combine(median_func=np.median)

    return nccd


def flatfield_science_order(hrs, flat_hrs, median_filter_size=None, interp=False):
    """Apply a flat field to an hrs order.  If median_filter_size is set,
       the process first removes the overall flat shape by dividing out
       a median filter of the data.   Then it is divided through the
       science data.  This is done on a row by row basis after removing
       the shape of the science data.

    Parameters
    ----------

    hrs: `~pyhrs.HRSOrder`
        Object describing a single HRS order

    flat_hrs: `~pyhrs.HRSOrder`
        Object describing the flatfield for a single HRS order

    median_filter_size: None or int
        Size for median filter to be run on the data and remove the general
        flat field shape

    Returns
    -------
    hrs: `~pyhrs.HRSOrder`
        An hrs object with the flatfield applied to the flux property.

    """

    # set through each line and correct for the shape of the flux
    # if needed
    if median_filter_size is not None:
        # set up the box for the flat hrs
        fbox, coef = flat_hrs.create_box(flat_hrs.flux, interp=interp)
        # for each line, correct for the shape of the flux
        y1 = 0
        y2 = len(fbox)
        for i in range(y1, y2):
            if np.any(fbox[i] > 0):
                fbox[i][fbox[i] == np.inf] = 0
                sf = nd.filters.median_filter(fbox[i], size=median_filter_size)
                sf[sf <= 0] = max(sf.mean(), 1)
                fbox[i] = fbox[i] / sf

        # return to the original shape
        flat_hrs.flux = flat_hrs.unravel_box(fbox)

    # step to prevent divide by zero
    flat_hrs.flux[flat_hrs.flux == 0] = 1

    # divide and normalize by the flux
    hrs.flux = hrs.flux * np.median(flat_hrs.flux) / flat_hrs.flux

    if hrs.error is not None:
        hrs.error = hrs.error * np.median(flat_hrs.flux) / flat_hrs.flux

    return hrs


def flatfield_science(ccd, flat_frame, order_frame, median_filter_size=None, interp=False):
    """Flatfield all of the orders in a science frame

    Parameters
    ----------
    ccd: ~ccdproc.CCDData
        Science frame to be flatfielded

    flar_frame: ~ccdproc.CCDData
        Frame containing the flat field for each of the orders

    order_frame: ~ccdproc.CCDData
        Frame containting the positions of each of the orders

     median_filter_size: None or int
        Size for median filter to be run on the data and remove the general
        flat field shape

    Returns
    -------
    ccd: ~ccdproc.CCDData
        Flatfielded science frame

    """
    # get a list of orders
    o1 = order_frame.data[order_frame.data > 0].min()
    o2 = order_frame.data.max()
    order_arr = np.arange(o1, o2, dtype=int)

    ndata = 0.0 * ccd.data
    if ccd.uncertainty is not None:
        edata = ccd.uncertainty.array
    else:
        edata = None

    for n_order in order_arr:
        hrs = HRSOrder(n_order)
        hrs.set_order_from_array(order_frame.data)
        hrs.set_flux_from_array(ccd.data, flux_unit=ccd.unit, error=edata)

        flat_hrs = HRSOrder(n_order)
        flat_hrs.set_order_from_array(order_frame.data)
        flat_hrs.set_flux_from_array(flat_frame.data, flux_unit=flat_frame.unit)

        hrs = flatfield_science_order(hrs, flat_hrs, median_filter_size=median_filter_size, interp=interp)

        ndata[hrs.region[0], hrs.region[1]] = hrs.flux
        if edata is not None: edata[hrs.region[0], hrs.region[1]] = hrs.error

    ccd.data = ndata
    return ccd
