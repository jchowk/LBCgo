"""Pixel masks and inverse-variance weight maps for LBC chip images.

Masks and weights are written as *sidecar* files next to each processed
image, with the same HDU layout as the image (an empty primary HDU followed
by one image extension per chip, in the same order):

    <image>.fits  ->  <image>.mask.fits    (uint8 bitmask, 0 = good)
                      <image>.weight.fits  (float32 inverse variance, 0 = bad)

The ``.weight.fits`` suffix matches the SExtractor/SWarp ``WEIGHT_SUFFIX``
default, so the files can be passed to them as ``MAP_WEIGHT`` images.

Sidecar primary headers deliberately omit IMAGETYP, FILTER and OBJECT so that
an ``ImageFileCollection`` over a working directory never mistakes them for
science frames.

Mask bits
---------
SATURATED  (1)   raw ADU >= saturation_fraction * SATURATE (grown by
                 ``saturation_grow`` pixels); flagged before overscan
                 subtraction.
BADPIX     (2)   normalized flat deviates from its local median by more
                 than ``badpix_threshold``, non-positive or non-finite flat,
                 or listed in a user bad-pixel region file.
VIGNETTED  (4)   normalized flat < ``vignette_threshold``.
NONFINITE  (8)   flat-fielded data value is NaN or inf.
COSMIC     (16)  pixel replaced by cosmic-ray cleaning.

Weights
-------
For flat-fielded data d = raw / f (f = flat normalized by its median, as in
``go_flatfield``) and a sky level S in flat-fielded ADU, the raw variance is
S*f/g + (RN/g)^2 [ADU^2], so the inverse variance of d is

    w = f^2 / (S*f/g + (RN/g)^2)

with g the gain (e-/ADU) and RN the read noise (e-). Source Poisson noise is
not included (the usual convention for coadd weights). w = 0 wherever the
mask is non-zero.
"""

import os
import shutil
import warnings

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy import ndimage

SATURATED = 1
BADPIX = 2
VIGNETTED = 4
NONFINITE = 8
COSMIC = 16

MASK_BITS = {'SATURATED': SATURATED, 'BADPIX': BADPIX,
             'VIGNETTED': VIGNETTED, 'NONFINITE': NONFINITE,
             'COSMIC': COSMIC}

# Fallbacks when headers lack the keywords (values from LBC-Blue headers).
DEFAULT_GAIN = 1.75       # e-/ADU
DEFAULT_RDNOISE = 12.0    # e-
DEFAULT_SATURATE = 65535.0


def sidecar_name(filename, kind):
    """Return the sidecar filename for ``kind`` ('mask' or 'weight')."""
    base, ext = os.path.splitext(filename)
    if ext != '.fits':
        raise ValueError("Expected a .fits filename, got {0}".format(filename))
    return '{0}.{1}.fits'.format(base, kind)


def header_value(key, headers, default, verbose=False):
    """Return the first value of ``key`` found in ``headers``, else default."""
    for hdr in headers:
        if hdr is not None and key in hdr:
            return float(hdr[key])
    if verbose:
        print('Warning: {0} not found in header; using {1}.'.format(key, default))
    return float(default)


def saturation_mask(raw_data, saturate, fraction=0.9, grow=1):
    """Boolean mask of saturated pixels in raw (pre-overscan) ADU.

    Parameters
    ----------
    raw_data : ndarray
        Raw chip data in ADU, before overscan subtraction.
    saturate : float
        Saturation level in ADU (``SATURATE`` header keyword).
    fraction : float, optional
        Flag pixels at or above ``fraction * saturate``. Default: 0.9
    grow : int, optional
        Dilate the mask by this many pixels. Default: 1

    Returns
    -------
    ndarray of bool
    """
    mask = raw_data >= fraction * saturate
    if grow > 0 and mask.any():
        mask = ndimage.binary_dilation(mask, iterations=grow)
    return mask


def flat_mask(flat_norm, vignette_threshold=0.5, badpix_threshold=0.2,
              badpix_box=5):
    """Mask bits derived from a median-normalized flat field.

    Parameters
    ----------
    flat_norm : ndarray
        Flat field divided by its median.
    vignette_threshold : float, optional
        Pixels with ``flat_norm`` below this are VIGNETTED. Default: 0.5
    badpix_threshold : float, optional
        Pixels whose ratio to the local ``badpix_box`` x ``badpix_box`` median
        differs from 1 by more than this are BADPIX. Catches bad columns and
        isolated dead/hot pixels, not features wider than ~badpix_box/2.
        Default: 0.2
    badpix_box : int, optional
        Size of the median filter box in pixels. Default: 5

    Returns
    -------
    ndarray of uint8
    """
    mask = np.zeros(flat_norm.shape, dtype=np.uint8)

    bad = ~np.isfinite(flat_norm) | (flat_norm <= 0)
    safe = np.where(bad, 1.0, flat_norm)

    local = ndimage.median_filter(safe, size=badpix_box, mode='nearest')
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = safe / local
    bad |= ~np.isfinite(ratio) | (np.abs(ratio - 1.0) > badpix_threshold)

    mask[bad] |= BADPIX
    mask[safe < vignette_threshold] |= VIGNETTED
    return mask


def sky_level(data, mask=None):
    """Sigma-clipped median of the unmasked pixels of ``data``."""
    good = np.isfinite(data)
    if mask is not None:
        good &= (mask == 0)
    if not good.any():
        return 0.0
    _, median, _ = sigma_clipped_stats(data[good], sigma=3.0, maxiters=5)
    return float(median)


def inverse_variance_weight(flat_norm, sky, gain, rdnoise, mask=None):
    """Background-limited inverse-variance weight for flat-fielded data.

    See the module docstring for the formula. ``sky`` is in flat-fielded ADU
    and is clipped at zero. Returns float32, 0 where ``mask`` is non-zero or
    the result is not finite.
    """
    sky = max(float(sky), 0.0)
    var_raw = sky * flat_norm / gain + (rdnoise / gain) ** 2
    with np.errstate(divide='ignore', invalid='ignore'):
        weight = flat_norm ** 2 / var_raw
    weight[~np.isfinite(weight) | (flat_norm <= 0)] = 0.0
    if mask is not None:
        weight[mask != 0] = 0.0
    return weight.astype(np.float32)


def read_badpix_regions(filename):
    """Read a bad-pixel region file.

    Each non-comment line is ``chip x1 x2 y1 y2``: 1-based, inclusive pixel
    ranges in trimmed (post-overscan) chip coordinates. Lines starting with
    ``#`` are ignored.

    Returns
    -------
    dict
        ``{chip: [(x1, x2, y1, y2), ...]}``
    """
    regions = {}
    with open(filename) as fh:
        for line in fh:
            line = line.split('#')[0].strip()
            if not line:
                continue
            chip, x1, x2, y1, y2 = (int(v) for v in line.split())
            regions.setdefault(chip, []).append((x1, x2, y1, y2))
    return regions


def apply_badpix_regions(mask, regions):
    """Set BADPIX in ``mask`` for a list of (x1, x2, y1, y2) regions."""
    for x1, x2, y1, y2 in regions:
        mask[y1 - 1:y2, x1 - 1:x2] |= BADPIX
    return mask


def _primary_header(source_file, kind):
    hdr = fits.Header()
    hdr['ORIGFILE'] = (os.path.basename(source_file), 'Image this file describes')
    hdr['SIDECAR'] = (kind, 'LBCgo sidecar type')
    if kind == 'mask':
        for name, bit in MASK_BITS.items():
            hdr['MASK_{0}'.format(name[:3])] = (bit, 'Mask bit: {0}'.format(name))
    return hdr


def write_sidecar(source_file, kind, arrays, extnames=None, headers=None,
                  output_file=None):
    """Write a mask or weight sidecar mirroring ``source_file``'s HDU layout.

    Parameters
    ----------
    source_file : str
        The image file the sidecar describes; the sidecar name is derived
        from it unless ``output_file`` is given.
    kind : {'mask', 'weight'}
    arrays : list of ndarray
        One array per image extension, in HDU order.
    extnames : list of str, optional
        EXTNAME for each extension.
    headers : list of fits.Header, optional
        Extra cards for each extension (e.g. SKYLEVEL for weights).
    output_file : str, optional
        Explicit output filename.

    Returns
    -------
    str
        The filename written.
    """
    if output_file is None:
        output_file = sidecar_name(source_file, kind)
    dtype = np.uint8 if kind == 'mask' else np.float32

    hdul = fits.HDUList([fits.PrimaryHDU(header=_primary_header(source_file, kind))])
    for i, arr in enumerate(arrays):
        hdr = fits.Header()
        if headers is not None and headers[i] is not None:
            hdr.extend(headers[i])
        hdu = fits.ImageHDU(data=np.asarray(arr, dtype=dtype), header=hdr)
        if extnames is not None and extnames[i]:
            hdu.header['EXTNAME'] = extnames[i]
        hdul.append(hdu)
    hdul.writeto(output_file, overwrite=True)
    return output_file


def read_sidecar(source_file, kind, index):
    """Return extension ``index`` of the sidecar for ``source_file``, or None."""
    name = sidecar_name(source_file, kind)
    if not os.path.exists(name):
        return None
    with fits.open(name) as hdul:
        return np.array(hdul[index].data)


def extract_sidecar_chips(source_file, target_file, index):
    """Write extension ``index`` of each sidecar of ``source_file`` as the
    single-extension sidecar of ``target_file`` (used when splitting chips).
    """
    for kind in ('mask', 'weight'):
        name = sidecar_name(source_file, kind)
        if not os.path.exists(name):
            continue
        with fits.open(name) as hdul:
            ext = hdul[index]
            hdr = ext.header.copy()
            for key in ('XTENSION', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2',
                        'PCOUNT', 'GCOUNT', 'EXTNAME', 'BSCALE', 'BZERO'):
                hdr.remove(key, ignore_missing=True)
            write_sidecar(target_file, kind, [ext.data],
                          extnames=[ext.header.get('EXTNAME')],
                          headers=[hdr])


def move_sidecars(source_file, destination):
    """Move any mask/weight sidecars of ``source_file`` to ``destination``.

    ``destination`` is a directory; existing files there are overwritten.
    """
    for kind in ('mask', 'weight'):
        name = sidecar_name(source_file, kind)
        if os.path.exists(name):
            dest = os.path.join(destination, os.path.basename(name))
            if os.path.exists(dest):
                os.remove(dest)
            shutil.move(name, destination)


def copy_sidecars(source_file, target_file):
    """Copy any mask/weight sidecars of ``source_file`` to those of ``target_file``."""
    for kind in ('mask', 'weight'):
        name = sidecar_name(source_file, kind)
        if os.path.exists(name):
            shutil.copyfile(name, sidecar_name(target_file, kind))
