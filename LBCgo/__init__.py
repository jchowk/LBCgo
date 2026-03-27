"""LBCgo: Data reduction package for the Large Binocular Telescope's Large Binocular Camera (LBC).

Requires astromatic.net software sextractor, scamp, swarp.
"""

from importlib.metadata import version, PackageNotFoundError
try:
    __version__ = version("LBCgo")
except PackageNotFoundError:
    # Package not installed (e.g., running from source without pip install)
    __version__ = "unknown"

# Re-export primary functions to make them directly available
from .lbcproc import lbcgo
from .lbcregister import go_sextractor, go_register

# Import everything for backward compatibility
from .lbcproc import *
from .lbcregister import *

# Define public API
__all__ = ['lbcgo', 'go_sextractor', 'go_register']
