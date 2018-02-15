#!/usr/bin/env python
from astropy.io import fits
import sys, os
from subprocess import call
import textwrap
from glob import glob

from barak.io import loadobj, saveobj, writetxt
from astropy.io import ascii


def makedir(dirname, clean=False):
    """ Make a directory only if it doesn't exist.

    If clean is True, remove any files in the directory if it already
    exists."""
    if not os.path.lexists(dirname):
        print(('Creating %s/' % dirname))
        os.makedirs(dirname)
    else:
        if clean:
            call('rm %s/*' % dirname, shell=1)



