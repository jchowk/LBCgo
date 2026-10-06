from ccdproc import ImageFileCollection
ic0b = ImageFileCollection('./', keywords=keywds,
                                 glob_include='lbcb*fits.gz')
ic0r = ImageFileCollection('./', keywords=keywds,
                                 glob_include='lbcr*fits.gz')