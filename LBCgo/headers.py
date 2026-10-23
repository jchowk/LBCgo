"""Header scanning that tolerates invalid FITS header cards.

Raw LBC frames can carry cards astropy cannot parse (e.g. an unquoted
``PA_PNT = nan``). ``ccdproc.ImageFileCollection``
parses every card of every header, so one such card stops the whole scan,
and ``HDUList.writeto`` refuses to write a header that still holds it.

- `HeaderCollection` is a drop-in for the part of the
  ``ImageFileCollection`` API that LBCgo uses (``location``, ``files``,
  ``summary``, ``files_filtered``, ``values``). Only the requested keywords
  are parsed; a keyword that is absent or unparsable is masked in
  ``summary``, as ``ImageFileCollection`` does for absent keywords.
- `clean_header` drops unparsable cards from a header before it is copied
  into an output file, and records what it dropped in HISTORY.
- `read_primary_header` is ``fits.getheader`` followed by `clean_header`.
"""

import os
import warnings
from glob import glob

import numpy as np
from astropy.io import fits
from astropy.io.fits.verify import VerifyError, VerifyWarning
from astropy.table import MaskedColumn, Table

__all__ = ['HeaderCollection', 'clean_header', 'read_primary_header',
           'header_values']


def _usable_card(card):
    """Copy of `card`, fixed where astropy can fix it; None if it cannot.

    A card whose value does not parse is not usable, even though
    ``verify('fix')`` would "fix" it: it turns ``PA_PNT = nan`` into the
    string ``'nan'``, which later code would read as a value.
    """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', VerifyWarning)
            # Value first: reading card.image already applies the fix.
            card.value
            new = fits.Card.fromstring(card.image)
            new.verify('fix+exception')
        return new
    except (VerifyError, ValueError):
        return None


def clean_header(header, verbose=False):
    """Copy of `header` without the cards astropy cannot parse.

    Cards astropy can fix (as ``verify('fix')`` would) are kept, fixed.
    Each dropped card is recorded as a HISTORY entry naming its keyword, so
    an output header shows that something was removed. With `verbose`, the
    dropped keywords are printed.
    """
    out = fits.Header()
    dropped = []
    for card in header.cards:
        new = _usable_card(card)
        if new is not None:
            out.append(new, end=True)
        else:
            dropped.append(card.keyword)
    for key in dropped:
        out.add_history('LBCgo: dropped unparsable card {0}'.format(key))
    if dropped and verbose:
        print('Dropped unparsable header card(s): {0}'.format(
            ', '.join(dropped)))
    return out


def read_primary_header(path, verbose=False):
    """Primary header of `path` with unparsable cards removed."""
    return clean_header(fits.getheader(path), verbose=verbose)


def header_values(header, keywords):
    """Dict of lower-case keyword -> value for the keywords in `header`.

    Keywords that are absent, unparsable or not scalar are left out.
    """
    row = {}
    for key in keywords:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', VerifyWarning)
                value = header[key]
        except (KeyError, VerifyError, ValueError):
            continue
        if isinstance(value, (str, bool, int, float, np.number)):
            row[key.lower()] = value
    return row


def _column(name, values):
    """Column from a list holding None for missing entries (masked)."""
    present = [v for v in values if v is not None]
    if any(isinstance(v, str) for v in present):
        present = [str(v) for v in present]
        filled = [str(v) if v is not None else '' for v in values]
    else:
        fill = present[0] if present else ''
        filled = [v if v is not None else fill for v in values]
    mask = [v is None for v in values]
    return MaskedColumn(filled, name=name, mask=mask)


class HeaderCollection:
    """Summary of the primary headers of a set of FITS files.

    Parameters
    ----------
    location : str
        Directory holding the files.
    keywords : list of str
        Header keywords to read; column names in `summary` are lower case.
    filenames : list of str, optional
        File names (relative to `location`). If not given, the files
        matching `glob_include` are used.
    glob_include : str, optional
        Glob pattern for the files in `location`. Default: ``'*.fits*'``.
    """

    def __init__(self, location='./', keywords=None, filenames=None,
                 glob_include=None):
        self.location = location
        self.keywords = [k.lower() for k in (keywords or [])]
        if filenames is None:
            pattern = glob_include or '*.fits*'
            filenames = sorted(os.path.basename(f) for f in
                               glob(os.path.join(location, pattern)))
        self._files = [str(f) for f in filenames]

        rows = [header_values(fits.getheader(os.path.join(location, f)),
                              self.keywords) for f in self._files]
        columns = [MaskedColumn(self._files, name='file', dtype=str)]
        for key in self.keywords:
            columns.append(_column(key, [r.get(key) for r in rows]))
        self.summary = Table(columns, masked=True)

    def __len__(self):
        return len(self._files)

    @property
    def files(self):
        """List of file names in the collection."""
        return list(self._files)

    def _match(self, key, value):
        """Boolean array: rows whose `key` matches `value`.

        Strings match case-insensitively, ignoring surrounding whitespace;
        ``'*'`` matches any row where the keyword is present. Missing
        keywords never match.
        """
        key = key.lower()
        n = len(self._files)
        if key not in self.summary.colnames:
            return np.zeros(n, dtype=bool)
        col = self.summary[key]
        present = ~np.ma.getmaskarray(col)
        if isinstance(value, str) and value == '*':
            return present
        data = col.data.data if hasattr(col.data, 'data') else col.data
        if isinstance(value, str):
            wanted = value.strip().lower()
            hits = np.array([str(v).strip().lower() == wanted for v in data],
                            dtype=bool)
        else:
            hits = np.array([v == value for v in data], dtype=bool)
        return hits & present

    def files_filtered(self, **kwd):
        """Array of file names whose headers match all keyword=value pairs."""
        keep = np.ones(len(self._files), dtype=bool)
        for key, value in kwd.items():
            keep &= self._match(key, value)
        return np.array(self._files, dtype=str)[keep]

    def values(self, keyword, unique=False):
        """Values of `keyword`, one per file (missing -> None).

        With `unique`, the sorted distinct values of the files that have the
        keyword.
        """
        col = self.summary[keyword.lower()]
        mask = np.ma.getmaskarray(col)
        data = col.data.data if hasattr(col.data, 'data') else col.data
        vals = [None if m else v.item() if hasattr(v, 'item') else v
                for v, m in zip(data, mask)]
        if unique:
            return sorted({v for v in vals if v is not None})
        return vals
