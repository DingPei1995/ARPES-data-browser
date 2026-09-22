"""
loader/native.py
====================
The loader for the program's own saved format -- the one every derived
dataset is written in (see :mod:`loader.session`) and the one "Save..."
produces.

It is a loader like any beamline's, with one difference that matters in
use: **it never needs to be chosen**. The format marks its entries with
``loader.nxs_file.NATIVE_ATTR``, so recognising one is reading an attribute
rather than guessing from a layout. That is why its priority sits above
every beamline loader: a file this program wrote should never be mistaken
for a beamline file whose structure it happens to resemble, and reopening
yesterday's session must not ask the user which synchrotron it came from.

The reading itself is :func:`loader.nxs_file.load_soleil_nxs`, which already
dispatches a native entry to its own parser; this class is the registry
entry that makes that reachable by name and by detection.
"""
from __future__ import annotations

import h5py

from loader.nxs_file import (NATIVE_ATTR, list_datasets, load_soleil_nxs)
from loader.registry import Loader, register


class NativeLoader(Loader):
    name = "This program's own format"
    description = ("Datasets saved by this program, and everything in a "
                   "session folder. Detected automatically -- no need to "
                   "pick a beamline.")
    #: Saving here writes .nxs, but a file that has been renamed, or moved
    #: through a system that prefers .h5, is the same file and opens the same
    #: way -- the format is recognised by a marker inside it, not by its name.
    patterns = ("*.nxs", "*.h5", "*.hdf5")
    #: Above every beamline loader: this one *knows*, the others guess.
    priority = 100

    def can_open(self, path: str) -> bool:
        """True if any top-level entry carries the native marker."""
        try:
            with h5py.File(path, "r") as f:
                return any(NATIVE_ATTR in f[name].attrs
                           for name in f.keys()
                           if isinstance(f[name], h5py.Group))
        except Exception:
            return False

    def list_entries(self, path: str) -> list:
        """The entries, each carrying the name the dataset was saved under.

        A beamline file's row is named after the file, because that is all
        there is to go on. A file written here knows better: the dataset had
        a name in the list when it was saved, and it is stored with it. That
        matters most for a session folder, whose filenames carry a counter
        and an extension that are this module's business and not the user's
        -- a recovered dataset should come back as "gold Fermi surface", not
        as "007-gold Fermi surface.nxs".
        """
        entries = list_datasets(path)
        try:
            with h5py.File(path, "r") as f:
                for info in entries:
                    group = f.get(info.get("entry") or "")
                    if group is None:
                        continue
                    saved = group.attrs.get("nxsloader_name")
                    if isinstance(saved, bytes):
                        saved = saved.decode("utf-8", "replace")
                    if saved:
                        info["name"] = str(saved)
        except Exception:
            pass        # the names are a nicety; the entries are the point
        return entries

    def load(self, path: str, entry: str = None):
        # Deliberately no ``progress``: declaring it is what marks a reader
        # as one that expects to take long enough to need a progress bar,
        # and routes it through the worker thread. An HDF5 open is one call
        # into the library with nothing to report.
        return load_soleil_nxs(path, entry=entry)


register(NativeLoader())
