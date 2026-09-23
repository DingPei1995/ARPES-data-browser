"""
loader/soleil.py
====================
The SOLEIL / ANTARES nano-ARPES loader: the beamline this program was
written for.

All of the actual parsing lives in :mod:`loader.nxs_file` -- the
``scan_data/data_NN`` numbering, the ``/ANTARES/...`` metadata tree, the
four entry layouts and their quirks, built against real files and covered by
its own tests. This module is the thin piece that presents it as one loader
among several, so that a second beamline can be added beside it without
anyone editing a 1300-line parser that already works.

**This is the file to copy when adding a beamline.** A new loader needs the
four things below -- a name, a structural ``can_open``, a listing, and a
read that returns an :class:`loader.nxs_file.NxsScan` -- and one
``register()`` call at the bottom. Nothing else in the program has to
change: the loader window builds its combo from the registry, and axis
order and the meaning of a scanned axis are handled once, for every loader,
in :mod:`loader.registry`.

What an NxsScan has to come back with
-------------------------------------
``kind``   one of ``cut`` (2-D), ``map`` (3-D), ``spem_1d`` (3-D),
           ``spem_4d`` (4-D).
axes       the vectors for that kind's dimensions, in the array's own
           order: ``x, y`` for a cut; ``x, k, z`` for a map;
           ``y, x, k, z`` for a 4-D spatial scan.
``value``  the array -- a numpy array, or a
           :class:`loader.nxs_file.LazyArray` wrapping the still-open
           dataset, which is preferable for anything large.
``labels`` per axis slot, *including units*; the viewers print these
           verbatim, so this is where units live.
``info``   free-form metadata, shown in the information panel.
"""
from __future__ import annotations


from loader.nxs_file import list_datasets, load_soleil_nxs, _classify_entry
from loader.nxs_file import HANDLES
from loader.registry import Loader, register


class SoleilAntaresLoader(Loader):
    name = "SOLEIL ANTARES"
    description = ("Nano-ARPES .nxs files from the ANTARES beamline at "
                   "SOLEIL: spatial scans, deflector maps and single cuts.")
    patterns = ("*.nxs",)
    priority = 10

    def can_open(self, path: str) -> bool:
        """Recognised by *structure*, not by filename.

        ``_classify_entry`` looks at what an entry actually contains -- which
        scan_data datasets exist and what shape the actuators are -- which is
        the same check the reader itself uses, so "can_open said yes" and
        "load worked" cannot disagree.
        """
        try:
            with HANDLES.borrow(path) as f:
                for name in f.keys():
                    try:
                        if _classify_entry(f, "/" + name) is not None:
                            return True
                    except Exception:
                        continue
        except Exception:
            return False
        return False

    def list_entries(self, path: str) -> list:
        return list_datasets(path)

    def load(self, path: str, entry: str = None):
        return load_soleil_nxs(path, entry=entry)


register(SoleilAntaresLoader())
