"""
loader/session.py
==============
Where a dataset computed in this session actually lives.

The problem this solves
-----------------------
Until now a k-map, an arbitrary cut, a processed map -- anything the program
worked out rather than read -- existed *only* as a numpy array held by the
browser's list. Two consequences followed, and both were routinely painful:

* **Killing the program lost all of it.** There was no copy anywhere, and
  the only way to make one was to remember to use "Save...". A run that
  froze (which, with every computation on the GUI thread, looked identical
  to a crash) took a morning's processing with it.
* **Nothing could be unloaded.** Each derived cube stayed resident for the
  whole session; twenty of them at ~75 MB each is 1.5 GB, and the program
  had no notion of a budget, an eviction policy or a ceiling.

So every derived dataset is written to a **session folder** the moment it
exists, in the program's own format (:func:`loader.nxs_file.save_dataset`,
which is compressed, keeps the array's dtype, and is recognised on sight
when read back). The list then holds a *path*, exactly as it already did
for datasets read from a beamline file, and the array itself is a cache
that can be dropped and re-read. A dataset read from a file and a dataset
computed here stop being two different kinds of thing.

What this module is not
-----------------------
It is not an undo history and not a project file. The session folder is a
working copy: it is pruned after :data:`KEEP_DAYS`, and "Save..." remains
how the user puts results somewhere they choose. What it guarantees is only
that nothing is ever lost *because the program stopped*.

No Qt in here, so it can be (and is) tested directly.
"""
from __future__ import annotations

import os
import re
import shutil
import time
from collections import OrderedDict

import numpy as np

from loader.nxs_file import save_dataset, list_datasets, load_soleil_nxs

__all__ = ["SessionStore", "MemoryBudget", "DEFAULT_ROOT", "KEEP_DAYS",
           "LOG_NAME", "dataset_dict", "scan_to_dict", "read_log",
           "leftover_sessions"]

#: Where session folders go. Beside the config file the program already
#: keeps, so there is one place to look (and to clear out).
DEFAULT_ROOT = os.path.join(os.path.expanduser("~"), ".arpes_viewer", "sessions")

#: Where sessions were kept when the program was called NxsLoader. Still
#: looked in when offering to recover unsaved work, because an upgrade must
#: not be the reason somebody's crashed session becomes unfindable -- the
#: whole point of the session folder is that stopping the program does not
#: lose anything, and renaming the program is a way of stopping it.
LEGACY_ROOTS = (os.path.join(os.path.expanduser("~"), ".nxsloader", "sessions"),)

#: Session folders older than this are removed at startup.
#:
#: Three days rather than a week: a working copy is only worth keeping until
#: the result has been saved somewhere of the user's own, and a copy that is
#: still around after three days is one nobody came back for. A dataset that
#: *is* saved has its copy removed at that moment (see
#: :meth:`SessionStore.released`), so this only governs the ones a crash
#: left behind.
KEEP_DAYS = 3

#: The operations log, one file for every session, in the sessions root so
#: that finding it after a crash is "open this one file" rather than
#: "work out which folder was mine".
LOG_NAME = "operations.log"

#: How much derived data may stay resident before the least recently used is
#: dropped (it can be re-read from its file). Not a hard cap on the process:
#: a dataset open in a viewer is held by that window too, and stays until
#: the window closes.
DEFAULT_BUDGET_BYTES = 1_500_000_000


def _sanitise(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._ -]+", "_", str(name)).strip()
    return cleaned[:80] or "dataset"


def scan_to_dict(kind: str, scan, name: str) -> dict:
    """One dataset in the plain ``{name, kind, axes, labels, value, info,
    motors}`` form :func:`loader.nxs_file.save_dataset` takes.

    Which axes a kind has comes from ``loader.nxs_file.AXIS_SLOTS``, so this
    works for any kind listed there and needs no edit when one is added --
    it used to be a fourth hand-written copy of that table.

    The **constructor** order is used, not the array order: it is what the
    saved file's ``axis_x``/``axis_y``/... are read back as, and for a 4-D
    spatial scan the two differ.
    """
    from loader.nxs_file import axis_slots
    slots = axis_slots(kind, "constructor")
    if not slots:
        raise ValueError(f"{kind!r} has no plain axes-and-array form")
    axes = tuple(getattr(scan, slot, None) for slot in slots)
    value = scan.array() if hasattr(scan, "array") else scan.value
    if any(axis is None for axis in axes) or value is None:
        raise ValueError(f"{kind!r} dataset is missing an axis or its values")
    return {"name": name, "kind": kind, "axes": axes, "labels": dict(scan.labels),
            "value": value, "info": dict(scan.info),
            "motors": dict(getattr(scan, "fourd_info", {}) or {})}


def dataset_dict(data, name: str = None) -> dict:
    """:func:`scan_to_dict` for anything with ``.kind`` and ``.scan`` --
    an :class:`NxsData` or a :class:`MemoryData`."""
    return scan_to_dict(data.kind, data.scan,
                        name or getattr(data, "source_label", "dataset"))


class MemoryBudget:
    """A least-recently-used cache of loaded arrays, with a byte ceiling.

    Eviction only drops *this* cache's reference. A dataset being looked at
    is also referenced by its viewer window, so evicting it frees nothing
    and breaks nothing -- it simply stops this cache from being the reason
    it is resident. That is the intended behaviour: the budget governs what
    is kept *on spec*, not what is in use.
    """

    def __init__(self, limit_bytes: int = DEFAULT_BUDGET_BYTES):
        self.limit_bytes = int(limit_bytes)
        self._entries: "OrderedDict[str, tuple]" = OrderedDict()

    @staticmethod
    def _size_of(data) -> int:
        value = getattr(getattr(data, "scan", None), "value", None)
        if value is None:
            return 0
        return int(getattr(value, "nbytes", 0) or 0)

    # The budget *owns a reference* to what it caches, like any other
    # holder: it takes one in put() and gives it back when the entry goes.
    # It used to hand its one object to every caller without counting them,
    # so the first viewer to close released the only reference and closed
    # the file -- and the budget went on handing out the dead object, which
    # failed on the next read with "identifier is not of specified type".
    @staticmethod
    def _retain(data):
        retain = getattr(data, "retain", None)
        return retain() if callable(retain) else data

    @staticmethod
    def _release(data):
        close = getattr(data, "close", None)
        if callable(close):
            try:
                close()
            except Exception:                               # noqa: BLE001
                pass

    def get(self, key: str):
        """The cached dataset, or None. The budget keeps its own reference;
        a caller that keeps the object should :meth:`retain` it."""
        entry = self._entries.get(key)
        if entry is None:
            return None
        alive = getattr(entry[0], "alive", None)
        if callable(alive) and not alive():
            self.discard(key)
            return None
        self._entries.move_to_end(key)
        return entry[0]

    def put(self, key: str, data):
        self.discard(key)
        self._entries[key] = (self._retain(data), self._size_of(data))
        self._evict()
        return data

    def discard(self, key: str):
        entry = self._entries.pop(key, None)
        if entry is not None:
            self._release(entry[0])

    def clear(self):
        for key in list(self._entries):
            self.discard(key)

    def total_bytes(self) -> int:
        return sum(size for _data, size in self._entries.values())

    def _evict(self):
        while len(self._entries) > 1 and self.total_bytes() > self.limit_bytes:
            _key, (data, _size) = self._entries.popitem(last=False)
            self._release(data)


class SessionStore:
    """The folder this run's derived datasets are written to.

    One file per dataset rather than one file for the session: a dataset can
    then be re-read, replaced or discarded on its own, and a half-written
    file costs one result instead of all of them.
    """

    #: Folder names handed out in this process, so that two stores made
    #: before either has written anything still differ.
    _taken = set()

    def __init__(self, root: str = None, budget_bytes: int = DEFAULT_BUDGET_BYTES):
        self.root = root or DEFAULT_ROOT
        self.name = self._unique_name()
        self.folder = os.path.join(self.root, self.name)
        self.budget = MemoryBudget(budget_bytes)
        self._counter = 0
        self._ready = False

    def _unique_name(self) -> str:
        """``<date>-<time>-<pid>``, with a suffix if that is taken.

        The timestamp has one-second resolution and the pid is reused by the
        operating system, so "restarted within the same second" and "two
        stores in one process" both collide -- and a collision means two
        sessions sharing a folder, each treating the other's unsaved work as
        its own to prune. Rare, and silently destructive, so it is checked
        rather than assumed away.
        """
        base = f"{time.strftime('%Y%m%d-%H%M%S')}-{os.getpid()}"
        name, suffix = base, 1
        while name in SessionStore._taken or os.path.exists(
                os.path.join(self.root, name)):
            suffix += 1
            name = f"{base}-{suffix}"
        SessionStore._taken.add(name)
        return name

    # -- the folder ---------------------------------------------------------
    def _ensure(self):
        if not self._ready:
            os.makedirs(self.folder, exist_ok=True)
            self._ready = True
        return self.folder

    # -- the log --------------------------------------------------------------
    @property
    def log_path(self) -> str:
        return os.path.join(self.root, LOG_NAME)

    def log(self, action: str, name: str = "", path: str = "", note: str = ""):
        """Append one line to the operations log.

        The log exists for one situation: the program stopped, and the work
        is somewhere under ``~/.arpes_viewer/sessions`` under a filename nobody
        chose. Rather than making the user go looking, every write, every
        save and every deletion is recorded with its path, so the answer is
        in one text file in the order things happened.

        Deliberately plain text with fixed ``|``-separated fields -- readable
        in any editor on a beamline machine with nothing installed, and still
        parseable (:func:`read_log`). Never raises: a log that cannot be
        written must not stop data being saved, which is the thing it is
        there to protect.
        """
        line = " | ".join([
            time.strftime("%Y-%m-%dT%H:%M:%S"),
            self.name,
            f"{action:<8}",
            str(name),
            str(path),
            str(note),
        ])
        try:
            os.makedirs(self.root, exist_ok=True)
            with open(self.log_path, "a", encoding="utf-8") as handle:
                handle.write(line + "\n")
        except OSError:
            pass
        return line

    def prune_old(self, keep_days: int = KEEP_DAYS) -> int:
        """Delete session folders older than ``keep_days``. Returns how many
        went. Never touches this session's own folder, and records what it
        removed in the log -- a sweep that silently deleted a folder
        somebody was about to go looking for would defeat the point of
        keeping a log at all."""
        if not os.path.isdir(self.root):
            return 0
        cutoff = time.time() - keep_days * 86400
        removed = 0
        for name in sorted(os.listdir(self.root)):
            path = os.path.join(self.root, name)
            if path == self.folder or not os.path.isdir(path):
                continue
            try:
                if os.path.getmtime(path) < cutoff:
                    files = len([f for f in os.listdir(path) if f.endswith(".nxs")])
                    shutil.rmtree(path, ignore_errors=True)
                    self.log("PRUNE", name, path,
                             f"older than {keep_days} days, {files} dataset(s)")
                    removed += 1
            except OSError:
                continue
        return removed

    def leftovers(self, include_legacy: bool = True):
        """What earlier sessions left behind: see :func:`leftover_sessions`.

        Folders from before the program was renamed are included, so that
        upgrading is not a way to lose track of a crashed session's work.
        """
        found = leftover_sessions(self.root, exclude=self.folder)
        if include_legacy and self.root == DEFAULT_ROOT:
            for legacy in LEGACY_ROOTS:
                found.extend(leftover_sessions(legacy))
        return sorted(found, key=lambda entry: -entry["modified"])

    def bytes_on_disk(self) -> int:
        if not os.path.isdir(self.folder):
            return 0
        total = 0
        for name in os.listdir(self.folder):
            try:
                total += os.path.getsize(os.path.join(self.folder, name))
            except OSError:
                pass
        return total

    # -- storing ------------------------------------------------------------
    def store(self, item: dict, progress=None) -> str:
        """Write one dataset (in :func:`scan_to_dict` form) and return its
        path. The name only shapes the filename; uniqueness comes from a
        counter, so two datasets called the same thing do not collide."""
        self._ensure()
        self._counter += 1
        path = os.path.join(self.folder,
                            f"{self._counter:03d}-{_sanitise(item.get('name'))}.nxs")
        save_dataset(path, [item], progress=progress)
        try:
            size = os.path.getsize(path)
        except OSError:
            size = 0
        self.log("STORE", item.get("name", ""), path, f"{size / 1e6:.1f} MB")
        return path

    def store_data(self, data, name: str = None, progress=None) -> str:
        return self.store(dataset_dict(data, name), progress=progress)

    def released(self, path: str, saved_to: str):
        """The dataset at ``path`` now lives in a file the user chose, so
        the working copy is not needed any more.

        Separate from :meth:`discard` so the log says *why* the copy went.
        Deleting it at the moment it becomes redundant -- rather than
        waiting for the retention sweep -- is what keeps the session folder
        from being a pile of duplicates of things already saved.
        """
        self.log("SAVED", os.path.basename(path), saved_to,
                 "working copy removed")
        self._remove(path)

    def discard(self, path: str, reason: str = "removed from the list"):
        self.log("DISCARD", os.path.basename(path), path, reason)
        self._remove(path)

    def _remove(self, path: str):
        self.budget.discard(path)
        try:
            os.remove(path)
        except OSError:
            pass

    # -- reading back --------------------------------------------------------
    def load_scan(self, path: str):
        """The :class:`NxsScan` stored at ``path``. Its array is a
        :class:`loader.nxs_file.LazyArray`, so this is cheap and stays cheap
        until something asks for the numbers."""
        return load_soleil_nxs(path)

    def entries(self, path: str):
        return list_datasets(path)


def read_log(root: str = None, limit: int = 500):
    """The operations log, newest last, as ``(when, session, action, name,
    path, note)`` tuples.

    Parsing back what :meth:`SessionStore.log` wrote. Lines that do not
    have the expected shape are skipped rather than raising: a log is a
    record, and a corrupted tail (a crash mid-write) must not stop the rest
    of it being read -- which is exactly the situation it is consulted in.
    """
    path = os.path.join(root or DEFAULT_ROOT, LOG_NAME)
    if not os.path.isfile(path):
        return []
    rows = []
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parts = [field.strip() for field in line.rstrip("\n").split(" | ")]
                if len(parts) >= 5:
                    rows.append(tuple(parts[:6] + [""] * (6 - len(parts))))
    except OSError:
        return []
    return rows[-limit:] if limit else rows


def leftover_sessions(root: str = None, exclude: str = None):
    """Earlier sessions whose folders still hold datasets.

    After a crash this is the list of "work that was never saved anywhere
    you chose" -- each entry is ``{folder, name, files, bytes, modified}``,
    newest first, so the program can offer to put them back in the list
    rather than leaving the user to find a folder named after a timestamp
    and a process id.

    A session that shut down cleanly leaves nothing here: saving a dataset
    removes its working copy, and removing its row does too.
    """
    root = root or DEFAULT_ROOT
    if not os.path.isdir(root):
        return []
    found = []
    for name in os.listdir(root):
        folder = os.path.join(root, name)
        if not os.path.isdir(folder) or folder == exclude:
            continue
        try:
            files = sorted(f for f in os.listdir(folder) if f.endswith(".nxs"))
        except OSError:
            continue
        if not files:
            continue
        total = 0
        for filename in files:
            try:
                total += os.path.getsize(os.path.join(folder, filename))
            except OSError:
                pass
        found.append({"folder": folder, "name": name, "files": files,
                      "bytes": total,
                      "modified": os.path.getmtime(folder)})
    return sorted(found, key=lambda entry: -entry["modified"])


def estimate_bytes(*arrays) -> int:
    """Total footprint of some arrays, counting a lazy one as zero (it is
    not resident). For reporting, not for allocation."""
    total = 0
    for array in arrays:
        if array is None:
            continue
        if isinstance(array, np.ndarray):
            total += array.nbytes
    return total
