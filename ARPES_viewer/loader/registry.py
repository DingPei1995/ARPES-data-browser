"""
loader/registry.py
==============
Which beamline a file came from, and who knows how to read it.

Until now "reading a file" meant one function, ``load_soleil_nxs``, whose
name was honest: it knows SOLEIL/ANTARES's ``scan_data/data_09`` numbering
and its ``/ANTARES/...`` metadata paths, and nothing else. Every other
beamline writes a different layout, so adding one meant editing that
function and hoping not to break the first.

This module makes the reader a choice instead. A loader is a small object
that answers three questions -- *can you open this?*, *what is in it?*, *read
this entry* -- and each one lives in its own file (``loader/soleil.py``,
``loader/native.py``), so a new beamline is a new file plus one
``register()`` call, touching nothing that already worked.

Detection versus choosing
-------------------------
:func:`detect` asks each loader whether it recognises a file, in
:attr:`Loader.priority` order. The program's own saved format is recognised
outright -- it is self-describing, which is the whole point of saving in it
-- so a dataset that came back from a session folder never needs the user to
say where it came from. A beamline file is a guess, and the loader dialog
shows which loader was guessed and lets it be overridden, because a wrong
guess on an unfamiliar layout should be a correction rather than a dead end.

Axis options
------------
Two things a reader cannot know but the user does, both carried in
:class:`LoadOptions`:

* **The axis order.** Two beamlines can record the same measurement with
  the array's dimensions in a different order. Rather than a per-loader
  transpose nobody can audit, the loader reports what it found and the user
  can permute it, once, at the door.
* **What the scanned axis actually is.** A "map" here means angle vs angle
  vs energy, but the same acquisition is used to scan photon energy,
  temperature or a gate voltage against the analyser's two axes. The data is
  identical in shape; only the label and the physics differ. Saying so at
  load time keeps the wrong unit from propagating into every later step --
  and stops a k-conversion being offered for an axis that is not an angle.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np

__all__ = ["Loader", "LoadOptions", "AXIS_ROLES", "register", "loaders",
           "get_loader", "detect", "load", "list_entries", "apply_options"]


#: What the first axis of a 3-D "map" can actually be. The key is what gets
#: stored in ``scan.info["axis0.role"]``; the value is the label the viewers
#: then show, with its unit.
#:
#: "deflector angle" is the default and the only one a k-conversion makes
#: sense for -- converting a temperature axis to momentum is meaningless,
#: and the viewer uses this to say so rather than producing a plausible
#: wrong answer.
AXIS_ROLES = {
    "angle": ("Angle (deflector)", "deg", True),
    "photon_energy": ("Photon energy", "eV", False),
    "temperature": ("Temperature", "K", False),
    "gate_voltage": ("Gate voltage", "V", False),
    "delay": ("Delay", "ps", False),
    "position": ("Position", "mm", False),
    "other": ("Scanned variable", "", False),
}


@dataclass
class LoadOptions:
    """What the user chose in the loader dialog, alongside the file itself.

    ``permutation`` is applied to the array's dimensions *after* the loader
    has read it, together with the matching axes, so that "swap the first
    two axes" means exactly that and needs no cooperation from the loader.
    ``None`` leaves the order alone.

    ``axis0_role`` of ``None`` means *whatever the loader worked out*. Some
    loaders can tell: a CASSIOPEE folder is a photon-energy scan exactly when
    the monochromator readback varies across its members, which the reader
    knows and the user should not have to restate. ``None`` keeps that; any
    other value is the user overriding it, which stays possible because a
    reader's guess is still a guess.
    """
    loader: str = None                 # None -> detect
    permutation: tuple = None          # e.g. (1, 0, 2)
    axis0_role: str = None             # a key of AXIS_ROLES; None -> as read
    axis0_label: str = ""              # overrides the role's own label
    extra: dict = field(default_factory=dict)


class Loader:
    """One beamline's reader.

    Subclasses (or any object with these four attributes) provide:

    ``name``         what the dialog's combo shows, e.g. "SOLEIL ANTARES".
    ``patterns``     filename globs, for the file dialog's filter.
    ``priority``     detection order; higher goes first. The self-describing
                     native format sits above every beamline guess.
    ``can_open(p)``  cheap structural check -- open the file, look, close.
    ``list_entries(p)``  one dict per openable dataset, as
                     ``loader.nxs_file.list_datasets`` returns them.
    ``load(p, entry, progress)``   an ``NxsScan``.

    ``progress`` is ``callable(done, total, label)`` or ``None``, and a
    loader may ignore it: an HDF5 read is one call into the library and has
    nothing to report. It exists for the readers where a load is genuinely
    long -- assembling a CASSIOPEE folder means parsing a hundred text files,
    fifteen seconds of it -- so the window can say what it is doing rather
    than appearing to hang. Accepting it is optional, and :func:`load` checks
    before passing it, so a loader written against the older three-argument
    signature keeps working.
    """

    name = "unnamed"
    description = ""
    patterns = ("*.nxs",)
    priority = 0

    def can_open(self, path: str) -> bool:
        raise NotImplementedError

    def list_entries(self, path: str) -> list:
        raise NotImplementedError

    def load(self, path: str, entry: str = None, progress=None):
        raise NotImplementedError


_REGISTRY: "dict[str, Loader]" = {}


def register(loader: Loader) -> Loader:
    """Add a loader. Importing a ``nxs_loader_*`` module is what calls this,
    so adding a beamline is a new file and an import, not an edit here."""
    _REGISTRY[loader.name] = loader
    return loader


def loaders() -> list:
    """Every registered loader, best-detection-chance first."""
    return sorted(_REGISTRY.values(), key=lambda l: -l.priority)


def get_loader(name: str) -> Loader:
    try:
        return _REGISTRY[name]
    except KeyError:
        raise ValueError(
            f"no loader called {name!r}; have "
            f"{sorted(_REGISTRY)}") from None


def detect(path: str) -> "Loader | None":
    """The first loader that recognises ``path``, or None.

    Never raises: a loader that throws while sniffing an unfamiliar file is
    treated as "not mine", because one broken loader must not stop the
    others from being asked.
    """
    for loader in loaders():
        try:
            if loader.can_open(path):
                return loader
        except Exception:
            continue
    return None


def list_entries(path: str, options: LoadOptions = None) -> list:
    """The datasets in ``path``, with the chosen (or detected) loader."""
    loader = _resolve(path, options)
    if loader is None:
        return []
    try:
        return loader.list_entries(path)
    except Exception:
        return []


def _resolve(path: str, options: LoadOptions = None) -> "Loader | None":
    if options is not None and options.loader:
        return get_loader(options.loader)
    return detect(path)


def load(path: str, entry: str = None, options: LoadOptions = None,
         progress=None):
    """Read one dataset, apply the load-time axis options, and return the
    :class:`loader.nxs_file.NxsScan`.

    ``progress`` is handed to the loader only if it accepts one, so that
    adding the argument here did not have to break every loader at once.
    """
    options = options or LoadOptions()
    loader = _resolve(path, options)
    if loader is None:
        raise ValueError(
            f"no loader recognises {os.path.basename(path)}. Choose one "
            f"explicitly in the Loader window if you know which it is.")
    if progress is not None and _accepts_progress(loader.load):
        scan = loader.load(path, entry, progress=progress)
    else:
        scan = loader.load(path, entry)
    scan.info.setdefault("loader.name", loader.name)
    return apply_options(scan, options)


def _accepts_progress(method) -> bool:
    import inspect
    try:
        parameters = inspect.signature(method).parameters
    except (TypeError, ValueError):                         # pragma: no cover
        return False
    return "progress" in parameters or any(
        p.kind is inspect.Parameter.VAR_KEYWORD for p in parameters.values())


# -- load-time axis handling ---------------------------------------------------
def _axis_slots(kind: str) -> tuple:
    """Which of NxsScan's axis attributes hold this kind's dimensions, in
    the array's own order -- from the one table that says so."""
    from loader.nxs_file import axis_slots
    return axis_slots(kind, "array")


def apply_options(scan, options: LoadOptions):
    """Permute the axes and label the scanned axis, as the dialog asked.

    Done here rather than inside each loader so that every beamline gets the
    same behaviour from the same code, and so a loader author has one less
    thing to get right.

    The user's choice wins where they made one. Where they left it at "as the
    file says", a role the loader established is kept -- and its label with
    it, because a loader that knows the axis is a photon energy has already
    written a better label than a generic one. Anything that said nothing is
    an angle, which is what every file read before roles existed was.
    """
    from loader.nxs_file import CUBE_KINDS

    if options.permutation:
        scan = permute_axes(scan, options.permutation)
    if scan.kind not in CUBE_KINDS:
        return scan

    chosen = options.axis0_role
    if not chosen:
        # No choice made: keep whatever the reader worked out, and leave its
        # label alone. Relabelling here would overwrite "Photon energy (eV)"
        # with the generic wording for the same role, which is worse.
        scan.info["axis0.role"] = str(scan.info.get("axis0.role", "angle"))
        return scan

    scan.info["axis0.role"] = chosen
    if chosen != "angle" or options.axis0_label:
        label, unit, _convertible = AXIS_ROLES.get(chosen, AXIS_ROLES["other"])
        text = options.axis0_label or (f"{label} ({unit})" if unit else label)
        scan.labels = dict(scan.labels or {})
        scan.labels["x"] = text
    return scan


def permute_axes(scan, permutation):
    """Reorder a scan's array dimensions, taking its axes with them.

    ``permutation`` is in numpy's ``transpose`` sense: ``(1, 0, 2)`` makes
    the old second dimension the new first. The axis vectors, and the axis
    labels, follow -- which is the part that is easy to get wrong by hand
    and the reason this is offered as one control rather than left to a
    per-beamline transpose.
    """
    slots = _axis_slots(scan.kind)
    permutation = tuple(int(p) for p in permutation)
    if not slots:
        raise ValueError(f"{scan.kind!r} datasets have no permutable axes")
    if sorted(permutation) != list(range(len(slots))):
        raise ValueError(
            f"permutation {permutation} is not an ordering of the "
            f"{len(slots)} axes of a {scan.kind}")
    if permutation == tuple(range(len(slots))):
        return scan

    array = scan.value4d if scan.kind == "spem_4d" else scan.value
    if array is None:
        raise ValueError("this dataset has no array to permute")
    # A lazy array has to be read to be transposed; that is the honest cost
    # of asking for a different order, and it happens once, at load.
    permuted = np.transpose(np.asarray(array), permutation)

    old_axes = [getattr(scan, slot) for slot in slots]
    old_labels = dict(scan.labels or {})
    for new_index, old_index in enumerate(permutation):
        setattr(scan, slots[new_index], old_axes[old_index])
    new_labels = dict(old_labels)
    for new_index, old_index in enumerate(permutation):
        if slots[old_index] in old_labels:
            new_labels[slots[new_index]] = old_labels[slots[old_index]]
        else:
            new_labels.pop(slots[new_index], None)
    scan.labels = new_labels

    if scan.kind == "spem_4d":
        scan.value4d = permuted
        scan.value = permuted
    else:
        scan.value = permuted
    scan.info["loader.permutation"] = str(permutation)
    return scan


def role_is_angle(scan) -> bool:
    """Whether this map's first axis is an emission angle -- i.e. whether a
    k conversion of it means anything. Anything that did not say is assumed
    to be an angle, which is what every file read before this existed was."""
    role = (getattr(scan, "info", {}) or {}).get("axis0.role", "angle")
    _label, _unit, convertible = AXIS_ROLES.get(str(role), AXIS_ROLES["other"])
    return bool(convertible)


def _install_default_loaders():
    """Import the loaders that ship with the program, so that importing this
    module is enough to have them. Kept in a function so the import order is
    explicit and a broken third-party loader cannot stop the rest loading."""
    for module in ("loader.native", "loader.soleil", "loader.cassiopee",
                   "loader.cassiopee_spin"):
        try:
            __import__(module)
        except Exception as exc:                        # noqa: BLE001
            import warnings
            warnings.warn(f"loader module {module} did not load: {exc}")


_install_default_loaders()
