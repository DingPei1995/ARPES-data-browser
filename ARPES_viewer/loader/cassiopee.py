"""
loader/cassiopee.py
===================
SOLEIL CASSIOPEE: Scienta text spectra, and folders of them.

Two shapes of measurement come off this beamline, and they are the same
file format twice over:

* **One spectrum** -- energy against analyser angle -- written by the Scienta
  SES software as a plain text file: ``[Region n]`` headers with the two
  axis scales spelled out number by number, ``[Info n]`` with the
  acquisition settings, ``[Data n]`` with one line per energy point.
* **A folder of them**, numbered, which together are one measurement: a
  Fermi-surface map (the sample's polar angle stepped between spectra) or a
  photon-energy scan (the monochromator stepped instead). Nothing inside any
  one file says which; it is the *set* that says so, by which of the two
  quantities varies across it.

Ported from the lab's ``load_Soleil_Cassiopee_struct.m`` and
``load_Soleil_Cassiopee_folder_struct.m``, with the differences noted where
they occur. The most important ones, in one place:

* **The stepped quantity is worked out and then declared.** The MATLAB
  version returns a struct whose ``x`` is a polar angle or a photon energy
  depending on the branch it took, and nothing downstream can tell which. A
  photon-energy series comes back here with its first axis marked
  ``axis0.role = "photon_energy"``, which is what stops the k conversion
  being offered for it -- converting a photon energy to an in-plane momentum
  is meaningless, and it would produce a perfectly plausible-looking picture.
* **The median filter is off by default.** The MATLAB loader ran a two-point
  median filter along both axes of every spectrum, unconditionally, inside
  the load. That is a processing step: on raw counts it removes single-pixel
  information and biases the statistics, and doing it silently at load time
  means nobody downstream can tell it happened. It is offered here as a
  load-time option (and recorded in ``info`` when used); the despiking in
  the Process panel is the same thing, visible.
* **A file name may contain underscores.** The MATLAB version splits names
  on ``_`` and takes the first piece as the base, which breaks for a folder
  called something like ``Map80eV_Phi18`` -- there is a commented-out
  hand-set ``name_base`` in it that says as much. The series here is found by
  matching the ``_<number>_ROI<n>_`` tail with a regular expression, so the
  base can be anything.

Nothing in here imports Qt.
"""
from __future__ import annotations

import os
import re
import warnings
from dataclasses import dataclass, field

import numpy as np

from loader.nxs_file import NxsScan
from loader.registry import Loader, register

__all__ = ["CassiopeeLoader", "parse_scienta", "read_parameter_file",
           "series_members", "work_function", "ANALYSER_WORK_FUNCTION"]

#: The marker that says "this is a Scienta text spectrum". The axis scales
#: are written out in full, which no other text format here does.
SCIENTA_MARKER = "Dimension 1 scale="

#: Measured analyser work function against photon energy, from
#: ``load_Soleil_Cassiopee_folder_struct.m``. Beamline calibration data, not
#: a formula: it is not monotonic (it dips between 95 and 111 eV), which is
#: why it is interpolated rather than fitted.
ANALYSER_WORK_FUNCTION = {
    "photon_energy_eV": (45.0, 75.0, 90.0, 95.0, 111.0, 120.0, 130.0, 141.0,
                         150.0, 170.0),
    "work_function_eV": (4.2348, 4.5356, 4.7295, 4.7785, 4.6109, 4.68, 4.7803,
                         4.9464, 5.0237, 5.279),
}

#: What the single-spectrum MATLAB loader used when it had no table to go on.
DEFAULT_WORK_FUNCTION = 4.45

#: Lines of the ``_i`` parameter file, and what they are called in ``info``.
#: Matched on the text before the colon, case-insensitively, rather than by
#: character position as the MATLAB version does -- the files are written by
#: hand-rolled formatting and the spacing around the colon varies.
PARAMETER_FIELDS = {
    "x (mm)": "sample_X_mm",
    "y (mm)": "sample_Y_mm",
    "z (mm)": "sample_Z_mm",
    "theta (deg)": "sample_theta_deg",
    "phi (deg)": "sample_phi_deg",
    "tilt (deg)": "sample_tilt_deg",
    "p(mbar)": "vacuum_pressure_mbar",
    "hv mono (ev)": "photon_energy_eV",
    "t_b (k)": "temperature_K",
    "t_a (k)": "temperature_sample_K",
}

#: Polarisation codes, from the commented-out block in the MATLAB loader.
POLARISATIONS = {0: "LV", 1: "LH", 2: "AV", 3: "AH", 4: "CR"}


@dataclass
class ScientaRegion:
    """One ``[Region n]`` of a Scienta text file.

    ``values`` is ``(angle, energy)`` -- the transpose of the file's own
    row-per-energy layout, and the orientation the rest of this program
    calls a ``cut``.
    """
    name: str
    energy: np.ndarray
    angle: np.ndarray
    values: np.ndarray
    energy_label: str = "Kinetic energy (eV)"
    angle_label: str = "Angle (deg)"
    info: dict = field(default_factory=dict)


# --------------------------------------------------------------------------
# Reading one spectrum
# --------------------------------------------------------------------------
def _read_lines(path: str):
    # latin-1 rather than utf-8: these are written by Windows acquisition
    # software and the comment fields occasionally carry a stray byte that
    # utf-8 refuses. Nothing here depends on the text beyond ASCII.
    with open(path, "r", encoding="latin-1", errors="replace") as handle:
        return handle.read().splitlines()


def _sections(lines):
    """``{"[Region 1]": (start, end), ...}`` -- line ranges of each section."""
    starts = [(i, line.strip()) for i, line in enumerate(lines)
              if line.startswith("[") and line.strip().endswith("]")]
    found = {}
    for position, (index, name) in enumerate(starts):
        end = starts[position + 1][0] if position + 1 < len(starts) else len(lines)
        found[name] = (index + 1, end)
    return found


def _fields(lines, span):
    out = {}
    for line in lines[span[0]:span[1]]:
        key, sep, value = line.partition("=")
        if sep:
            out[key.strip()] = value.strip()
    return out


def _numbers(text: str) -> np.ndarray:
    return np.array(text.split(), dtype=float) if text.strip() else np.zeros(0)


def _maybe_numbers(text: str):
    """The numbers in ``text``, or ``None`` if it is not numeric.

    Used for the parameter file, where a field that is written as words --
    a sample name, a comment the operator typed into a numeric slot -- must
    not take down the load of the whole map with it.
    """
    try:
        return _numbers(text)
    except ValueError:
        return None


def _regularise(axis: np.ndarray, name: str, path: str,
                tolerance: float = 0.05) -> np.ndarray:
    """Replace an axis by the straight line through it.

    The scales are printed to fixed precision, so the steps carry rounding
    jitter; the viewers draw an image and want an exactly even axis.
    ``load_Soleil_Cassiopee_struct.m`` does the same (``polyfit`` over the
    point index), and on a real file the correction is around 1e-14 of a
    step -- it is defence, not a fix.

    Unlike the original, a scale that is *genuinely* not linear is reported
    rather than quietly straightened: past ``tolerance`` of one step the
    difference is no longer rounding, and silently linearising it would move
    data onto the wrong energies.
    """
    if axis.size < 2:
        return axis
    index = np.arange(axis.size, dtype=float)
    slope, intercept = np.polyfit(index, axis, 1)
    straight = index * slope + intercept
    step = abs(slope) or 1.0
    worst = float(np.max(np.abs(axis - straight))) / step
    if worst > tolerance:
        warnings.warn(
            f"{os.path.basename(path)}: the {name} scale is not linear "
            f"(off by up to {worst:.2f} of a step); it has been replaced by "
            f"the straight line through it, which may misplace the data.")
    return straight


def _median_filter_2(values: np.ndarray) -> np.ndarray:
    """The two-point median filter the MATLAB loader applied to every
    spectrum, along both axes.

    MATLAB's ``medfilt1`` takes the median of an even-length window as the
    mean of its two middle values, so order 2 is the mean of each point and
    the one before it: a mild low-pass, not a despike. That part is
    reproduced.

    The edge is not. ``medfilt1`` zero-pads, so its first row and first
    column come out at *half* their true value -- a visible dark line down
    two sides of every spectrum, and one that would be averaged into a map
    over and over. The first point is left alone here instead.
    """
    out = np.asarray(values, dtype=float)
    for axis in (0, 1):
        shifted = np.roll(out, 1, axis=axis)
        index = [slice(None)] * out.ndim
        index[axis] = 0
        shifted[tuple(index)] = out[tuple(index)]      # edge: keep the value
        out = 0.5 * (out + shifted)
    return out


def parse_scienta(path: str, median_filter: bool = False,
                  regularise: bool = True) -> "list[ScientaRegion]":
    """Every region of a Scienta text file.

    Usually one, but the format numbers them (``[Region 1]``, ``[Data 1]``,
    ...) and a multi-region file is a legitimate thing to be handed, so all
    of them are read rather than the first.
    """
    lines = _read_lines(path)
    sections = _sections(lines)
    regions = []

    for name, span in sections.items():
        match = re.fullmatch(r"\[Region (\d+)\]", name)
        if not match:
            continue
        number = match.group(1)
        header = _fields(lines, span)
        info_span = sections.get(f"[Info {number}]")
        info = _fields(lines, info_span) if info_span else {}
        data_span = sections.get(f"[Data {number}]")
        if data_span is None:
            continue

        energy = _numbers(header.get("Dimension 1 scale", ""))
        angle = _numbers(header.get("Dimension 2 scale", ""))
        if regularise:
            energy = _regularise(energy, "Dimension 1", path)
            angle = _regularise(angle, "Dimension 2", path)

        block = "\n".join(lines[data_span[0]:data_span[1]]).strip()
        flat = np.array(block.split(), dtype=float) if block else np.zeros(0)
        # One row per energy point: the point's own energy, then one count
        # per angle channel. The leading column is dropped -- it repeats the
        # Dimension 1 scale, which has already been read (and regularised).
        columns = angle.size + 1
        if flat.size % columns:
            raise ValueError(
                f"{os.path.basename(path)}: the [Data {number}] block has "
                f"{flat.size} numbers, which is not a whole number of rows "
                f"of {columns} (1 + Dimension 2 size)")
        table = flat.reshape(-1, columns)
        if table.shape[0] != energy.size:
            raise ValueError(
                f"{os.path.basename(path)}: [Data {number}] has "
                f"{table.shape[0]} rows but Dimension 1 has {energy.size} "
                f"points")
        values = table[:, 1:].T                     # -> (angle, energy)
        if median_filter:
            values = _median_filter_2(values)

        regions.append(ScientaRegion(
            name=header.get("Region Name", f"Region {number}"),
            energy=energy, angle=angle, values=values,
            energy_label=_axis_label(header.get("Dimension 1 name"),
                                     "Kinetic energy (eV)"),
            angle_label=_axis_label(header.get("Dimension 2 name"),
                                    "Angle (deg)"),
            info=_spectrum_info(info, header, path, median_filter)))
    return regions


def _axis_label(raw: str, fallback: str) -> str:
    """"Kinetic Energy [eV]" -> "Kinetic Energy (eV)". The viewers print the
    label verbatim, and the rest of this program writes units in brackets."""
    if not raw:
        return fallback
    return raw.replace("[", "(").replace("]", ")").strip()


def _spectrum_info(info: dict, header: dict, path: str,
                   median_filter: bool) -> dict:
    """The acquisition settings, under names the information panel can show.

    Everything in ``[Info n]`` is kept, prefixed, rather than a chosen few:
    a field the beamline adds later then shows up by itself instead of being
    silently dropped.
    """
    out = {f"scienta.{key}": value for key, value in info.items()}
    out.update({f"scienta.{key}": value for key, value in header.items()
                if not key.startswith("Dimension") or key.endswith("name")})
    for key, name in (("Excitation Energy", "photon_energy_eV"),
                      ("Pass Energy", "pass_energy_eV"),
                      ("Lens Mode", "lens_mode"),
                      ("Location", "beamline"),
                      ("Sample", "sample"),
                      ("Date", "date"), ("Time", "time")):
        if key in info:
            try:
                out[name] = float(info[key])
            except (TypeError, ValueError):
                out[name] = info[key]
    out["cassiopee.median_filter"] = int(bool(median_filter))
    out["cassiopee.source_file"] = os.path.basename(path)
    return out


# --------------------------------------------------------------------------
# The parameter file that sits beside each spectrum
# --------------------------------------------------------------------------
def parameter_path(spectrum_path: str) -> "str | None":
    """The ``_i`` file beside a ``_ROI<n>_`` spectrum, if it is there.

    The pair is ``<base>_<number>_ROI1_.txt`` and ``<base>_<number>_i.txt``:
    the spectrum and the motor positions it was taken at. Only the second
    knows the sample angles and the monochromator's own readback, which is
    what decides whether a folder is a Fermi-surface map or a photon-energy
    scan.
    """
    folder, name = os.path.split(spectrum_path)
    stem, extension = os.path.splitext(name)
    candidate = re.sub(r"ROI\d+_?$", "i", stem)
    if candidate == stem:
        return None
    path = os.path.join(folder, candidate + extension)
    return path if os.path.isfile(path) else None


def _split_parameter_line(line: str):
    """``"Polarisation [0:LV, 1:LH] : 1"`` -> ``("Polarisation [...]", "1")``.

    The separator is the first colon *outside* any brackets. Splitting on the
    plain first colon lands inside the polarisation field's own label, which
    lists its codes as ``[0:LV, 1:LH, 2:AV, 3:AH, 4:CR]``; splitting on the
    last one would eat a value like a timestamp.
    """
    depth = 0
    for position, character in enumerate(line):
        if character == "[":
            depth += 1
        elif character == "]":
            depth = max(0, depth - 1)
        elif character == ":" and depth == 0:
            return line[:position].strip(), line[position + 1:].strip()
    return None, None


def read_parameter_file(path: str) -> dict:
    """Motor positions and beamline readbacks from an ``_i`` file.

    Split on the colon rather than at fixed character positions (which is
    what the MATLAB version does): the files are hand-formatted and the
    spacing varies, and a field that moved by a character would otherwise
    come back as a silent NaN.
    """
    if not path or not os.path.isfile(path):
        return {}
    out = {}
    for line in _read_lines(path):
        key, value = _split_parameter_line(line)
        if key is None:
            continue
        lowered = key.strip().lower()
        name = PARAMETER_FIELDS.get(lowered)
        if name:
            numbers = _maybe_numbers(value)
            if numbers is None:
                out[name] = value            # a field written as text
            elif numbers.size == 1:
                out[name] = float(numbers[0])
            elif numbers.size:
                out[name] = numbers
        elif lowered.startswith("polarisation"):
            numbers = _maybe_numbers(value)
            if numbers is not None and numbers.size:
                out["polarisation"] = POLARISATIONS.get(int(numbers[-1]),
                                                        str(int(numbers[-1])))
    return out


# --------------------------------------------------------------------------
# The work-function calibration
# --------------------------------------------------------------------------
def work_function(photon_energy_eV: float) -> float:
    """The analyser work function at this photon energy.

    Interpolated through :data:`ANALYSER_WORK_FUNCTION` with a cubic spline,
    matching ``interp1(..., 'spline', 'extrap')`` in the MATLAB loader.

    Outside the calibrated 45-170 eV the spline is *extrapolated*, which the
    MATLAB version also does silently. A cubic through a non-monotonic table
    leaves its range fast, so it is warned about here: a wrong work function
    shifts the binding-energy axis of a whole photon-energy scan, and the
    picture looks fine either way.
    """
    hv = np.asarray(ANALYSER_WORK_FUNCTION["photon_energy_eV"], dtype=float)
    phi = np.asarray(ANALYSER_WORK_FUNCTION["work_function_eV"], dtype=float)
    value = float(photon_energy_eV)
    if value < hv[0] or value > hv[-1]:
        warnings.warn(
            f"photon energy {value:g} eV is outside the calibrated "
            f"{hv[0]:g}-{hv[-1]:g} eV range of the CASSIOPEE work-function "
            f"table; the value is extrapolated and may be well off.")
    try:
        from scipy.interpolate import CubicSpline
        spline = CubicSpline(hv, phi, bc_type="not-a-knot", extrapolate=True)
        return float(spline(value))
    except ImportError:                                     # pragma: no cover
        return float(np.interp(value, hv, phi))


# --------------------------------------------------------------------------
# A folder of numbered spectra
# --------------------------------------------------------------------------
#: ``<base>_<number>_ROI<n>_`` -- anchored at the end so the base may contain
#: underscores, which is where the MATLAB version's ``split('_')`` gives up.
SERIES_PATTERN = re.compile(r"^(?P<base>.+)_(?P<index>\d+)_ROI(?P<roi>\d+)_?$")


def series_members(path: str):
    """Every spectrum of the numbered series ``path`` belongs to, in order.

    Returns ``[(index, spectrum path), ...]``, or ``[]`` if this file is not
    part of one. Membership is by *name*: same folder, same base, same ROI
    number, differing only in the index.
    """
    folder, name = os.path.split(os.path.abspath(path))
    stem, extension = os.path.splitext(name)
    match = SERIES_PATTERN.match(stem)
    if not match:
        return []
    base, roi = match.group("base"), match.group("roi")

    found = []
    try:
        names = os.listdir(folder)
    except OSError:
        return []
    for other in names:
        other_stem, other_extension = os.path.splitext(other)
        if other_extension.lower() != extension.lower():
            continue
        other_match = SERIES_PATTERN.match(other_stem)
        if not other_match:
            continue
        if other_match.group("base") == base and other_match.group("roi") == roi:
            found.append((int(other_match.group("index")),
                          os.path.join(folder, other)))
    return sorted(found)


def _member_scan(path: str, median_filter: bool):
    """One member of a series: its region, plus the parameters beside it."""
    regions = parse_scienta(path, median_filter=median_filter)
    if not regions:
        raise ValueError(f"{os.path.basename(path)}: no spectrum in it")
    region = regions[0]
    parameters = read_parameter_file(parameter_path(path))
    return region, parameters


def _stepped_axis(thetas, photon_energies, count):
    """Which quantity was stepped across the series, and so what the first
    axis of the assembled map *is*.

    The file format does not say; the set does. This is the decision the
    MATLAB loader makes with two ``if`` statements and then forgets -- here
    it comes back as a role the rest of the program can act on (the k
    conversion refuses a photon-energy axis, for one).
    """
    def varies(values):
        finite = [v for v in values if v is not None and np.isfinite(v)]
        return len(finite) == len(values) and len(set(np.round(finite, 6))) > 1

    if varies(thetas):
        return (np.asarray(thetas, dtype=float), "angle",
                "Polar angle theta (deg)", "theta")
    if varies(photon_energies):
        return (np.asarray(photon_energies, dtype=float), "photon_energy",
                "Photon energy (eV)", "hv")
    return (np.arange(1, count + 1, dtype=float), "other", "Index (cut)",
            "index")


def load_series(path: str, median_filter: bool = False, progress=None) -> NxsScan:
    """Assemble the numbered series ``path`` belongs to into one map.

    The result is an ordinary ``map``: first axis whatever was stepped,
    then the analyser angle, then energy -- the same shape a deflector map
    from another beamline has, so every viewer, cut, fit and figure works on
    it without knowing where it came from.

    ``progress`` is an optional ``callable(done, total, label)``; reading a
    hundred 6 MB text files is slow enough to be worth reporting.
    """
    members = series_members(path)
    if not members:
        raise ValueError(
            f"{os.path.basename(path)} is not part of a numbered "
            f"<base>_<n>_ROI<n>_ series")

    frames, thetas, photon_energies, first = [], [], [], None
    for position, (index, member) in enumerate(members):
        if progress is not None:
            progress(position, len(members), os.path.basename(member))
        region, parameters = _member_scan(member, median_filter)
        if first is None:
            first = (region, parameters)
        elif region.values.shape != first[0].values.shape:
            raise ValueError(
                f"{os.path.basename(member)} is {region.values.shape}, but "
                f"the first spectrum of the series is "
                f"{first[0].values.shape}; they cannot be one map")
        frames.append(region.values)
        thetas.append(parameters.get("sample_theta_deg"))
        photon_energies.append(parameters.get(
            "photon_energy_eV", region.info.get("photon_energy_eV")))

    region, parameters = first
    stepped, role, label, kind = _stepped_axis(thetas, photon_energies,
                                               len(members))
    cube = np.stack(frames, axis=0)              # (step, angle, energy)

    scan = NxsScan(kind="map", filename_prefix="")
    scan.x, scan.k = stepped, region.angle
    scan.value = cube
    energy_label = region.energy_label

    if kind == "hv":
        # A photon-energy scan is only comparable across its members once the
        # energy axis is referred to the Fermi level: each spectrum was taken
        # at a different hv, so the same kinetic energy is a different state.
        # E - E_F = E_kin - hv + phi, as the MATLAB loader does.
        reference_hv = float(photon_energies[0])
        phi = work_function(reference_hv)
        scan.z = region.energy - reference_hv + phi
        energy_label = "E - E_F (eV)"
        scan.info["cassiopee.work_function_eV"] = phi
        scan.info["cassiopee.reference_photon_energy_eV"] = reference_hv
    else:
        scan.z = region.energy
        if kind == "theta":
            reference_hv = photon_energies[0]
            if reference_hv is not None and np.isfinite(reference_hv):
                scan.info["cassiopee.work_function_eV"] = \
                    work_function(float(reference_hv))

    scan.labels = {"x": label, "k": region.angle_label, "z": energy_label}
    scan.info.update(region.info)
    scan.info.update({f"cassiopee.{key}": value
                      for key, value in parameters.items()})
    scan.info["axis0.role"] = role
    scan.info["cassiopee.series"] = kind
    scan.info["cassiopee.members"] = len(members)
    if kind == "theta":
        scan.info["cassiopee.theta_range_deg"] = \
            f"{stepped.min():g} to {stepped.max():g}"
    elif kind == "hv":
        scan.info["cassiopee.photon_energy_range_eV"] = \
            f"{stepped.min():g} to {stepped.max():g}"
    return scan


def load_cut(path: str, region_index: int = 0,
             median_filter: bool = False) -> NxsScan:
    """One spectrum as a ``cut``: analyser angle against energy."""
    regions = parse_scienta(path, median_filter=median_filter)
    if not regions:
        raise ValueError(f"{os.path.basename(path)}: no spectrum in it")
    region = regions[min(region_index, len(regions) - 1)]

    scan = NxsScan(kind="cut", filename_prefix="")
    scan.x, scan.y = region.angle, region.energy
    scan.value = region.values
    scan.labels = {"x": region.angle_label, "y": region.energy_label}
    scan.info.update(region.info)
    parameters = read_parameter_file(parameter_path(path))
    scan.info.update({f"cassiopee.{key}": value
                      for key, value in parameters.items()})
    photon_energy = parameters.get("photon_energy_eV",
                                   region.info.get("photon_energy_eV"))
    if photon_energy is not None:
        # The monochromator's own readback, where there is one, beats the
        # value written into the spectrum header -- which is the setpoint.
        scan.info["photon_energy_eV"] = float(photon_energy)
    scan.info.setdefault("cassiopee.work_function_eV", DEFAULT_WORK_FUNCTION)
    return scan


# --------------------------------------------------------------------------
# The loader itself
# --------------------------------------------------------------------------
class CassiopeeLoader(Loader):
    name = "SOLEIL CASSIOPEE"
    description = ("Scienta text spectra from the CASSIOPEE beamline at "
                   "SOLEIL, and folders of numbered spectra assembled into "
                   "a Fermi-surface map or a photon-energy scan.")
    patterns = ("*.txt",)
    priority = 20

    def can_open(self, path: str) -> bool:
        """Recognised by content: a Scienta axis scale in the first few lines.

        The parameter files that sit beside the spectra are also ``.txt`` and
        are deliberately *not* recognised -- they hold motor positions, not a
        spectrum, and listing them would put rows in the browser that cannot
        be opened.

        Binary is rejected first. Every loader gets asked about every file,
        this one runs before the HDF5 readers, and ``latin-1`` decodes
        arbitrary bytes without complaint -- so without the check a compressed
        chunk that happened to contain the marker's bytes would be claimed as
        a text spectrum. A NUL byte in the first block is the usual quick
        test, and no text file has one.
        """
        try:
            with open(path, "rb") as handle:
                head = handle.read(4096)
        except OSError:
            return False
        if b"\x00" in head:
            return False
        return SCIENTA_MARKER.encode("latin-1") in head

    def list_entries(self, path: str) -> list:
        """One entry per region, plus the whole series if this file is in one.

        Offering both is the point: the same file is a spectrum in its own
        right and a slice of a map, and which one is wanted depends on what
        is being looked at.
        """
        entries = []
        try:
            regions = parse_scienta(path)
        except Exception:
            return []
        for index, region in enumerate(regions):
            entries.append({
                "entry": f"region {index + 1}" if len(regions) > 1 else "cut",
                "kind": "Cut",
                "name": region.name or os.path.splitext(os.path.basename(path))[0],
                "title": region.name,
                "start_time": region.info.get("date"),
            })

        members = series_members(path)
        if len(members) > 1:
            base = SERIES_PATTERN.match(
                os.path.splitext(os.path.basename(path))[0]).group("base")
            entries.append({
                "entry": "series",
                "kind": "Map",
                "name": f"{base} ({len(members)} cuts)",
                "title": f"{len(members)} spectra assembled",
                "start_time": None,
            })
        return entries

    def __init__(self, median_filter: bool = False):
        #: Off, for the reason in the module docstring: filtering raw counts
        #: invisibly at load time is a processing step nobody downstream can
        #: see. The Process panel's despiking is the same operation, visible
        #: and undoable. A script that wants the MATLAB behaviour can call
        #: :func:`load_cut` or :func:`load_series` with ``median_filter=True``,
        #: or register an instance of this class made with it on.
        self.median_filter = bool(median_filter)

    def load(self, path: str, entry: str = None, progress=None):
        if entry == "series":
            return load_series(path, median_filter=self.median_filter,
                               progress=progress)
        index = 0
        if entry and entry.startswith("region "):
            index = int(entry.split()[1]) - 1
        return load_cut(path, region_index=index,
                        median_filter=self.median_filter)


register(CassiopeeLoader())
