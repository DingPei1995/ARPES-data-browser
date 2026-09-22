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
    """``{"[Region 1]": (start, end), ...}`` -- line ranges of each section.

    ``strip()`` on both ends because the older writer puts a tab after the
    header: the line is ``"[Region 1]\\t "``, not ``"[Region 1]"``.
    """
    starts = [(i, line.strip()) for i, line in enumerate(lines)
              if line.lstrip().startswith("[") and line.strip().endswith("]")]
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


def _is_numeric_row(line: str, least: int = 2) -> bool:
    """Is this line a row of the data block?

    A row is ``energy`` followed by one count per angle channel, all
    whitespace- or tab-separated. Everything else in these files is either a
    ``key=value`` line, a ``[Section]`` header, or blank -- none of which
    parse as several numbers, so "does it parse" is a sound test and does
    not depend on knowing where the block starts.
    """
    parts = line.split()
    if len(parts) < least:
        return False
    for part in parts:
        try:
            float(part)
        except ValueError:
            return False
    return True


def _data_block(lines, sections, number: str):
    """The line range holding region ``number``'s counts.

    Two writers, two layouts, and the difference is not announced anywhere
    in the file:

    * **SES 1.3.1** puts the counts in a ``[Data n]`` section of their own.
    * **SES 1.2.5** has no such section. The counts simply follow the
      ``[Run Mode Information n]`` block, after a stray ``inputA=`` line and
      a blank one.

    ``load_Soleil_Cassiopee_struct.m`` handles the second by counting three
    lines on from ``[Run Mode Information`` and trusting the offset. That
    works for the files it was written against and breaks silently on the
    first file with one more or one fewer trailing field -- it would read
    counts as an axis, or an axis as counts, and still produce a picture.
    Here the block is found by what it *is*: the run of lines that parse as
    rows of numbers. That recognises both layouts with one rule.
    """
    explicit = sections.get(f"[Data {number}]")
    if explicit is not None:
        return explicit

    after = sections.get(f"[Run Mode Information {number}]")
    if after is None:
        after = sections.get(f"[Info {number}]") or sections.get(f"[Region {number}]")
    if after is None:
        return None

    start = after[0]
    while start < len(lines) and not _is_numeric_row(lines[start]):
        # Only blank lines and key=value leftovers may be skipped; a section
        # header means the region ended without any counts.
        if lines[start].strip().startswith("["):
            return None
        start += 1
    if start >= len(lines):
        return None
    end = start
    while end < len(lines) and _is_numeric_row(lines[end]):
        end += 1
    return (start, end)


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


def _decimals_needed(value: float, most: int = 9) -> int:
    """The fewest decimal places that write ``value`` exactly."""
    for decimals in range(most + 1):
        scaled = value * (10.0 ** decimals)
        if abs(scaled - round(scaled)) < 1e-6 * max(1.0, abs(scaled)):
            return decimals
    return most


def _print_quantum(values: np.ndarray) -> float:
    """The coarsest grid the numbers could have been *printed* on.

    A value written out to a finite precision cannot be further than half a
    printing step from where it really was, so this is how much of an axis's
    raggedness is the file's formatting rather than the measurement.

    The formatting is *significant figures*, not decimal places, which
    matters more than it sounds. A real energy scale here crosses 100 eV,
    and the writer -- printing five significant figures -- switches from
    ``99.998`` to ``100.01`` as it does: three decimals below the boundary,
    two above it. So the rounding grid gets ten times coarser part-way
    through one axis. Counting decimals instead of significant figures takes
    the fine half of the axis for the whole of it, and then reports the
    coarse half as a non-linear scale -- which is exactly what a real file
    from this beamline did.

    Taking the *widest* significant-figure count in the array, rather than
    the narrowest, is deliberate: a value that lands on a round number
    (``100.00``) looks like it was printed to fewer figures than it was, and
    believing that would make the allowance enormous. At least one value in
    a real scale shows the full precision.
    """
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite) & (finite != 0.0)]
    if finite.size == 0:
        return 0.0
    figures = 0
    for value in finite:
        integer_digits = int(np.floor(np.log10(abs(value)))) + 1
        figures = max(figures, _decimals_needed(float(value)) + integer_digits)
    if figures <= 0:
        return 0.0
    largest = int(np.floor(np.log10(np.max(np.abs(finite)))))
    return 10.0 ** (largest - figures + 1)


def _regularise(axis: np.ndarray, name: str, path: str,
                tolerance: float = 0.05) -> np.ndarray:
    """Replace an axis by the straight line through it.

    The scales are printed to fixed precision, so the steps carry rounding
    jitter; the viewers draw an image and want an exactly even axis.
    ``load_Soleil_Cassiopee_struct.m`` does the same (``polyfit`` over the
    point index).

    Unlike the original, a scale that is *genuinely* not linear is reported
    rather than quietly straightened -- silently linearising a curved scale
    moves data onto the wrong energies.

    What counts as "genuinely" has to allow for how the numbers were
    printed, which is the part that is easy to get wrong. A real energy
    scale here steps by 0.004 eV and is written to three decimals, so
    rounding alone displaces a point by up to 0.0005 eV -- an eighth of a
    step. A flat "5% of a step" rule calls every one of those files
    non-linear, which is a warning on correct data, and a warning that
    always fires is one nobody reads. The allowance is therefore the larger
    of a fraction of a step and the printing grid itself.
    """
    if axis.size < 2:
        return axis
    index = np.arange(axis.size, dtype=float)
    slope, intercept = np.polyfit(index, axis, 1)
    straight = index * slope + intercept
    step = abs(slope) or 1.0
    worst = float(np.max(np.abs(axis - straight)))
    allowed = max(tolerance * step, _print_quantum(axis))
    if worst > allowed:
        warnings.warn(
            f"{os.path.basename(path)}: the {name} scale is not linear "
            f"(off by up to {worst / step:.2f} of a step); it has been "
            f"replaced by the straight line through it, which may misplace "
            f"the data.")
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
        data_span = _data_block(lines, sections, number)
        if data_span is None:
            continue

        # The point counts come from the scales, never from the
        # "Dimension n size" fields. Those are unreliable: a real file from
        # the older writer says ``Dimension 1 size=264.64`` for a scale of
        # 629 points -- not merely wrong but not even a whole number.
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
        elif lowered.startswith("polarisation") and value:
            # Written either as the code its own label documents
            # (``Polarisation [0:LV, 1:LH, ...] : 1``) or, by the newer
            # station software, as the word itself (``Polarisation : LH``).
            numbers = _maybe_numbers(value)
            if numbers is not None and numbers.size:
                out["polarisation"] = POLARISATIONS.get(int(numbers[-1]),
                                                        str(int(numbers[-1])))
            else:
                out["polarisation"] = value
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


#: What the browser calls a folder, by what was stepped across it.
KIND_FOR_SERIES = {"theta": "Map", "hv": "kz map", "index": "Map"}

SERIES_TITLES = {
    "theta": "{n} spectra, polar angle stepped",
    "hv": "{n} spectra, photon energy stepped",
    "index": "{n} spectra, nothing stepped",
}


def series_kind(path: str) -> str:
    """What the folder ``path`` belongs to is: ``"theta"``, ``"hv"`` or
    ``"index"`` -- without reading a single spectrum.

    The browser has to label a folder before anyone opens it, and at that
    moment reading the spectra is out of the question: they are a megabyte
    each and there are sixty of them. It does not have to. Which quantity
    was stepped is recorded in the ``_i`` parameter files, which are a few
    hundred bytes, so the whole folder can be characterised for about the
    cost of one spectrum.

    Falls back to ``"index"`` if the parameter files are missing, which is
    also what an unstepped folder is -- both mean "nothing here says this is
    more than a pile of spectra".
    """
    members = series_members(path)
    if len(members) < 2:
        return "index"
    thetas, photon_energies = [], []
    for _index, member in members:
        parameters = read_parameter_file(parameter_path(member))
        thetas.append(parameters.get("sample_theta_deg"))
        photon_energies.append(parameters.get("photon_energy_eV"))
    return _stepped_axis(thetas, photon_energies, len(members))[3]


def _member_scan(path: str, median_filter: bool):
    """One member of a series: its region, plus the parameters beside it."""
    regions = parse_scienta(path, median_filter=median_filter)
    if not regions:
        raise ValueError(f"{os.path.basename(path)}: no spectrum in it")
    region = regions[0]
    parameters = read_parameter_file(parameter_path(path))
    return region, parameters


def _even_steps(values: np.ndarray) -> np.ndarray:
    """The evenly stepped axis the series was *meant* to be measured on.

    A scan steps one motor by a fixed amount, so the axis is even by
    intention. What the parameter files record is where the motor actually
    stopped, which carries the mechanism's own error: a real 61-point polar
    scan here asks for 0.5 deg a step and reads back 0.504, 0.504, 0.495,
    0.504 ... That jitter is the instrument, not the sample.

    It is straightened rather than kept, for two reasons. The viewers draw a
    cube as an image, which places rows at even intervals whatever the axis
    says -- so an uneven axis is not honoured, it is merely mislabelled. And
    a k conversion differentiating along a jittering angle turns a 0.01 deg
    readback error into a visible ripple.

    The readbacks are kept in ``info`` (``cassiopee.axis0_measured``), so
    nothing is thrown away and the spread can be checked.
    """
    values = np.asarray(values, dtype=float)
    if values.size < 2:
        return values
    index = np.arange(values.size, dtype=float)
    slope, intercept = np.polyfit(index, values, 1)
    return index * slope + intercept


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
        return (_even_steps(thetas), "angle",
                "Polar angle theta (deg)", "theta")
    if varies(photon_energies):
        return (_even_steps(photon_energies), "photon_energy",
                "Photon energy (eV)", "hv")
    return (np.arange(1, count + 1, dtype=float), "other", "Index (cut)",
            "index")


def load_series(path: str, median_filter: bool = False, progress=None,
                energy_reference: str = "common") -> NxsScan:
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

    frames, energies, thetas, photon_energies, first = [], [], [], [], None
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
        energies.append(region.energy)
        thetas.append(parameters.get("sample_theta_deg"))
        photon_energies.append(parameters.get(
            "photon_energy_eV", region.info.get("photon_energy_eV")))

    region, parameters = first
    stepped, role, label, series = _stepped_axis(thetas, photon_energies,
                                                 len(members))
    cube = np.stack(frames, axis=0)              # (step, angle, energy)
    energy_label = region.energy_label
    notes = {}

    if series == "hv":
        cube, energy_axis, energy_label, notes = _refer_to_fermi_level(
            cube, energies, photon_energies, energy_reference)
        kind = "kz_map"
    else:
        energy_axis = region.energy
        kind = "map"
        if series == "theta":
            reference_hv = photon_energies[0]
            if reference_hv is not None and np.isfinite(reference_hv):
                notes["cassiopee.work_function_eV"] = \
                    work_function(float(reference_hv))

    scan = NxsScan(kind=kind, filename_prefix="")
    scan.x, scan.k, scan.z = stepped, region.angle, energy_axis
    scan.value = cube
    scan.labels = {"x": label, "k": region.angle_label, "z": energy_label}
    scan.info.update(region.info)
    scan.info.update({f"cassiopee.{key}": value
                      for key, value in parameters.items()})
    scan.info.update(notes)
    scan.info["axis0.role"] = role
    scan.info["cassiopee.series"] = series
    scan.info["cassiopee.members"] = len(members)

    # What the motors actually read back, beside the even axis that replaced
    # them (see _even_steps). Keeping both means the straightening can be
    # checked, and undone, by whoever needs to.
    measured = thetas if series == "theta" else (
        photon_energies if series == "hv" else None)
    if measured is not None and all(
            value is not None and np.isfinite(value) for value in measured):
        measured = np.asarray(measured, dtype=float)
        scan.info["cassiopee.axis0_measured"] = measured
        wobble = float(np.max(np.abs(measured - stepped)))
        scan.info["cassiopee.axis0_readback_spread"] = f"{wobble:.4g}"

    if series == "theta":
        scan.info["cassiopee.theta_range_deg"] = \
            f"{stepped.min():g} to {stepped.max():g}"
    elif series == "hv":
        scan.info["cassiopee.photon_energy_range_eV"] = \
            f"{stepped.min():g} to {stepped.max():g}"
    return scan


#: How a photon-energy scan's members are put on one energy axis. A cube
#: has one energy axis and the members were each measured at a different
#: photon energy, so this choice has to be made and there is no option that
#: is right without an assumption.
#:
#: ``"common"``      **the default.** Assume every member covers the same
#:                   binding-energy window, and take the first member's
#:                   referred axis for all of them. Nothing is interpolated
#:                   and nothing is trimmed. This is what the MATLAB loader
#:                   does, and it is right whenever the operator set each
#:                   spectrum's window from the Fermi level they measured --
#:                   which is how an hv scan is normally taken, since the
#:                   point is to have E_F in frame at every photon energy.
#: ``"per_member"``  Refer each member by its own hv and phi(hv), then
#:                   resample onto the range they all share. Trusts the
#:                   nominal photon energies and the work-function table to
#:                   place E_F. Correct when they are; when they are not, it
#:                   trims the cube by however wrong they are -- and the part
#:                   it trims off the top is E_F.
#: ``"kinetic"``     No referencing at all: the first member's kinetic-energy
#:                   axis, for calibrating by hand afterwards.
#:
#: Whichever is chosen, the two are compared and the disagreement reported,
#: because that disagreement is the measurement telling you its photon-energy
#: calibration needs a Fermi edge fitted per hv.
ENERGY_REFERENCES = ("common", "per_member", "kinetic")


def _refer_to_fermi_level(cube, energies, photon_energies,
                          reference="common"):
    """Put every member of a photon-energy scan on one binding-energy axis.

    Each spectrum of an hv scan is taken at a different photon energy, so
    the same *kinetic* energy is a different state in each one. Only
    ``E - E_F = E_kin - hv + phi(hv)`` is comparable across the set, and a
    cube needs one energy axis, so the members have to be brought onto a
    common one.

    ``load_Soleil_Cassiopee_folder_struct.m`` does this by taking the
    **first** member's axis, shifting it by the **first** member's photon
    energy, and using that for the whole cube. That is only correct if every
    member was measured over the same kinetic-energy window. On the real
    scans from this beamline they are not: across a 40-120 eV scan the
    window drifts by about 1 eV -- 250 energy channels -- because the
    operator sets it from each photon energy in turn. Using member one's
    axis for member forty puts its data a whole electronvolt from where it
    belongs, and nothing about the picture says so.

    So each member is referred by its own hv, and then resampled onto the
    common axis. The resampling is linear and along energy only. The axis
    spans the *overlap* of the members rather than their union: outside it
    some members have no data at all, and a cube half full of holes is
    worse than a slightly shorter one.

    What this cannot fix is that ``phi(hv)`` comes from a calibration table,
    not from these measurements. If the members still do not line up after
    referring them, the scan needs a Fermi edge fitted per photon energy --
    which is a real step, not a detail, so the residual is measured here and
    reported rather than left to be discovered in the picture.
    """
    if reference not in ENERGY_REFERENCES:
        raise ValueError(
            f"energy_reference must be one of {ENERGY_REFERENCES}, "
            f"not {reference!r}")
    if reference == "kinetic":
        return (cube, np.asarray(energies[0], dtype=float),
                "Kinetic energy (eV)",
                {"cassiopee.energy_reference":
                     "kinetic energy, first member's axis (not referred to "
                     "the Fermi level)"})

    hv = np.asarray(photon_energies, dtype=float)
    if not np.all(np.isfinite(hv)):
        # Nothing to refer them by. Fall back to the raw kinetic axis of the
        # first member and say so, rather than inventing a Fermi level.
        warnings.warn(
            "this photon-energy scan has members with no photon energy "
            "recorded; the energy axis is left as kinetic energy.")
        return (cube, np.asarray(energies[0], dtype=float),
                "Kinetic energy (eV)",
                {"cassiopee.energy_reference": "kinetic (hv unknown)"})

    with warnings.catch_warnings():
        # One warning for the scan, not one per member, if the table is
        # being extrapolated.
        warnings.simplefilter("ignore", UserWarning)
        phi = np.array([work_function(value) for value in hv])
    outside = (hv < ANALYSER_WORK_FUNCTION["photon_energy_eV"][0]) | \
              (hv > ANALYSER_WORK_FUNCTION["photon_energy_eV"][-1])
    if outside.any():
        warnings.warn(
            f"{int(outside.sum())} of {hv.size} photon energies "
            f"({hv[outside].min():g}-{hv[outside].max():g} eV) are outside "
            f"the calibrated "
            f"{ANALYSER_WORK_FUNCTION['photon_energy_eV'][0]:g}-"
            f"{ANALYSER_WORK_FUNCTION['photon_energy_eV'][-1]:g} eV "
            f"work-function table; their binding energies are extrapolated.")

    binding = [np.asarray(e, dtype=float) - h + p
               for e, h, p in zip(energies, hv, phi)]
    spread = float(max(float(b[0]) for b in binding)
                   - min(float(b[0]) for b in binding))
    step = float(np.median(np.abs(np.diff(binding[0])))) if binding[0].size > 1 \
        else 0.0

    if reference == "common":
        # Stack them as they are, on the first member's axis. Nothing is
        # interpolated, nothing is trimmed, and every count stays where the
        # file put it -- which is what makes this the honest default for a
        # cube that is going to be calibrated properly later.
        #
        # It is an assumption, and worth naming: that each member's window
        # was set from the Fermi level measured at that photon energy, so
        # they already share a binding-energy scale. That is how an hv scan
        # is normally taken. Where it does not hold, the fix is a Fermi edge
        # fitted per photon energy -- a real step with its own tool, not
        # something to guess at during a load. The numbers that step needs
        # are recorded below rather than warned about, since a warning on
        # every load of every hv scan is noise.
        return cube, binding[0], "E - E_F (eV)", {
            "cassiopee.energy_reference":
                "E - E_F, first member's axis for all (members stacked as "
                "measured; calibrate with a Fermi edge per photon energy)",
            "cassiopee.work_function_eV": float(phi[0]),
            "cassiopee.member_window_spread_eV": float(spread),
            # Enough to rebuild every member's own kinetic-energy axis
            # exactly -- they share a step and a length -- so a later
            # calibration can place each one for itself without re-reading
            # sixty megabytes of text.
            "cassiopee.member_photon_energy_eV": hv,
            "cassiopee.member_work_function_eV": phi,
            "cassiopee.member_kinetic_start_eV": np.array(
                [float(e[0]) for e in energies]),
            "cassiopee.member_kinetic_step_eV": float(step),
        }

    low = max(float(b[0]) for b in binding)
    high = min(float(b[-1]) for b in binding)
    if not (high > low):
        raise ValueError(
            "the members of this photon-energy scan have no binding-energy "
            "range in common once referred to the Fermi level, so they "
            "cannot be one cube. Check that the photon energies in the "
            "parameter files are right, or load it with "
            "energy_reference='common'.")

    points = max(int(b.size) for b in binding)
    axis = np.linspace(low, high, points)

    out = np.empty((cube.shape[0], cube.shape[1], points), dtype=float)
    for index, b in enumerate(binding):
        if np.array_equal(b, axis):
            out[index] = cube[index]
            continue
        # np.interp needs an increasing sample axis; a scale written
        # high-to-low is legal and does happen.
        order = slice(None) if b[0] <= b[-1] else slice(None, None, -1)
        source = b[order]
        for channel in range(cube.shape[1]):
            out[index, channel] = np.interp(axis, source,
                                            cube[index, channel][order])

    widest = (min(float(b[0]) for b in binding),
              max(float(b[-1]) for b in binding))
    trimmed = (low - widest[0]) + (widest[1] - high)
    notes = {
        "cassiopee.energy_reference": "E - E_F, per member (hv and phi(hv))",
        "cassiopee.work_function_eV": f"{phi.min():.4f} to {phi.max():.4f}",
        "cassiopee.binding_energy_range_eV": f"{low:.4f} to {high:.4f}",
        "cassiopee.member_window_spread_eV": f"{spread:.4f}",
        "cassiopee.energy_trimmed_eV": f"{trimmed:.4f}",
    }
    if step and spread > 5 * step:
        # The members should land on top of each other once referred. That
        # they do not means the nominal photon energies and the tabulated
        # work function together do not describe this scan -- and since the
        # common range is cut down by exactly that disagreement, the trim is
        # being driven by the calibration error rather than by the data.
        # Worth saying plainly, because the part trimmed off the top is the
        # Fermi level, which is usually the whole point of the measurement.
        notes["cassiopee.alignment_warning"] = (
            "members disagree by more than a few energy steps; fit a Fermi "
            "edge per photon energy to calibrate")
        warnings.warn(
            f"after referring each member to E_F by its own photon energy, "
            f"their energy windows still differ by {spread:.3f} eV "
            f"({spread / step:.0f} energy steps), and keeping only the range "
            f"they all share cost {trimmed:.3f} eV. The nominal photon "
            f"energies and the tabulated work function do not describe this "
            f"scan between them. Fit a Fermi edge per photon energy to "
            f"calibrate it, or load it with energy_reference='common' to "
            f"keep the full window on the assumption that every member "
            f"covers the same binding energies.")
    return out, axis, "E - E_F (eV)", notes


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
        """The whole folder first, then this one spectrum.

        A numbered folder here is *one measurement* -- the map is the point
        and the individual spectra are its rows -- so the assembled series is
        offered first and is what opening the folder gives you. The single
        spectrum stays available below it, because checking one cut of a map
        is a normal thing to want.

        The series entry carries a ``path`` of its own: the folder's first
        member, whichever member was actually clicked. That is what makes
        selecting five files out of a folder add *one* map rather than five
        copies of it -- the entry identifies the measurement, not the file
        that happened to be picked.
        """
        entries = []
        members = series_members(path)
        if len(members) > 1:
            base = SERIES_PATTERN.match(
                os.path.splitext(os.path.basename(path))[0]).group("base")
            stepped = series_kind(path)
            entries.append({
                "entry": "series",
                "path": members[0][1],
                "kind": KIND_FOR_SERIES.get(stepped, "Map"),
                "name": f"{base} ({len(members)} cuts)",
                "title": SERIES_TITLES.get(stepped, "").format(n=len(members)),
                "start_time": None,
            })

        try:
            regions = parse_scienta(path)
        except Exception:
            return entries
        for index, region in enumerate(regions):
            entries.append({
                "entry": f"region {index + 1}" if len(regions) > 1 else "cut",
                "kind": "Cut",
                "name": region.name or os.path.splitext(os.path.basename(path))[0],
                "title": region.name,
                "start_time": region.info.get("date"),
            })
        return entries

    def __init__(self, median_filter: bool = False,
                 energy_reference: str = "common"):
        #: Off, for the reason in the module docstring: filtering raw counts
        #: invisibly at load time is a processing step nobody downstream can
        #: see. The Process panel's despiking is the same operation, visible
        #: and undoable. A script that wants the MATLAB behaviour can call
        #: :func:`load_cut` or :func:`load_series` with ``median_filter=True``,
        #: or register an instance of this class made with it on.
        self.median_filter = bool(median_filter)
        #: How a photon-energy scan's members are put on one energy axis;
        #: see :data:`ENERGY_REFERENCES`.
        self.energy_reference = energy_reference

    def load(self, path: str, entry: str = None, progress=None):
        if entry == "series":
            return load_series(path, median_filter=self.median_filter,
                               progress=progress,
                               energy_reference=self.energy_reference)
        index = 0
        if entry and entry.startswith("region "):
            index = int(entry.split()[1]) - 1
        return load_cut(path, region_index=index,
                        median_filter=self.median_filter)


register(CassiopeeLoader())
