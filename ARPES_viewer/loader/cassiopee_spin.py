"""
loader/cassiopee_spin.py
========================
SOLEIL CASSIOPEE, spin-resolved end station: MBS A-1 analyser files.

The MBS acquisition software (``MBS A1Soft``) writes two formats, and both
come off this end station:

``.krx``  binary. What the analyser actually saved.
``.txt``  a text export of the same thing.

and three kinds of measurement, which the header -- not the file name --
tells apart:

=========  ==============================================  ==============
Kind       What is in the file                             Opens as
=========  ==============================================  ==============
Cut        one image, angle × energy (main detector)       ``cut``
Map        one image per deflector step                    ``map``
Spin       one 1 × N spectrum per spin channel             ``spin_edc``
=========  ==============================================  ==============

The ``.krx`` layout
-------------------
Worked out from the files themselves, checked against the ``.txt`` export of
the same cut (identical to the count) and against the header's own sizes:

* The file opens with a **pointer table** of 64-bit integers. The first is
  ``3 × n_images``; then, for each image, ``(offset, n_y, n_x)`` -- the
  offset counted in 32-bit words from the start of the file. (Older MBS
  versions write the same table in 32-bit integers; both are read.)
* Each image is ``n_y × n_x`` little-endian ``int32`` counts, row-major:
  ``n_y`` is the analyser's angle channels (``NoS``), ``n_x`` the energy
  steps (``No. Steps``).
* Directly after each image comes its **header**: an ``int32`` byte count,
  then that many bytes of ``key<TAB>value`` lines ending in ``DATA:``.

Axes
----
* **Energy** is kinetic, ``Start K.E. + i × Step Size`` for
  ``i < No. Steps``. In the swept cut this reproduces the ``.txt`` file's
  own energy column to the last digit. In fixed mode it centres the window
  on ``Center K.E.`` to 1e-5 eV, which is how a fixed-mode window is
  defined, so the same formula is used for both.
* **Angle** is ``ScaleMin + j × ScaleMult``. ``(ScaleMax − ScaleMin) /
  ScaleMult + 1`` equals ``NoS`` in every file seen, and that is checked.
* **Map axis**: ``MapStartX … MapEndX`` in ``MapNoXSteps`` steps (or the Y
  equivalents, as ``MapCoordinate`` says).

None of these files records the photon energy, the sample angles, the
temperature or the work function. The energy axis therefore stays
**kinetic**, labelled as such, and anything that needs E_F finds it from the
data (the Fermi fit).

Spin channels
-------------
A spin file carries one spectrum per entry in ``SpinComp#n``, for example::

    SpinComp#0   <0,0> +X  (GUI+Z)
    SpinComp#1   <90,180> +X  (GUI-Z)
    SpinComp#2   <0,0> -X  (GUI-Z)
    SpinComp#3   <90,180> -X  (GUI+Z)

-- the spin-manipulator setting, the target's magnetisation direction, and
the spin direction (in the GUI's frame) that the channel counts as
positive. They are kept verbatim, one per channel, and
:mod:`tools.spin` reads them. Nothing about the spin analysis is decided
here: this file only reads what was measured.

Nothing in here imports Qt.
"""
from __future__ import annotations

import os
import re

import numpy as np

from loader.nxs_file import NxsScan
from loader.registry import Loader, register

__all__ = ["CassiopeeSpinLoader", "read_krx", "read_mbs_text", "load_mbs",
           "parse_header", "classify", "MBS_MARKERS"]

#: Lines an MBS header always carries. A file needs all of them to be ours.
MBS_MARKERS = ("Start K.E.", "Step Size", "No. Steps", "NoS")

#: What the list shows for each kind this reader produces.
KIND_NAMES = {"cut": "Cut", "map": "Map", "spin_edc": "Spin EDC"}

_HEADER_LIMIT = 1 << 20          # no real header is anywhere near a megabyte


# --------------------------------------------------------------------------
# Headers
# --------------------------------------------------------------------------
def parse_header(text: str) -> dict:
    """``key<TAB>value`` lines into a dict.

    The first occurrence of a key wins -- ``RegNo`` appears twice in every
    header, with the same value. ``SpinComp#n`` lines are also gathered, in
    order, under ``"_spin_components"``. Keys are stripped and the value is
    kept as the text it is; :func:`_number` converts where a number is
    wanted.
    """
    header = {}
    components = []
    for line in text.replace("\r\n", "\n").replace("\r", "\n").split("\n"):
        if "\t" not in line:
            continue
        key, value = line.split("\t", 1)
        key, value = key.strip(), value.strip()
        if not key or key == "DATA:":
            continue
        match = re.match(r"SpinComp#(\d+)$", key)
        if match:
            components.append((int(match.group(1)), value))
        header.setdefault(key, value)
    header["_spin_components"] = [value for _, value in sorted(components)]
    return header


def _number(header: dict, key: str, default=None):
    try:
        return float(str(header[key]).strip())
    except (KeyError, TypeError, ValueError):
        return default


def _yes(header: dict, key: str) -> bool:
    return str(header.get(key, "")).strip().upper() in ("YES", "1", "TRUE", "ON")


def _pretty_scale(name: str, fallback: str) -> str:
    """``"Y Angle(Degrees)"`` -> ``"Y angle (deg)"``; the analyser's own
    wording, with the unit where the rest of the program puts it."""
    name = str(name or "").strip()
    if not name:
        return fallback
    name = re.sub(r"\s*\(\s*Degrees\s*\)", " (deg)", name, flags=re.I)
    name = re.sub(r"\s*\(\s*eV\s*\)", " (eV)", name)
    words = name.split()
    # Sentence case, keeping units and one-letter axis names (X, Y) as are.
    return " ".join([words[0]] + [w if w.startswith("(") or len(w) == 1
                                  else w.lower() for w in words[1:]])


# --------------------------------------------------------------------------
# The two file formats, down to (header, image) pairs
# --------------------------------------------------------------------------
def _pointer_table(path: str):
    """``(table, word_bits)``: the ``(offset, n_y, n_x)`` rows and whether
    the table was written in 64- or 32-bit integers. Raises ValueError if
    the start of the file is not a plausible table."""
    size = os.path.getsize(path)
    with open(path, "rb") as handle:
        head = handle.read(8)
    if len(head) < 8:
        raise ValueError("too short to be a .krx file")
    for dtype, width in (("<i8", 8), ("<i4", 4)):
        count = int(np.frombuffer(head[:width], dtype)[0])
        if count <= 0 or count % 3 or count > 3 * 100_000:
            continue
        with open(path, "rb") as handle:
            handle.seek(width)
            table = np.frombuffer(handle.read(width * count), dtype)
        if table.size != count:
            continue
        table = table.reshape(-1, 3).astype(np.int64)
        offsets, ny, nx = table.T
        if np.any(offsets * 4 < width * (count + 1)) or np.any(ny <= 0) \
                or np.any(nx <= 0):
            continue
        if np.any((offsets + ny * nx) * 4 + 4 > size):
            continue
        return table, width * 8
    raise ValueError("no MBS pointer table at the start of the file")


def read_krx(path: str):
    """``[(header_dict, image), ...]`` from an MBS ``.krx`` file.

    ``image`` is ``(n_y, n_x)`` -- angle channels × energy steps -- as a
    read-only memory map of the file, so a 100 MB map costs nothing until
    it is stacked.
    """
    table, _bits = _pointer_table(path)
    words = np.memmap(path, dtype="<i4", mode="r")
    out = []
    with open(path, "rb") as handle:
        for offset, ny, nx in table:
            end = int(offset + ny * nx)
            handle.seek(end * 4)
            length = int(np.frombuffer(handle.read(4), "<i4")[0])
            if not 0 < length < _HEADER_LIMIT:
                raise ValueError(f"image at word {offset}: header length "
                                 f"{length} is not plausible")
            text = handle.read(length).decode("latin-1")
            image = words[int(offset):end].reshape(int(ny), int(nx))
            out.append((parse_header(text), image))
    return out


def _numeric_row(line: str):
    fields = [f for f in line.replace(",", ".").split("\t") if f.strip()]
    if len(fields) < 2:
        return None
    try:
        return [float(f) for f in fields]
    except ValueError:
        return None


def read_mbs_text(path: str):
    """``[(header_dict, image, energy), ...]`` from an MBS ``.txt`` export.

    One block per header: the rows under ``DATA:`` are ``energy, c1, c2 …``,
    one row per energy step. A block's columns become the image's rows, so
    ``image`` is ``(n_columns, n_energy)`` -- the same orientation as the
    ``.krx`` -- and ``energy`` is the file's own first column.
    """
    with open(path, "rb") as handle:
        text = handle.read().decode("latin-1")
    blocks = []
    header_lines, rows = [], []

    def close_block():
        if rows:
            array = np.asarray(rows, dtype=float)
            blocks.append((parse_header("\n".join(header_lines)),
                           np.ascontiguousarray(array[:, 1:].T), array[:, 0]))

    for line in text.replace("\r\n", "\n").replace("\r", "\n").split("\n"):
        row = _numeric_row(line)
        if row is not None and (rows or any(l.startswith("DATA:")
                                            for l in header_lines)):
            if rows and len(row) != len(rows[0]):
                raise ValueError(f"{os.path.basename(path)}: a data row has "
                                 f"{len(row)} fields, the first had "
                                 f"{len(rows[0])}")
            rows.append(row)
            continue
        if rows:                    # a header line after data: a new block
            close_block()
            header_lines, rows = [], []
        if line.strip():
            header_lines.append(line)
    close_block()
    if not blocks:
        raise ValueError(f"{os.path.basename(path)}: no DATA block")
    return blocks


# --------------------------------------------------------------------------
# What the file is
# --------------------------------------------------------------------------
def classify(header: dict, n_images: int) -> str:
    """``"spin_edc"``, ``"map"`` or ``"cut"``, from the header.

    The spin system being on with the main detector off is a spin
    measurement, however many images it has; otherwise more than one image
    with a map step count above one is a map, and one image is a cut.
    """
    if _yes(header, "SpinSystemOn") and not _yes(header, "MainDetOn"):
        return "spin_edc"
    if n_images > 1:
        return "map"
    return "cut"


def energy_axis(header: dict, n: int) -> np.ndarray:
    start = _number(header, "Start K.E.")
    step = _number(header, "Step Size")
    if start is None or step is None:
        raise ValueError("the header has no Start K.E. / Step Size")
    return start + step * np.arange(int(n), dtype=float)


def angle_axis(header: dict, n: int) -> np.ndarray:
    start = _number(header, "ScaleMin")
    step = _number(header, "ScaleMult")
    if start is None or not step:
        return np.arange(int(n), dtype=float)
    return start + step * np.arange(int(n), dtype=float)


def _map_axis(header: dict, n_images: int):
    """The deflector values of a map, its label, and a note if the file
    holds fewer images than the header planned (an interrupted map)."""
    coordinate = str(header.get("MapCoordinate", "")).upper()
    nx = int(_number(header, "MapNoXSteps", 1) or 1)
    ny = int(_number(header, "MapNoYSteps", 1) or 1)
    if nx > 1 and ny > 1:
        raise ValueError(
            f"a {nx} × {ny} two-direction deflector map: this is a 4-D "
            f"measurement, and there is no viewer for it yet")
    use_y = ny > 1 or coordinate.startswith("Y")
    prefix = "Y" if use_y else "X"
    start = _number(header, f"MapStart{prefix}", 0.0)
    end = _number(header, f"MapEnd{prefix}", 0.0)
    planned = ny if use_y else nx
    note = ""
    if planned <= 1:
        planned = n_images
        values = np.arange(n_images, dtype=float)
        label = "Image index"
        return values, label, "the header gives no map steps; images numbered"
    values = np.linspace(start, end, planned)
    if n_images < planned:
        note = (f"interrupted: {n_images} of the {planned} planned "
                f"{prefix} steps are in the file")
        values = values[:n_images]
    elif n_images > planned:
        raise ValueError(f"{n_images} images but the header plans only "
                         f"{planned} map steps")
    label = _pretty_scale(header.get(f"{prefix}MapScaleName"),
                          f"Deflector {prefix} (deg)")
    return values, label, note


def _info(header: dict, path: str) -> dict:
    info = {f"mbs.{key}": value for key, value in header.items()
            if not key.startswith("_")}
    info["title"] = header.get("RegName") or os.path.splitext(
        os.path.basename(path))[0]
    if header.get("STim"):
        info["start_time"] = " ".join(str(header["STim"]).split())
    match = re.match(r"PE0*(\d+(?:\.\d+)?)", str(header.get("Pass Energy", "")))
    if match:
        info["pass_energy_eV"] = float(match.group(1))
    if header.get("Lens Mode"):
        info["lens_mode"] = header["Lens Mode"]
    for key, name in (("DeflX", "mbs.deflector_x_deg"),
                      ("DeflY", "mbs.deflector_y_deg")):
        value = _number(header, key)
        if value is not None:
            info[name] = value
    info["energy_reference"] = "kinetic"
    return info


def _oriented(header: dict, image: np.ndarray) -> np.ndarray:
    """The image as (angle, energy). The energy steps are ``No. Steps``;
    if that is the first dimension the file was written the other way
    round, which no file seen does, but it costs one comparison to be sure."""
    steps = _number(header, "No. Steps")
    if steps is not None and image.shape[1] != int(steps) \
            and image.shape[0] == int(steps):
        return image.T
    return image


# --------------------------------------------------------------------------
# To an NxsScan
# --------------------------------------------------------------------------
def _read(path: str):
    """``[(header, image, energy_or_None)]`` from either format."""
    if path.lower().endswith(".krx"):
        return [(h, img, None) for h, img in read_krx(path)]
    return read_mbs_text(path)


def load_mbs(path: str, progress=None) -> NxsScan:
    """Read one MBS file as a cut, a map or a spin EDC."""
    frames = _read(path)
    header = frames[0][0]
    kind = classify(header, len(frames))
    info = _info(header, path)
    info["mbs.format"] = os.path.splitext(path)[1].lower().lstrip(".")

    if kind == "spin_edc":
        return _spin_scan(path, frames, header, info)

    images = [np.asarray(_oriented(h, img)) for h, img, _ in frames]
    shapes = {img.shape for img in images}
    if len(shapes) != 1:
        raise ValueError(f"the images are not all one size: {sorted(shapes)}")
    n_angle, n_energy = images[0].shape
    energy = frames[0][2] if frames[0][2] is not None else energy_axis(
        header, n_energy)
    if energy.size != n_energy:
        energy = energy_axis(header, n_energy)
    angle = angle_axis(header, n_angle)
    scale_max = _number(header, "ScaleMax")
    if scale_max is not None and n_angle > 1 \
            and abs(angle[-1] - scale_max) > 0.5 * abs(angle[1] - angle[0]):
        info["mbs.note"] = (f"the angle scale ends at {angle[-1]:.4f}, the "
                            f"header says {scale_max:.4f}")
    angle_label = _pretty_scale(header.get("ScaleName"), "Angle (deg)")
    energy_label = _pretty_scale(header.get("EScaleName"), "Kinetic energy (eV)")

    if kind == "cut":
        scan = NxsScan(kind="cut", filename_prefix="")
        scan.x, scan.y = angle, np.asarray(energy, dtype=float)
        scan.value = np.asarray(images[0], dtype=float)
        scan.labels = {"x": angle_label, "y": energy_label}
        scan.info.update(info)
        return scan

    values, map_label, note = _map_axis(header, len(images))
    cube = np.empty((len(images), n_angle, n_energy), dtype=np.float32)
    for index, image in enumerate(images):
        cube[index] = image
        if progress is not None:
            progress(index + 1, len(images), f"image {index + 1}")
    scan = NxsScan(kind="map", filename_prefix="")
    scan.x, scan.k, scan.z = values, angle, np.asarray(energy, dtype=float)
    scan.y = scan.k
    scan.value = cube
    scan.labels = {"x": map_label, "k": angle_label, "z": energy_label}
    scan.info.update(info)
    scan.info["axis0.role"] = "angle" if map_label != "Image index" else "other"
    if note:
        scan.info["mbs.map_note"] = note
    return scan


def _spin_scan(path, frames, header, info) -> NxsScan:
    """One spectrum per spin channel, as a (energy, channel) table."""
    components = list(header.get("_spin_components") or [])
    spectra = []
    energy = None
    for h, image, file_energy in frames:
        image = np.asarray(image, dtype=float)
        if image.ndim != 2:
            raise ValueError("a spin image is not two-dimensional")
        if image.shape[0] == 1 or image.shape[1] == 1:
            spectra.append(image.reshape(-1))
            if energy is None and file_energy is not None:
                energy = np.asarray(file_energy, dtype=float)
        elif len(frames) == 1 and components and \
                image.shape[0] == len(components):
            # A text export with one column per channel in a single block.
            spectra.extend(image)
            if energy is None and file_energy is not None:
                energy = np.asarray(file_energy, dtype=float)
        else:
            raise ValueError(
                f"a spin image is {image.shape[0]} × {image.shape[1]}; a "
                f"spin EDC should be one angle channel wide")
    lengths = {s.size for s in spectra}
    if len(lengths) != 1:
        raise ValueError(f"the spin channels differ in length: {sorted(lengths)}")
    n = lengths.pop()
    if energy is None or energy.size != n:
        energy = energy_axis(header, n)
    table = np.column_stack(spectra)
    n_channels = table.shape[1]
    if len(components) != n_channels:
        components = (components + [f"channel {i}" for i in
                                    range(len(components), n_channels)])[:n_channels]

    scan = NxsScan(kind="spin_edc", filename_prefix="")
    scan.x = np.asarray(energy, dtype=float)
    scan.y = np.arange(n_channels, dtype=float)
    scan.value = table
    scan.labels = {"x": _pretty_scale(header.get("EScaleName"),
                                      "Kinetic energy (eV)"),
                   "y": "Spin channel"}
    scan.info.update(info)
    for index, label in enumerate(components):
        scan.info[f"spin.component.{index}"] = " ".join(str(label).split())
    scan.info["curve.channels"] = "|".join(
        f"C{index} {' '.join(str(label).split())}"
        for index, label in enumerate(components))
    scan.info["curve.value_label"] = "Counts"
    scan.info["spin.n_channels"] = n_channels
    if header.get("ManMagDir"):
        scan.info["spin.magnetisation_now"] = header["ManMagDir"]
    return scan


# --------------------------------------------------------------------------
# The loader
# --------------------------------------------------------------------------
class CassiopeeSpinLoader(Loader):
    name = "SOLEIL CASSIOPEE spin (MBS)"
    description = ("MBS A-1 analyser files from CASSIOPEE's spin-resolved "
                   "end station: .krx or .txt, holding a cut, a deflector "
                   "map, or spin-resolved EDCs (one per spin channel).")
    patterns = ("*.krx", "*.txt")
    priority = 21

    def can_open(self, path: str) -> bool:
        """A ``.krx`` whose pointer table and first header check out, or a
        text file with the MBS header lines. Never the Scienta text files
        the other CASSIOPEE reader takes: those have none of these keys."""
        lower = path.lower()
        if lower.endswith(".krx"):
            try:
                table, _ = _pointer_table(path)
                offset, ny, nx = (int(v) for v in table[0])
                with open(path, "rb") as handle:
                    handle.seek((offset + ny * nx) * 4)
                    length = int(np.frombuffer(handle.read(4), "<i4")[0])
                    if not 0 < length < _HEADER_LIMIT:
                        return False
                    head = handle.read(min(length, 4096)).decode("latin-1")
            except (OSError, ValueError, IndexError):
                return False
            return all(marker in head for marker in MBS_MARKERS[:3])
        if lower.endswith(".txt"):
            try:
                with open(path, "rb") as handle:
                    head = handle.read(4096)
            except OSError:
                return False
            if b"\x00" in head:
                return False
            text = head.decode("latin-1")
            return all(f"\n{marker}\t" in "\n" + text.replace("\r", "")
                       for marker in MBS_MARKERS)
        return False

    def list_entries(self, path: str) -> list:
        """One entry per file; the header says which kind."""
        try:
            if path.lower().endswith(".krx"):
                frames = read_krx(path)
                header, n = frames[0][0], len(frames)
            else:
                frames = read_mbs_text(path)
                header, n = frames[0][0], len(frames)
        except Exception:
            return []
        kind = classify(header, n)
        stem = os.path.splitext(os.path.basename(path))[0]
        region = header.get("RegName") or ""
        return [{
            "entry": kind,
            "kind": KIND_NAMES[kind],
            "name": f"{stem} ({region})" if region and region != stem else stem,
            "title": region,
            "start_time": " ".join(str(header.get("STim", "")).split()) or None,
        }]

    def load(self, path: str, entry: str = None, progress=None):
        return load_mbs(path, progress=progress)


register(CassiopeeSpinLoader())
