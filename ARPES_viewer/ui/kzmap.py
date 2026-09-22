"""
ui/kzmap.py
===========
The kz map processing window: pick one box, calibrate the whole scan.

A photon-energy scan arrives from the loader stacked as measured, on the
first spectrum's energy axis, because nothing in the files says where the
Fermi level of each spectrum actually is (see ``loader/cassiopee.py``). This
is where that gets fixed, using the only reliable reference there is -- the
edge in each spectrum.

The shape of it
---------------
The slit cut -- angle against energy, summed over photon energy -- opens
alongside, and the box dragged on it says two things at once: which angle
channels to sum into an EDC, and over what energy range to fit. That box is
then used *by index* in every spectrum of the scan. Deliberately: an index
range is the same detector channels and the same analyser window everywhere,
while an energy range would mean something different in each spectrum, since
these are precisely the spectra whose energy scales do not yet agree. Using
energies here would make the fit depend on the misalignment it is meant to
measure.

Everything the fitting, shifting and cropping actually does lives in
:mod:`tools.kzmap`, which imports no Qt and can be run from a script. This
file is the window around it.
"""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (QCheckBox, QDialog, QDialogButtonBox,
                             QDoubleSpinBox, QFormLayout, QGroupBox,
                             QHBoxLayout, QLabel, QPushButton, QVBoxLayout,
                             QWidget)

from tools import kzmap
from tools.process import Step, record_step
from ui import jobs as nxs_jobs
from ui.widgets import MemoryData

__all__ = ["KzMapProcessDialog"]


class KzMapProcessDialog(QDialog):
    """Fit the Fermi edge of every spectrum of a kz map and align them.

    Modeless: the box is dragged on another window, so this one must not
    block it.
    """

    datasetsCreated = pyqtSignal(list)

    def __init__(self, contour):
        super().__init__(contour)
        self.setWindowTitle("kz map processing")
        self.setModal(False)
        self.resize(560, 620)
        self.contour = contour
        self.result = None

        self.hv, self.angles, self.energy, _cube = contour.data.angle_cube
        self.hv = np.asarray(self.hv, dtype=float)
        self.angles = np.asarray(self.angles, dtype=float)
        self.energy = np.asarray(self.energy, dtype=float)

        layout = QVBoxLayout(self)
        intro = QLabel(
            "Drag a box on the slit cut that has just opened. It needs to "
            "cover a Fermi edge, with points clearly above and clearly below "
            "it: the angle range says which channels are summed into the EDC, "
            "and the energy range is what the fit runs over.\n\n"
            "That same box is then used, by detector index, on every "
            "spectrum of the scan.")
        intro.setWordWrap(True)
        layout.addWidget(intro)

        # -- what is currently picked ---------------------------------------
        picked = QGroupBox("The box")
        picked_layout = QVBoxLayout(picked)
        self.region_label = QLabel("Nothing picked yet.")
        self.region_label.setWordWrap(True)
        picked_layout.addWidget(self.region_label)
        row = QHBoxLayout()
        self.edge_button = QPushButton("Back to the edge")
        self.edge_button.setToolTip(
            "Every angle channel, and the energy range around where the edge "
            "appears to be in the scan as a whole. This is where the box "
            "starts.")
        self.edge_button.clicked.connect(self.select_edge)
        self.whole_button = QPushButton("Use the whole cut")
        self.whole_button.setToolTip(
            "Every angle channel and the full energy range. Slower, and on a "
            "wide window the far tails pull the fit a little; useful when the "
            "edge is not where the automatic guess put it.")
        self.whole_button.clicked.connect(self.select_everything)
        row.addWidget(self.edge_button)
        row.addWidget(self.whole_button)
        row.addStretch(1)
        picked_layout.addLayout(row)
        layout.addWidget(picked)

        # -- the fit --------------------------------------------------------
        options = QGroupBox("Fitting")
        form = QFormLayout(options)
        self.temperature = QDoubleSpinBox()
        self.temperature.setRange(0.1, 1000.0)
        self.temperature.setDecimals(1)
        self.temperature.setSuffix(" K")
        self.temperature.setValue(_temperature_of(contour.data))
        self.temperature.setToolTip(
            "Held during the fit. On a single edge the temperature and the "
            "resolution are not separable -- the width is "
            "sqrt((3.53kT)^2 + FWHM^2) -- so fitting both gives two numbers "
            "that mean nothing apart.")
        form.addRow("Temperature", self.temperature)

        self.hold_temperature = QCheckBox("Hold the temperature")
        self.hold_temperature.setChecked(True)
        form.addRow("", self.hold_temperature)

        self.normalise = QCheckBox(
            "Normalise each spectrum by its total intensity")
        self.normalise.setChecked(True)
        self.normalise.setToolTip(
            "The photon flux and the analyser transmission both vary by a "
            "large factor across a wide photon-energy scan, so raw counts in "
            "a kz map are largely a picture of the beamline. Dividing each "
            "spectrum by its own total leaves the shape of each, which is "
            "what a kz map is read for.\n\n"
            "Done after aligning and cropping, so every total is a sum over "
            "the same energy range.")
        form.addRow("", self.normalise)
        layout.addWidget(options)

        # -- the answer -----------------------------------------------------
        found = QGroupBox("Fermi level against photon energy")
        found_layout = QVBoxLayout(found)
        self.plot = pg.PlotWidget()
        self.plot.setLabel("bottom", contour.data.scan.labels.get("x", "hv (eV)"))
        self.plot.setLabel("left", "E_F (eV)")
        self.plot.showGrid(x=True, y=True, alpha=0.3)
        self.plot.setMinimumHeight(170)
        found_layout.addWidget(self.plot)
        self.curve = self.plot.plot([], [], pen=pg.mkPen("#1f77b4", width=2),
                                    symbol="o", symbolSize=5,
                                    symbolBrush="#1f77b4")
        self.failed_curve = self.plot.plot(
            [], [], pen=None, symbol="x", symbolSize=9,
            symbolPen=pg.mkPen("#d62728", width=2))
        self.note = QLabel(
            "A smooth curve here means the fits are sound. Scatter means "
            "some of them found something other than the edge.")
        self.note.setWordWrap(True)
        found_layout.addWidget(self.note)
        layout.addWidget(found, stretch=1)

        buttons = QDialogButtonBox()
        self.fit_button = buttons.addButton("Fit and assemble",
                                            QDialogButtonBox.AcceptRole)
        self.fit_button.clicked.connect(self.run)
        buttons.addButton("Close", QDialogButtonBox.RejectRole)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self._open_the_cut()
        self.finished.connect(lambda *_: self._release())

    # -- the slit cut and its box ------------------------------------------
    def _open_the_cut(self):
        """Open (or raise) the slit cut and turn its selection box on."""
        self.contour.open_cut("slit")
        self.cut_window = self.contour.cut_windows.get("slit")
        if self.cut_window is None:
            self.region_label.setText(
                "The slit cut could not be opened, so there is nothing to "
                "pick a box on.")
            self.fit_button.setEnabled(False)
            return
        view = self.cut_window.view
        view.set_selection_visible(True)
        view.selectionChanged.connect(self._box_moved)
        self.select_edge()

    #: How much energy either side of the apparent edge the starting box
    #: covers. Wide enough to hold the edge and a stretch of the flat regions
    #: on both sides, which is what the fit needs, and no wider: on the real
    #: 41-spectrum scan this takes the fitting from 37 s to under 4, and the
    #: answer is slightly *better*, because a window stretching far past the
    #: edge lets the tails pull the background terms around.
    EDGE_HALF_WIDTH_EV = 0.35

    def select_edge(self):
        """Put the box around wherever the edge looks to be.

        Guessed from the whole scan summed into one EDC -- every spectrum has
        the same edge to within the misalignment being measured, so the sum
        locates it well enough to start from even though it is exactly what
        blurs it.
        """
        window = getattr(self, "cut_window", None)
        if window is None:
            return
        lo, hi = float(self.energy.min()), float(self.energy.max())
        try:
            from tools.fermi import initial_guess
            edc = np.nansum(np.asarray(self.contour.full_cube(), dtype=float),
                            axis=(0, 1))
            centre = float(initial_guess(self.energy, edc,
                                         self.temperature.value())["ef"])
        except Exception:                                   # noqa: BLE001
            centre = 0.5 * (lo + hi)        # no guess: the middle will do
        half = self.EDGE_HALF_WIDTH_EV
        low = max(lo, centre - half)
        high = min(hi, centre + half)
        if high - low < 8 * abs(self.energy[1] - self.energy[0]):
            low, high = lo, hi              # too narrow to fit in: take it all
        window.view.set_selection_corners(
            float(self.angles.min()), low, float(self.angles.max()), high)
        self._box_moved()

    def _release(self):
        window = getattr(self, "cut_window", None)
        if window is None:
            return
        try:
            window.view.selectionChanged.disconnect(self._box_moved)
        except TypeError:
            pass          # already gone, which is the normal case on close

    def select_everything(self):
        window = getattr(self, "cut_window", None)
        if window is None:
            return
        window.view.set_selection_corners(
            float(self.angles.min()), float(self.energy.min()),
            float(self.angles.max()), float(self.energy.max()))
        self._box_moved()

    def _box_moved(self, *_):
        region = self.index_region()
        if region is None:
            self.region_label.setText("Nothing picked yet.")
            self.fit_button.setEnabled(False)
            return
        (a0, a1), (e0, e1) = region
        self.region_label.setText(
            f"Angle channels {a0}-{a1}  "
            f"({self.angles[a0]:.3f} to {self.angles[a1]:.3f}), "
            f"{a1 - a0 + 1} summed into the EDC.\n"
            f"Fitting over energy {min(self.energy[e0], self.energy[e1]):.4f} "
            f"to {max(self.energy[e0], self.energy[e1]):.4f} eV "
            f"({abs(e1 - e0) + 1} points).")
        # A fit needs points on both sides of the edge; eight is what
        # tools.fermi itself insists on.
        self.fit_button.setEnabled(abs(e1 - e0) + 1 >= 8)

    def index_region(self):
        """The box as inclusive index bounds, or None."""
        window = getattr(self, "cut_window", None)
        if window is None:
            return None
        return window.view.selection_index_ranges()

    # -- doing it ----------------------------------------------------------
    def run(self):
        region = self.index_region()
        if region is None:
            return None
        cube = np.asarray(self.contour.full_cube(), dtype=float)
        energy = self.energy
        temperature = float(self.temperature.value())
        fixed = ("temperature",) if self.hold_temperature.isChecked() else ()
        normalise = self.normalise.isChecked()
        name = self.contour.filename

        def work(report):
            def tick(done, total):
                report(done / max(1, total),
                       f"Fitting the Fermi edge, spectrum {done + 1} of {total}")
            return kzmap.process_kz_map(
                cube, energy, region, temperature=temperature, fixed=fixed,
                normalise=normalise, progress=tick)

        started = nxs_jobs.run_job(
            self, "kz map processing", work,
            on_done=self._finish,
            on_error=lambda exc: self._failed(exc),
            on_cancel=lambda: self.note.setText("Cancelled."))
        if not started:
            # Another job is already running; run_job says so itself.
            return None
        return True

    def _failed(self, exc):
        self.note.setText(f"<span style='color:#c00'>{exc}</span>")

    def _finish(self, result):
        self.result = result
        self._draw(result)
        dataset = self._dataset(result)
        self.datasetsCreated.emit([dataset])
        self.note.setText(
            result.summary().replace("\n", "<br>")
            + f"<br><b>Added to the list as &ldquo;{dataset.source_label}"
              f"&rdquo;.</b>")

    def _draw(self, result):
        ok = result.ok
        self.curve.setData(self.hv[ok], result.ef[ok])
        self.failed_curve.setData(self.hv[~ok], result.ef[~ok])

    def _dataset(self, result):
        """The aligned cube as a dataset the main list can take."""
        data = self.contour.data
        scan = data.scan
        labels = dict(scan.labels or {})
        labels["z"] = "E - E_F (eV)"

        parameters = {
            "angle_index_from": int(self.index_region()[0][0]),
            "angle_index_to": int(self.index_region()[0][1]),
            "energy_index_from": int(self.index_region()[1][0]),
            "energy_index_to": int(self.index_region()[1][1]),
            "temperature_K": float(self.temperature.value()),
            "temperature_held": bool(self.hold_temperature.isChecked()),
            "normalised_by_total": bool(result.normalised),
            "fermi_level_spread_eV": float(result.spread),
            "energy_trimmed_eV": float(result.trimmed),
            "spectra_fitted": int(result.ok.sum()),
            "spectra_interpolated": int((~result.ok).sum()),
        }
        info = record_step(dict(scan.info),
                           Step("kz_align", parameters, source=self.contour.filename))
        # The per-spectrum Fermi levels are the calibration this produced.
        # They are what makes the result checkable later, and what a second
        # scan of the same sample would be compared against.
        info["kz.fermi_level_eV"] = np.asarray(result.ef, dtype=float)
        info["kz.fermi_level_fitted"] = np.asarray(result.ok, dtype=bool)

        name = _unique(f"{self.contour.filename} [E-Ef]",
                       self.contour.existing_names())
        return MemoryData(
            data.kind, (self.hv, self.angles, result.energy), result.cube,
            labels, source_label=name, parameters=parameters,
            prefix="kz_align", source_path=getattr(data, "path", ""),
            source_info=info, source_motors=dict(scan.fourd_info))


def _temperature_of(data, default: float = 30.0) -> float:
    from ui.windows import metadata_temperature
    return metadata_temperature(data, default)


def _unique(name: str, taken) -> str:
    taken = set(taken or ())
    if name not in taken:
        return name
    index = 2
    while f"{name} ({index})" in taken:
        index += 1
    return f"{name} ({index})"
