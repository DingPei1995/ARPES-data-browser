"""
ui/kzconv.py
============
The kz-to-momentum window: convert a photon-energy scan, and settle V0.

The conversion itself is arithmetic (:mod:`tools.kzconv`). What needs a
window is the inner potential, because the experiment does not measure it:
it is chosen so that the pattern the data shows repeats with the lattice's
own period. So the window is built around that choice --

* a preview of one energy slice, reconverted as the parameters change,
* the Brillouin-zone boundaries drawn over it, from the space group and
  lattice constants,
* a scan of V0 against the period the data shows,
* and a two-point measurement for reading a period straight off the picture
  and asking which lattice planes could produce it.

Each preview is one energy slice, a few milliseconds, which is what lets the
inner potential be a thing you drag rather than a thing you re-run.
"""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (QCheckBox, QComboBox, QDialog, QDialogButtonBox,
                             QDoubleSpinBox, QFormLayout, QGroupBox,
                             QHBoxLayout, QLabel, QPushButton, QSpinBox,
                             QTabWidget, QTextEdit, QVBoxLayout, QWidget)

from tools import cleavage, kzconv
from tools.lattice import LatticeParams, validate_lattice_parameters
from tools.process import Step, record_step
from ui import jobs as nxs_jobs
from tools import colormaps
from ui.widgets import add_button, MemoryData

__all__ = ["KzConversionDialog"]

#: How bent the Fermi edge may be across the analyser angle before the
#: conversion is worth warning about, in eV. A tenth of an electronvolt of
#: bend turns into a visible bend in k_z, and the conversion cannot tell it
#: from dispersion.
FLATNESS_WARNING_EV = 0.03


class KzConversionDialog(QDialog):
    """Convert a kz map to momentum, choosing the inner potential."""

    datasetsCreated = pyqtSignal(list)

    def __init__(self, contour):
        super().__init__(contour)
        self.setWindowTitle("kz map -> momentum")
        self.setModal(False)
        self.resize(940, 820)
        self.contour = contour
        self._pick_points = []
        self._picking = False

        data = contour.data
        self.hv, self.angle, self.energy, _cube = data.angle_cube
        self.hv = np.asarray(self.hv, dtype=float)
        self.angle = np.asarray(self.angle, dtype=float)
        self.energy = np.asarray(self.energy, dtype=float)

        layout = QVBoxLayout(self)
        layout.addWidget(self._flatness_banner())

        body = QHBoxLayout()
        body.addWidget(self._settings_panel(), stretch=0)
        body.addWidget(self._preview_panel(), stretch=1)
        layout.addLayout(body, stretch=1)

        buttons = QDialogButtonBox()
        self.convert_button = add_button(buttons,"Convert the whole cube",
                                                 QDialogButtonBox.AcceptRole)
        self.convert_button.clicked.connect(self.convert)
        add_button(buttons,"Close", QDialogButtonBox.RejectRole)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self.refresh()

    # -- is the Fermi surface flat enough? ---------------------------------
    def _flatness_banner(self) -> QLabel:
        """Say so, up front, if the Fermi edge still bends across the slit.

        Not a refusal: there are reasons to convert an uncorrected cube, and
        a scan taken well away from E_F has no edge to flatten. But the
        conversion reads the energy axis literally, and a bend it inherits
        is a bend in k_z that looks exactly like dispersion -- which is the
        one artefact nobody spots afterwards.
        """
        banner = QLabel()
        banner.setWordWrap(True)
        try:
            cube = np.asarray(self.contour.full_cube(), dtype=float)
            spread, _positions, _angles = kzconv.edge_flatness(
                cube, self.angle, self.energy)
        except Exception:                                   # noqa: BLE001
            banner.setVisible(False)
            return banner
        if spread > FLATNESS_WARNING_EV:
            banner.setStyleSheet(
                "QLabel { background: #fff3cd; color: #664d03; padding: 8px; "
                "border: 1px solid #ffe69c; border-radius: 4px; }")
            banner.setText(
                f"<b>The Fermi edge still bends by {spread:.3f} eV across the "
                f"analyser angle.</b> Straighten it first -- open the slit "
                f"cut and use Fermi-surface correction -- or that bend is "
                f"converted into k_z and comes out looking like dispersion. "
                f"Converting anyway is allowed; this is a reminder, not a "
                f"refusal.")
        else:
            banner.setStyleSheet(
                "QLabel { background: #e7f5e9; color: #14532d; padding: 8px; "
                "border: 1px solid #b7e0c0; border-radius: 4px; }")
            banner.setText(
                f"The Fermi edge is flat across the analyser angle to "
                f"{spread:.3f} eV. Good to convert.")
        return banner

    # -- the settings ------------------------------------------------------
    def _settings_panel(self) -> QWidget:
        holder = QWidget()
        holder.setMaximumWidth(360)
        column = QVBoxLayout(holder)
        column.setContentsMargins(0, 0, 0, 0)

        physics = QGroupBox("Final state")
        form = QFormLayout(physics)
        self.inner_potential = self._spin(12.0, 0.0, 200.0, 0.5, " eV")
        self.inner_potential.setToolTip(
            "V0, the inner potential: how far the final-state band bottom "
            "sits below the vacuum level. Not measured by the experiment -- "
            "chosen so the k_z pattern repeats with the lattice period.")
        form.addRow("Inner potential V0", self.inner_potential)

        self.effective_mass = self._spin(1.0, 0.05, 10.0, 0.05, "")
        self.effective_mass.setToolTip(
            "Final-state effective mass, in units of the free electron mass. "
            "1.0 is the free-electron final state; leave it there unless you "
            "have a reason.")
        form.addRow("Effective mass m*", self.effective_mass)

        self.work_function = self._spin(self._work_function_from_file(),
                                        0.0, 20.0, 0.05, " eV")
        self.work_function.setToolTip(
            "The analyser work function, used to turn a photon energy into a "
            "kinetic energy. Read from the file where the loader recorded "
            "one; change it if you know better.")
        form.addRow("Work function", self.work_function)
        column.addWidget(physics)

        geometry = QGroupBox("Geometry")
        form = QFormLayout(geometry)
        self.angle_offset = self._spin(0.0, -90.0, 90.0, 0.1, " deg")
        self.angle_offset.setToolTip(
            "Where normal emission sits on the analyser slit. Subtracted "
            "from the slit angle before anything else.")
        form.addRow("Normal emission at", self.angle_offset)

        self.theta_position = self._spin(self._theta_from_file(),
                                         -90.0, 90.0, 0.1, " deg")
        self.theta_position.setToolTip(
            "The manipulator's polar angle. It tilts the detector arc out "
            "of the plane holding the slit, which adds a second in-plane "
            "momentum component -- and that one has to come out of k_z too.")
        form.addRow("Manipulator theta", self.theta_position)
        column.addWidget(geometry)

        grid = QGroupBox("Output grid")
        form = QFormLayout(grid)
        self.n_kz = QSpinBox(); self.n_kz.setRange(16, 2048); self.n_kz.setValue(256)
        self.n_kpar = QSpinBox(); self.n_kpar.setRange(16, 2048); self.n_kpar.setValue(256)
        form.addRow("k_z points", self.n_kz)
        form.addRow("k_par points", self.n_kpar)
        self.preview_energy = self._spin(0.0, float(self.energy.min()),
                                          float(self.energy.max()), 0.01, " eV")
        self.preview_energy.setValue(float(np.clip(0.0, self.energy.min(),
                                                   self.energy.max())))
        self.preview_energy.setToolTip(
            "Which energy slice the preview shows. The Fermi level is "
            "usually where the k_z periodicity is clearest.")
        form.addRow("Preview at E", self.preview_energy)
        column.addWidget(grid)

        column.addWidget(self._lattice_group())
        column.addStretch(1)

        for box in (self.inner_potential, self.effective_mass,
                    self.work_function, self.angle_offset,
                    self.theta_position, self.preview_energy):
            box.valueChanged.connect(self.refresh)
        return holder

    def _lattice_group(self) -> QGroupBox:
        group = QGroupBox("Lattice (for the zone lines and the period)")
        form = QFormLayout(group)
        self.space_group = QSpinBox()
        self.space_group.setRange(1, 230)
        self.space_group.setValue(194)
        self.space_group.setToolTip(
            "The space group, which fixes the centering -- and the centering "
            "decides the k_z period: a body-centred lattice repeats every "
            "4*pi/a along [001], not 2*pi/a.")
        form.addRow("Space group", self.space_group)
        self.lat_a = self._spin(3.2, 0.1, 100.0, 0.01, " A")
        self.lat_c = self._spin(6.0, 0.1, 100.0, 0.01, " A")
        form.addRow("a", self.lat_a)
        form.addRow("c", self.lat_c)
        self.normal_index = QComboBox()
        self.normal_index.setToolTip(
            "Which planes the crystal cleaved along. It sets the k_z period "
            "the zone lines are drawn with, and what the V0 scan matches "
            "against. Measure it with the two-point tool if unsure.")
        form.addRow("Surface normal", self.normal_index)
        self.zone_lines = QCheckBox("Draw the zone boundaries")
        self.zone_lines.setChecked(True)
        form.addRow("", self.zone_lines)
        for widget in (self.space_group, self.lat_a, self.lat_c):
            widget.valueChanged.connect(self._lattice_changed)
        self.normal_index.currentIndexChanged.connect(self.refresh)
        self.zone_lines.toggled.connect(self.refresh)
        self._lattice_changed()
        return group

    @staticmethod
    def _spin(value, low, high, step, suffix):
        box = QDoubleSpinBox()
        box.setRange(low, high)
        box.setDecimals(3)
        box.setSingleStep(step)
        box.setValue(value)
        box.setSuffix(suffix)
        box.setKeyboardTracking(False)
        return box

    # -- what the file already knows ---------------------------------------
    def _work_function_from_file(self) -> float:
        info = self.contour.data.scan.info or {}
        for key in ("cassiopee.work_function_eV", "work_function_eV",
                    "analyser_work_function_eV"):
            value = info.get(key)
            if isinstance(value, (int, float)) and np.isfinite(value):
                return float(value)
            if isinstance(value, np.ndarray) and value.size:
                return float(np.nanmean(value))
        return 4.5

    def _theta_from_file(self) -> float:
        info = self.contour.data.scan.info or {}
        value = info.get("cassiopee.sample_theta_deg")
        if isinstance(value, (int, float)) and np.isfinite(value):
            return float(value)
        return 0.0

    # -- the preview -------------------------------------------------------
    def _preview_panel(self) -> QWidget:
        tabs = QTabWidget()

        page = QWidget()
        column = QVBoxLayout(page)
        self.plot = pg.PlotWidget()
        self.plot.setLabel("bottom", "k_z (A^-1)")
        self.plot.setLabel("left", "k_par (A^-1)")
        self.plot.setAspectLocked(True)
        self.image = pg.ImageItem()
        self.plot.addItem(self.image)
        self.zone_item = pg.PlotDataItem(pen=pg.mkPen("#e8e8e8", width=1,
                                                       style=Qt.DashLine))
        self.plot.addItem(self.zone_item)
        self.pick_item = pg.PlotDataItem(
            pen=pg.mkPen("#ff7f0e", width=2), symbol="o", symbolSize=9,
            symbolBrush="#ff7f0e")
        self.plot.addItem(self.pick_item)
        column.addWidget(self.plot, stretch=1)

        row = QHBoxLayout()
        self.pick_button = QPushButton("Measure a k_z period...")
        self.pick_button.setCheckable(True)
        self.pick_button.setToolTip(
            "Click two points on the preview that are the same feature one "
            "zone apart. Their k_z separation is then matched against every "
            "lattice plane, which says which way the crystal cleaved.")
        self.pick_button.toggled.connect(self._toggle_picking)
        row.addWidget(self.pick_button)
        self.scan_button = QPushButton("Scan V0")
        self.scan_button.setToolTip(
            "Convert one energy slice at each of a range of inner "
            "potentials and measure the k_z period each gives. The right V0 "
            "is where it equals the lattice period.")
        self.scan_button.clicked.connect(self.scan_v0)
        row.addWidget(self.scan_button)
        row.addStretch(1)
        column.addLayout(row)

        self.report = QTextEdit()
        self.report.setReadOnly(True)
        self.report.setMaximumHeight(170)
        column.addWidget(self.report)
        tabs.addTab(page, "Preview")

        curve_page = QWidget()
        curve_column = QVBoxLayout(curve_page)
        self.scan_plot = pg.PlotWidget()
        self.scan_plot.setLabel("bottom", "Inner potential V0 (eV)")
        self.scan_plot.setLabel("left", "k_z period (A^-1)")
        self.scan_plot.showGrid(x=True, y=True, alpha=0.3)
        self.scan_curve = self.scan_plot.plot([], [],
                                              pen=pg.mkPen("#1f77b4", width=2))
        self.target_curve = self.scan_plot.plot(
            [], [], pen=pg.mkPen("#d62728", width=2, style=Qt.DashLine))
        curve_column.addWidget(self.scan_plot)
        self.scan_note = QLabel("Press “Scan V0” on the preview tab.")
        self.scan_note.setWordWrap(True)
        curve_column.addWidget(self.scan_note)
        tabs.addTab(curve_page, "V0 scan")
        return tabs

    # -- settings in one place ---------------------------------------------
    def settings(self) -> dict:
        return dict(
            inner_potential=float(self.inner_potential.value()),
            work_function=float(self.work_function.value()),
            effective_mass=float(self.effective_mass.value()),
            angle_offset=float(self.angle_offset.value()),
            theta_position=float(self.theta_position.value()))

    def lattice(self) -> LatticeParams:
        a = float(self.lat_a.value())
        return LatticeParams(a=a, b=a, c=float(self.lat_c.value()),
                             space_group=int(self.space_group.value()))

    def surface_period(self) -> float:
        """The k_z period the chosen surface normal implies, A^-1."""
        data = self.normal_index.currentData()
        if data is None:
            return 2.0 * np.pi / float(self.lat_c.value())
        return float(data[1])

    def _lattice_changed(self, *_):
        """Refill the surface-normal choices for the current lattice."""
        previous = self.normal_index.currentText()
        self.normal_index.blockSignals(True)
        self.normal_index.clear()
        complaints = []
        try:
            params = self.lattice()
            # Returns warnings rather than raising: which of the two the
            # user got wrong -- the lengths or the space group -- is their
            # call, so it is reported beside the choices, not enforced.
            complaints = validate_lattice_parameters(params)
            for hkl, length, _direction in cleavage.reciprocal_lengths(
                    params, max_index=2)[:12]:
                name = "(" + " ".join(str(int(v)) for v in hkl) + ")"
                self.normal_index.addItem(
                    f"{name}   period {length:.4f} A^-1", (hkl, length))
        except Exception as exc:                            # noqa: BLE001
            self.normal_index.addItem(f"lattice not usable: {exc}", None)
        index = self.normal_index.findText(previous)
        if index >= 0:
            self.normal_index.setCurrentIndex(index)
        self.normal_index.blockSignals(False)
        self._lattice_complaints = list(complaints)
        self.refresh()

    # -- drawing -----------------------------------------------------------
    def refresh(self, *_):
        """Reconvert the one preview slice and redraw.

        Does nothing until the preview exists: the settings panel is built
        first and its widgets are already connected to this, so the lattice
        box fires one refresh before there is anything to draw on.
        """
        if not hasattr(self, "image"):
            return
        try:
            kz_axis, kpar_axis, plane = self._preview_slice()
        except Exception as exc:                            # noqa: BLE001
            self.report.setHtml(f"<span style='color:#c00'>{exc}</span>")
            return
        self._preview = (kz_axis, kpar_axis, plane)
        shown = np.where(np.isfinite(plane), plane, np.nan)
        self.image.setImage(shown, autoLevels=True)
        self.image.setRect(pg.QtCore.QRectF(
            kz_axis[0], kpar_axis[0],
            kz_axis[-1] - kz_axis[0], kpar_axis[-1] - kpar_axis[0]))
        self._apply_colormap()
        self._draw_zone_lines(kz_axis, kpar_axis)
        self._write_report(kz_axis, kpar_axis, plane)

    def _apply_colormap(self):
        """The contour's own colormap on the preview.

        ``ui.widgets.apply_colormap`` drives a ``pg.ImageView``; this is a
        bare ``ImageItem``, whose colour table is set directly. It is the
        same lookup table either way, so the preview and the window that
        opened it match.
        """
        try:
            lut = colormaps.get_lut(self.contour.colormap,
                                    flip=self.contour.flip)
            self.image.setLookupTable(np.asarray(lut))
        except Exception:                                   # noqa: BLE001
            pass        # a preview in the default colours is no disaster

    def _preview_slice(self):
        cube = np.asarray(self.contour.full_cube(), dtype=float)
        index = int(np.argmin(np.abs(self.energy
                                     - float(self.preview_energy.value()))))
        one = np.asarray([float(self.energy[index])])
        kz_axis, kpar_axis, _e, out = kzconv.to_kz_cube(
            self.hv, self.angle, one, cube[:, :, index][:, :, None],
            n_kz=int(self.n_kz.value()), n_kpar=int(self.n_kpar.value()),
            **self.settings())
        return kz_axis, kpar_axis, out[:, :, 0]

    def _draw_zone_lines(self, kz_axis, kpar_axis):
        """Zone boundaries as vertical lines at multiples of the period.

        A k_z scan cuts the Brillouin zone along one direction, so what the
        zone contributes here is a set of planes perpendicular to k_z --
        drawn as one polyline with gaps so a single item carries all of them.
        """
        if not self.zone_lines.isChecked():
            self.zone_item.setData([], [])
            return
        period = self.surface_period()
        if not np.isfinite(period) or period <= 0:
            self.zone_item.setData([], [])
            return
        first = int(np.floor(kz_axis[0] / period))
        last = int(np.ceil(kz_axis[-1] / period))
        xs, ys = [], []
        for n in range(first, last + 1):
            position = n * period
            if not (kz_axis[0] <= position <= kz_axis[-1]):
                continue
            xs += [position, position, np.nan]
            ys += [kpar_axis[0], kpar_axis[-1], np.nan]
        self.zone_item.setData(np.array(xs), np.array(ys), connect="finite")

    def _write_report(self, kz_axis, kpar_axis, plane):
        period = self.surface_period()
        coverage = 100.0 * float(np.isfinite(plane).mean())
        lines = [
            f"k_z {kz_axis[0]:.3f} to {kz_axis[-1]:.3f} A<sup>-1</sup> "
            f"({(kz_axis[-1] - kz_axis[0]) / period:.2f} zones of "
            f"{period:.4f} A<sup>-1</sup>)",
            f"k_par {kpar_axis[0]:.3f} to {kpar_axis[-1]:.3f} A<sup>-1</sup>",
            f"{coverage:.0f}% of the grid is covered by the measurement",
        ]
        for complaint in getattr(self, "_lattice_complaints", ()):
            lines.append(f"<span style='color:#a60'>{complaint}</span>")
        if self._pick_points:
            lines.append(self._pick_summary())
        self.report.setHtml("<br>".join(lines))

    # -- picking two points ------------------------------------------------
    def _toggle_picking(self, on: bool):
        self._picking = bool(on)
        self.pick_button.setText("Click two points..." if on
                                 else "Measure a k_z period...")
        if on:
            self._pick_points = []
            self.pick_item.setData([], [])
            self.plot.scene().sigMouseClicked.connect(self._clicked)
        else:
            try:
                self.plot.scene().sigMouseClicked.disconnect(self._clicked)
            except TypeError:
                pass
        self.refresh()

    def _clicked(self, event):
        if not self._picking:
            return
        point = self.plot.getPlotItem().vb.mapSceneToView(event.scenePos())
        self._pick_points.append((float(point.x()), float(point.y())))
        if len(self._pick_points) > 2:
            self._pick_points = self._pick_points[-2:]
        xs = [p[0] for p in self._pick_points]
        ys = [p[1] for p in self._pick_points]
        self.pick_item.setData(xs, ys)
        if len(self._pick_points) == 2:
            self.pick_button.setChecked(False)      # also stops the picking
        else:
            self.refresh()

    def _pick_summary(self) -> str:
        if len(self._pick_points) < 2:
            return "One point picked; click the same feature one zone along."
        (x0, _y0), (x1, _y1) = self._pick_points
        distance = abs(x1 - x0)
        if distance <= 0:
            return ("<span style='color:#c00'>Those two points are at the "
                    "same k_z.</span>")
        lines = [f"<b>Picked k_z separation {distance:.4f} A<sup>-1</sup></b> "
                 f"(real-space repeat {2 * np.pi / distance:.3f} A)"]
        try:
            found = cleavage.candidates(distance, self.lattice(),
                                        tolerance=0.15)
        except Exception as exc:                            # noqa: BLE001
            return "<br>".join(lines + [f"<span style='color:#c00'>{exc}</span>"])
        if not found:
            lines.append("No lattice plane matches that within 15%. Either "
                         "the lattice constants are wrong, or the two points "
                         "are not one zone apart.")
        else:
            lines.append("Planes that match within 15%:")
            for candidate in found[:6]:
                lines.append("&nbsp;&nbsp;" + candidate.describe())
        return "<br>".join(lines)

    # -- scanning V0 -------------------------------------------------------
    def scan_v0(self):
        cube = np.asarray(self.contour.full_cube(), dtype=float)
        period = self.surface_period()
        spacing = 2.0 * np.pi / period
        settings = self.settings()
        settings.pop("inner_potential")
        binding = float(self.preview_energy.value())

        def work(report):
            def tick(done, total):
                report(done / max(1, total),
                       f"Trying inner potentials ({done + 1} of {total})")
            return kzconv.scan_inner_potential(
                self.hv, self.angle, self.energy, cube, spacing=spacing,
                binding_energy=binding, progress=tick, **settings)

        started = nxs_jobs.run_job(
            self, "Scanning the inner potential", work,
            on_done=self._show_scan,
            on_error=lambda exc: self.scan_note.setText(
                f"<span style='color:#c00'>{exc}</span>"))
        if not started:
            self.scan_note.setText("Another operation is running.")

    def _show_scan(self, result):
        good = np.isfinite(result.periods)
        self.scan_curve.setData(result.inner_potentials[good],
                                result.periods[good])
        self.target_curve.setData(
            [result.inner_potentials[0], result.inner_potentials[-1]],
            [result.target, result.target])
        if result.best is None:
            self.scan_note.setText(
                "The period never matches the lattice over this range of "
                "inner potentials. Check the surface normal and the lattice "
                "constants, or widen the range.")
            return
        uncertainty = result.uncertainty()
        self.scan_note.setText(
            f"<b>The k_z period matches {result.target:.4f} A<sup>-1</sup> at "
            f"V0 = {result.best:.1f} eV.</b><br>"
            f"Sensitivity {result.sensitivity:.5f} A<sup>-1</sup> per eV, so "
            f"measuring the period to 2% pins V0 only to about "
            f"±{uncertainty:.1f} eV. That is this scan's photon-energy range "
            f"talking, not the fit: covering more zones would narrow it. "
            f"Use it to start from, then look at whether the pattern lines "
            f"up with the zone boundaries.")
        self.inner_potential.setValue(float(result.best))

    # -- converting --------------------------------------------------------
    def convert(self):
        cube = np.asarray(self.contour.full_cube(), dtype=float)
        settings = self.settings()
        n_kz, n_kpar = int(self.n_kz.value()), int(self.n_kpar.value())
        hv, angle, energy = self.hv, self.angle, self.energy

        def work(report):
            def tick(done, total):
                report(done / max(1, total),
                       f"Converting energy slice {done + 1} of {total}")
            return kzconv.to_kz_cube(hv, angle, energy, cube, n_kz=n_kz,
                                     n_kpar=n_kpar, progress=tick, **settings)

        started = nxs_jobs.run_job(
            self, "Converting to momentum", work, on_done=self._finish,
            on_error=lambda exc: self.report.setHtml(
                f"<span style='color:#c00'>{exc}</span>"))
        if not started:
            self.report.setHtml("Another operation is running.")

    def _finish(self, converted):
        kz_axis, kpar_axis, energy, out = converted
        data = self.contour.data
        scan = data.scan
        settings = self.settings()
        parameters = dict(settings)
        parameters.update({
            "n_kz": int(self.n_kz.value()),
            "n_kpar": int(self.n_kpar.value()),
            "space_group": int(self.space_group.value()),
            "lattice_a": float(self.lat_a.value()),
            "lattice_c": float(self.lat_c.value()),
            "surface_period_invA": float(self.surface_period()),
        })
        info = record_step(dict(scan.info),
                           Step("kz_to_k", parameters,
                                source=self.contour.filename))
        info["kz.inner_potential_eV"] = float(settings["inner_potential"])

        name = _unique(f"{self.contour.filename} [kz]",
                       self.contour.existing_names())
        result = MemoryData(
            "kz_map_k", (kz_axis, kpar_axis, energy), out,
            {"x": "k_z (Å⁻¹)", "k": "k_par (Å⁻¹)",
             "z": scan.labels.get("z", "E - E_F (eV)")},
            source_label=name, parameters=parameters, prefix="kz_to_k",
            source_path=getattr(data, "path", ""), source_info=info,
            source_motors=dict(scan.fourd_info))
        self.datasetsCreated.emit([result])
        coverage = 100.0 * float(np.isfinite(out).mean())
        self.report.setHtml(
            f"<b>Converted to (k_z, k_par, E) at V0 = "
            f"{settings['inner_potential']:.1f} eV.</b><br>"
            f"k_z {kz_axis[0]:.3f} to {kz_axis[-1]:.3f} Å⁻¹, "
            f"k_par {kpar_axis[0]:.3f} to {kpar_axis[-1]:.3f} Å⁻¹, "
            f"{coverage:.0f}% covered.<br>"
            f"Added to the list as “{name}”.")


def _unique(name: str, taken) -> str:
    taken = set(taken or ())
    if name not in taken:
        return name
    index = 2
    while f"{name} ({index})" in taken:
        index += 1
    return f"{name} ({index})"
