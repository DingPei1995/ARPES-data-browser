"""
ui/panels/fsview.py
===================
A 3-D Fermi surface on screen: the surfaces of
:func:`tools.crystal.fermisurface.fermi_surface`, the Brillouin zone's
edges and its named points, turned with the mouse.

Drawn with QPainter, triangles sorted back to front and shaded by their
normal -- the same no-OpenGL approach as the program's other 3-D views
(:mod:`ui.panels.volume`), so it looks the same on every machine and the
picture for the figure composer is exactly what is on screen. While the
mouse turns it, a coarser mesh is drawn, so turning stays smooth; the fine
one comes back on release.
"""
from __future__ import annotations

import numpy as np
from PyQt5.QtCore import QPointF, Qt, pyqtSignal
from PyQt5.QtGui import (QColor, QFont, QIcon, QImage, QPainter, QPen, QPixmap,
                         QPolygonF)
from PyQt5.QtWidgets import (QAbstractItemView, QCheckBox, QColorDialog,
                             QComboBox, QDialog, QDoubleSpinBox, QFileDialog,
                             QHBoxLayout, QLabel, QPushButton, QSizePolicy,
                             QTableWidget, QTableWidgetItem, QVBoxLayout,
                             QWidget)

__all__ = ["FermiSurfaceView", "FermiSurfaceWindow", "BandStyleTable",
           "ZoneWindow", "BAND_COLOURS", "VIEWS"]

#: the bands' colours, in order (the theory side's MATLAB scripts start
#: with blue, sky blue, green and gold)
BAND_COLOURS = ("#0000ff", "#00bfff", "#32cd32", "#ffd700", "#d62728",
                "#9467bd", "#8c564b", "#e377c2", "#17becf", "#7f7f7f")

#: named views: (azimuth, elevation), degrees
VIEWS = {"3D": (60.0, 20.0), "top (along k_z)": (0.0, 90.0),
         "side (along k_x)": (0.0, 0.0), "side (along k_y)": (90.0, 0.0)}

#: where the light comes from, in screen coordinates (x right, y up,
#: z towards the viewer)
_LIGHT = np.array([-0.35, 0.45, 0.82])
_LIGHT /= np.linalg.norm(_LIGHT)


def _rotation(azimuth: float, elevation: float) -> np.ndarray:
    """Rows: the screen's x (right), y (up) and z (towards the viewer) in
    the crystal's Cartesian frame, looking from ``azimuth`` (about k_z,
    from k_x) and ``elevation`` (above the k_x-k_y plane)."""
    a, e = np.deg2rad(azimuth), np.deg2rad(elevation)
    towards = np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])
    up_guess = np.array([0.0, 0.0, 1.0]) if abs(np.sin(e)) < 0.999 else \
        np.array([-np.cos(a), -np.sin(a), 0.0])
    right = np.cross(up_guess, towards)
    right /= np.linalg.norm(right)
    up = np.cross(towards, right)
    return np.array([right, up, towards])


class FermiSurfaceView(QWidget):
    """Draws a :class:`~tools.crystal.fermisurface.FermiSurface`."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(360, 320)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.surface = None
        self.coarse = None
        self.points = []            # [(label, Cartesian k[, colour])]
        #: band index -> {"colour", "opacity", "visible"}; a band not listed
        #: takes the next colour of BAND_COLOURS, opaque
        self.styles = {}
        #: polygons drawn in the zone, translucent: [(points (N, 3), colour)]
        self.planes = []
        self.azimuth, self.elevation = VIEWS["3D"]
        self.zoom = 1.0
        self.show_zone = True
        self.show_labels = True
        self.message = "Nothing drawn yet."
        self._drag = None
        self.setMouseTracking(False)

    # -- content ------------------------------------------------------------
    def set_surface(self, surface, coarse=None, points=()):
        self.surface = surface
        self.coarse = coarse
        self.points = list(points)
        self.update()

    def style_of(self, order: int, band: int) -> dict:
        """How band ``band`` (the ``order``-th surface) is drawn."""
        style = {"colour": BAND_COLOURS[order % len(BAND_COLOURS)],
                 "opacity": 1.0, "visible": True}
        style.update(self.styles.get(band, {}))
        return style

    def set_view(self, name_or_angles):
        if isinstance(name_or_angles, str):
            name_or_angles = VIEWS[name_or_angles]
        self.azimuth, self.elevation = name_or_angles
        self.update()

    # -- mouse --------------------------------------------------------------
    def mousePressEvent(self, event):
        self._drag = (event.pos(), self.azimuth, self.elevation)

    def mouseMoveEvent(self, event):
        if self._drag is None:
            return
        start, azimuth, elevation = self._drag
        delta = event.pos() - start
        self.azimuth = azimuth - 0.5 * delta.x()
        self.elevation = float(np.clip(elevation + 0.5 * delta.y(), -90, 90))
        self.update()

    def mouseReleaseEvent(self, event):
        self._drag = None
        self.update()

    def mouseDoubleClickEvent(self, event):
        self.zoom = 1.0
        self.set_view("3D")

    def wheelEvent(self, event):
        self.zoom = float(np.clip(self.zoom * (1.1 ** (event.angleDelta().y() / 120)),
                                  0.3, 6.0))
        self.update()

    # -- drawing ------------------------------------------------------------
    def paintEvent(self, _event):
        painter = QPainter(self)
        self.paint(painter, self.width(), self.height(),
                   coarse=self._drag is not None)
        painter.end()

    def render_image(self, width: int = 1400, height: int = 1200) -> QImage:
        """The picture as on screen, at ``width`` x ``height`` pixels, on
        white -- for the figure composer."""
        image = QImage(width, height, QImage.Format_ARGB32)
        image.fill(QColor("white"))
        painter = QPainter(image)
        self.paint(painter, width, height, coarse=False, scale_text=width / max(self.width(), 1))
        painter.end()
        return image

    def paint(self, painter, width, height, coarse=False, scale_text=1.0):
        painter.setRenderHint(QPainter.Antialiasing, True)
        painter.fillRect(0, 0, width, height, QColor("white"))
        surface = self.coarse if (coarse and self.coarse is not None) else self.surface
        if surface is None or surface.zone is None:
            painter.setPen(QColor("#666666"))
            painter.drawText(0, 0, width, height, Qt.AlignCenter | Qt.TextWordWrap,
                             self.message)
            return
        rot = _rotation(self.azimuth, self.elevation)
        radius = float(np.linalg.norm(surface.zone.vertices, axis=1).max())
        scale = 0.45 * min(width, height) / radius * self.zoom
        cx, cy = width / 2.0, height / 2.0

        def screen(points):
            p = np.asarray(points, float) @ rot.T
            return cx + scale * p[..., 0], cy - scale * p[..., 1], p[..., 2]

        edges = self._zone_edges(surface.zone) if self.show_zone else []
        back = [e for e in edges if np.mean(screen(e)[2]) < 0]
        front = [e for e in edges if np.mean(screen(e)[2]) >= 0]
        self._draw_edges(painter, back, screen, QColor(0, 0, 0, 90), Qt.DashLine)

        triangles, colours, normals = [], [], []
        styles = [self.style_of(order, item.band)
                  for order, item in enumerate(surface.surfaces)]
        for order, item in enumerate(surface.surfaces):
            if not styles[order]["visible"]:
                continue
            triangles.append(item.triangles)
            normals.append(item.normals)
            colours.append(np.full(len(item.triangles), order))
        if triangles:
            triangles = np.concatenate(triangles)
            normals = np.concatenate(normals) @ rot.T
            which = np.concatenate(colours)
            xs, ys, zs = screen(triangles)
            depth = zs.mean(1)
            lit = np.abs(normals @ _LIGHT)
            shade = 0.3 + 0.65 * lit
            gloss = 110.0 * lit ** 14
            base = np.array([QColor(st["colour"]).getRgb()[:3] for st in styles], float)
            alphas = np.array([int(255 * np.clip(st["opacity"], 0.05, 1.0))
                               for st in styles])[which]
            rgb = np.clip(base[which] * shade[:, None] + gloss[:, None], 0, 255).astype(int)
            painter.setPen(Qt.NoPen)
            for t in np.argsort(depth):
                r, g, b = rgb[t]
                alpha = int(alphas[t])
                colour = QColor(int(r), int(g), int(b), alpha)
                painter.setBrush(colour)
                # an outline in the fill colour closes the hairline gaps
                # antialiasing leaves between opaque neighbours
                painter.setPen(QPen(colour, 0.6) if alpha == 255 else Qt.NoPen)
                painter.drawPolygon(QPolygonF([QPointF(xs[t, 0], ys[t, 0]),
                                               QPointF(xs[t, 1], ys[t, 1]),
                                               QPointF(xs[t, 2], ys[t, 2])]))
        for polygon, plane_colour in self.planes:
            xs, ys, _z = screen(polygon)
            fill = QColor(plane_colour)
            edge = QColor(plane_colour)
            fill.setAlpha(70)
            painter.setBrush(fill)
            painter.setPen(QPen(edge, 1.6))
            painter.drawPolygon(QPolygonF([QPointF(float(x), float(y))
                                           for x, y in zip(xs, ys)]))
        self._draw_edges(painter, front, screen, QColor(0, 0, 0, 200), Qt.SolidLine)
        if self.show_labels and self.points:
            font = QFont()
            font.setPointSizeF(11 * scale_text)
            painter.setFont(font)
            for entry in self.points:
                label, k = entry[0], entry[1]
                dot = QColor(entry[2]) if len(entry) > 2 and entry[2] else QColor("black")
                x, y, _z = screen(np.asarray(k)[None])
                painter.setPen(Qt.NoPen)
                painter.setBrush(dot)
                radius = (4.5 if len(entry) > 2 and entry[2] else 3) * scale_text
                painter.drawEllipse(QPointF(float(x[0]), float(y[0])), radius, radius)
                painter.setPen(dot)
                painter.drawText(QPointF(float(x[0]) + 5 * scale_text,
                                         float(y[0]) - 5 * scale_text), label)

    def _zone_edges(self, zone):
        seen, edges = set(), []
        for face in zone.faces:
            for p, q in zip(face[:-1], face[1:]):
                key = tuple(sorted((tuple(np.round(p, 6)), tuple(np.round(q, 6)))))
                if key not in seen:
                    seen.add(key)
                    edges.append(np.array([p, q]))
        return edges

    @staticmethod
    def _draw_edges(painter, edges, screen, colour, style):
        pen = QPen(colour, 1.2)
        pen.setStyle(style)
        painter.setPen(pen)
        for edge in edges:
            xs, ys, _z = screen(edge)
            painter.drawLine(QPointF(xs[0], ys[0]), QPointF(xs[1], ys[1]))


class FermiSurfaceWindow(QDialog):
    """The same Fermi surface in a window of its own, as large as the
    screen allows: the surfaces, the bands hidden, the opacity and the
    view are those of ``source`` (a :class:`FermiSurfaceView`) when it
    opens; the parameters are set in the Calculated bands window."""

    def __init__(self, source: FermiSurfaceView, title: str = "3D Fermi surface",
                 parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(1000, 860)
        layout = QVBoxLayout(self)
        self.view = FermiSurfaceView()
        for name in ("surface", "coarse", "points", "planes", "azimuth",
                     "elevation", "zoom", "show_zone", "show_labels", "message"):
            setattr(self.view, name, getattr(source, name))
        self.view.styles = {k: dict(v) for k, v in source.styles.items()}
        layout.addWidget(self.view, 1)
        row = QHBoxLayout()
        row.addWidget(QLabel("<i>drag to turn, wheel to zoom</i>   View:"))
        views = QComboBox()
        views.addItems(list(VIEWS))
        views.activated.connect(lambda *_: self.view.set_view(views.currentText()))
        row.addWidget(views)
        row.addStretch(1)
        save = QPushButton("Save picture...")
        save.setAutoDefault(False)
        save.clicked.connect(lambda: self.save_picture())
        row.addWidget(save)
        close = QPushButton("Close")
        close.setAutoDefault(False)
        close.clicked.connect(self.close)
        row.addWidget(close)
        layout.addLayout(row)

    def save_picture(self, path: str = None):
        if path is None:
            path, _ = QFileDialog.getSaveFileName(self, "Save the picture",
                                                  "fermi_surface.png",
                                                  "PNG (*.png)")
            if not path:
                return None
        if not path.lower().endswith(".png"):
            path += ".png"
        self.view.render_image(2000, 1700).save(path)
        return path


#: the colours offered for a band, besides "Other colour..."
OFFERED = BAND_COLOURS + ("#000000", "#808080")


class BandStyleTable(QTableWidget):
    """One row per band: its colour (any, and the same for several bands --
    to group electron and hole pockets, say), its opacity (for a surface)
    and whether it is shown. Writes into ``styles`` (band index ->
    ``{"colour", "opacity", "visible"}``), which several tables can share,
    so a band keeps its colour from the contours to the 3-D view."""

    changed = pyqtSignal()

    def __init__(self, styles: dict, opacity: bool = True, parent=None):
        super().__init__(0, 4 if opacity else 3, parent)
        self.styles = styles
        self.with_opacity = opacity
        self.setHorizontalHeaderLabels(["Band", "Colour"] +
                                       (["Opacity"] if opacity else []) + ["Show"])
        self.verticalHeader().setVisible(False)
        self.horizontalHeader().setStretchLastSection(True)
        self.setSelectionMode(QAbstractItemView.NoSelection)
        self.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.setMaximumHeight(170)
        self.bands = []

    def style(self, band: int, order: int) -> dict:
        """``band``'s style, given the next free colour when it has none."""
        if band not in self.styles:
            self.styles[band] = {"colour": BAND_COLOURS[order % len(BAND_COLOURS)],
                                 "opacity": 1.0, "visible": True}
        return self.styles[band]

    def set_bands(self, entries):
        """``entries``: ``[(band index, label)]``."""
        self.bands = [band for band, _label in entries]
        self.setRowCount(len(entries))
        for row, (band, label) in enumerate(entries):
            style = self.style(band, row)
            item = QTableWidgetItem(label)
            self.setItem(row, 0, item)
            combo = QComboBox()
            for colour in OFFERED:
                combo.addItem(self._swatch(colour), colour, colour)
            combo.addItem("Other colour...", None)
            self._select(combo, style["colour"])
            combo.activated.connect(lambda index, b=band, c=combo: self._colour_chosen(b, c, index))
            self.setCellWidget(row, 1, combo)
            column = 2
            if self.with_opacity:
                spin = QDoubleSpinBox()
                spin.setRange(0.1, 1.0)
                spin.setSingleStep(0.1)
                spin.setDecimals(2)
                spin.setKeyboardTracking(False)
                spin.setValue(float(style["opacity"]))
                spin.setToolTip("1 opaque; lower lets what is behind show through.")
                spin.valueChanged.connect(lambda v, b=band: self._set(b, "opacity", float(v)))
                self.setCellWidget(row, column, spin)
                column += 1
            box = QCheckBox()
            box.setChecked(bool(style["visible"]))
            box.toggled.connect(lambda on, b=band: self._set(b, "visible", bool(on)))
            self.setCellWidget(row, column, box)
        self.resizeColumnsToContents()

    @staticmethod
    def _swatch(colour: str) -> QIcon:
        pixmap = QPixmap(14, 14)
        pixmap.fill(QColor(colour))
        return QIcon(pixmap)

    def _select(self, combo, colour: str):
        index = combo.findData(colour)
        if index < 0:
            combo.insertItem(combo.count() - 1, self._swatch(colour), colour, colour)
            index = combo.findData(colour)
        combo.setCurrentIndex(index)

    def _colour_chosen(self, band, combo, index):
        colour = combo.itemData(index)
        if colour is None:
            picked = QColorDialog.getColor(QColor(self.styles[band]["colour"]), self,
                                           "The band's colour")
            if not picked.isValid():
                self._select(combo, self.styles[band]["colour"])
                return
            colour = picked.name()
            self._select(combo, colour)
        self._set(band, "colour", colour)

    def set_colour(self, band: int, colour: str):
        """Give ``band`` ``colour`` (as choosing it in its row does)."""
        if band in self.bands:
            self._select(self.cellWidget(self.bands.index(band), 1), colour)
        self._set(band, "colour", colour)

    def _set(self, band, key, value):
        self.styles.setdefault(band, {"colour": BAND_COLOURS[0], "opacity": 1.0,
                                      "visible": True})[key] = value
        self.changed.emit()


class ZoneWindow(QDialog):
    """Where in k-space a cut lies: the Brillouin zone with its named
    points, the cut plane in it (and the two faces of its k⊥ window, along
    the plane's normal), and the points the plane holds, ringed in red.
    Turned with the mouse; :meth:`set_plane` follows the cut as it changes."""

    PLANE = "#d62728"
    SLAB = "#ff9896"

    def __init__(self, zone, points, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Where the cut lies in the Brillouin zone")
        self.resize(640, 600)
        self.zone = zone
        self.all_points = list(points)
        layout = QVBoxLayout(self)
        self.view = FermiSurfaceView()
        from arpes_viewer.tools.crystal.fermisurface import FermiSurface

        self.view.set_surface(FermiSurface(energy=0.0, zone=zone))
        layout.addWidget(self.view, 1)
        self.label = QLabel("")
        self.label.setWordWrap(True)
        layout.addWidget(self.label)
        row = QHBoxLayout()
        row.addWidget(QLabel("<i>drag to turn, wheel to zoom</i>   View:"))
        views = QComboBox()
        views.addItems(list(VIEWS))
        views.activated.connect(lambda *_: self.view.set_view(views.currentText()))
        row.addWidget(views)
        row.addStretch(1)
        close = QPushButton("Close")
        close.setAutoDefault(False)
        close.clicked.connect(self.close)
        row.addWidget(close)
        layout.addLayout(row)

    def _polygon(self, normal, offset):
        from arpes_viewer.tools.crystal.bz3d import cut_points_3d, plane_cut

        u, v = plane_cut(self.zone.faces, normal, offset * normal)
        if u is None:
            return None
        return cut_points_3d(u, v, normal, offset * normal)

    def set_plane(self, normal, offset: float, window: float = 0.0, hkl=None):
        normal = np.asarray(normal, float)
        planes = []
        if window > 0:
            for shift in (-window / 2, window / 2):
                polygon = self._polygon(normal, offset + shift)
                if polygon is not None:
                    planes.append((polygon, self.SLAB))
        polygon = self._polygon(normal, offset)
        if polygon is not None:
            planes.append((polygon, self.PLANE))
        self.view.planes = planes
        on = [p.label for p in self.all_points
              if abs(float(np.dot(p.k, normal)) - offset) < 1e-4]
        self.view.points = [(p.label, p.k, self.PLANE if
                             abs(float(np.dot(p.k, normal)) - offset) < 1e-4 else None)
                            for p in self.all_points]
        name = f"({' '.join(str(int(v)) for v in hkl)}) " if hkl is not None else ""
        self.label.setText(
            f"The plane {name}at k⊥ = {offset:.4f} Å⁻¹ along its normal"
            + (f", integrated over ±{window / 2:.3f} Å⁻¹ (the light faces)"
               if window > 0 else "")
            + (". It holds " + ", ".join(dict.fromkeys(on)) if on else
               ". It holds no named point")
            + ("." if polygon is not None else " -- and misses the first zone."))
        self.view.update()
