"""
ui/launcher/main_window.py
==========================
The launcher window's layout: a file picker, the global colormap controls,
and the metadata table for whichever file is selected. Hand-written PyQt5,
no Qt-Designer ``.ui``/``.qrc`` needed.

Two boxes side by side in one window: the tool panel on the left
(``ToolPanel``: everything that analyses or processes the selected datasets,
grouped by what dimension of data it takes, then the memory controls), and
the main panel on the right (the list, the operations that change the data
directly, and its summary). The tool panel is part of the window, so it
moves, resizes and closes with it and has no title bar of its own.

Visualisation is NOT here -- double-clicking a file in the list opens a
viewer window from ``ui/viewers/windows.py``. Keeping the launcher small
means several files (and several cuts of one file) can be open side by side
and arranged freely, which the single stacked layout this replaced could not
do.

Wiring lives in ``ARPES_viewer.py`` (``SetConnect()``), the same split the
original app used.
"""
from PyQt5 import QtCore, QtWidgets

import version


class Ui_MainWindow:
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 900)
        MainWindow.setWindowTitle(version.title())
        # The panel's colours. The disabled rules are not optional: a colour
        # set for every widget is also the colour of a *disabled* widget and
        # of a disabled menu entry, so without them a greyed-out entry of the
        # list's right-click menu looked exactly like an enabled one.
        MainWindow.setStyleSheet(
            "* { background-color: rgb(245, 241, 249);"
            " font: \"Lucida Sans Unicode\";"
            " color: rgb(102, 126, 161); }"
            "\n*:disabled { color: rgb(190, 190, 200); }"
            "\nQMenu::item { padding: 3px 24px 3px 20px; }"
            "\nQMenu::item:selected { background-color: rgb(214, 224, 240);"
            " color: rgb(40, 60, 95); }"
            "\nQMenu::item:disabled { color: rgb(190, 190, 200); }"
            "\nQMenu::item:disabled:selected { background-color: rgb(236, 234, 242);"
            " color: rgb(175, 175, 188); }"
            "\nQMenu::separator { height: 1px; background: rgb(210, 214, 226);"
            " margin: 4px 8px; }"
            # The two panels: each a box of its own, on a slightly darker
            # ground so the gap between them reads as a gap.
            "\nQWidget#centralwidget { background-color: rgb(232, 228, 240); }"
            "\nQFrame#ToolPanel, QFrame#MainPanel {"
            " border: 1px solid rgb(200, 204, 222); border-radius: 6px; }"
            "\nQFrame#ToolPanel QGroupBox { font-weight: bold; }"
            "\nQFrame#ToolPanel QPushButton { text-align: left;"
            " padding: 3px 8px; font-weight: normal; }"
            "\nQLabel#ToolPanelTitle { font-weight: bold; }"
        )

        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        MainWindow.setCentralWidget(self.centralwidget)
        window_layout = QtWidgets.QHBoxLayout(self.centralwidget)
        window_layout.setContentsMargins(6, 6, 6, 6)
        window_layout.setSpacing(8)

        self.ToolPanel = self._tool_panel(self.centralwidget)
        window_layout.addWidget(self.ToolPanel)

        self.MainPanel = QtWidgets.QFrame(self.centralwidget)
        self.MainPanel.setObjectName("MainPanel")
        window_layout.addWidget(self.MainPanel, stretch=1)
        menu_layout = QtWidgets.QVBoxLayout(self.MainPanel)

        # -- folder + explicit file selection --
        folder_row = QtWidgets.QHBoxLayout()
        self.FilePathLE = QtWidgets.QLineEdit(self.centralwidget)
        self.FilePathLE.setPlaceholderText("Folder containing .nxs files")
        self.BrowseFolderButton = QtWidgets.QPushButton("Folder...", self.centralwidget)
        folder_row.addWidget(self.FilePathLE, stretch=1)
        folder_row.addWidget(self.BrowseFolderButton)
        menu_layout.addLayout(folder_row)

        pick_row = QtWidgets.QHBoxLayout()
        # Loading has its own window (ui.launcher.loader_dialog): which beamline wrote
        # the file, what order its axes are in and what its scanned axis
        # physically is are all questions that belong at the door, and there
        # is no room for them on a button.
        self.LoadDataButton = QtWidgets.QPushButton("Load data...", self.centralwidget)
        self.LoadDataButton.setToolTip(
            "Open the loader: choose files, choose which beamline's layout "
            "to read them as, and set the axis options before they are "
            "added to the list.")
        self.PickFilesButton = QtWidgets.QPushButton("Quick add...", self.centralwidget)
        self.PickFilesButton.setToolTip(
            "Add .nxs files straight to the list, with the reader detected "
            "and the axes left as the file has them.")
        self.ClearListButton = QtWidgets.QPushButton("Clear list", self.centralwidget)
        pick_row.addWidget(self.LoadDataButton, stretch=1)
        pick_row.addWidget(self.PickFilesButton, stretch=1)
        pick_row.addWidget(self.ClearListButton)
        menu_layout.addLayout(pick_row)

        self.OpenHint = QtWidgets.QLabel(
            "<b>Double-click a file to open it</b> in its own window; right-click "
            "for its cuts, saving, renaming and its information. The tools on "
            "the left light up for whatever the selection can take.",
            self.centralwidget)
        self.OpenHint.setWordWrap(True)
        menu_layout.addWidget(self.OpenHint)

        # A plain list, not icons: the data kind is spelled out after each
        # name, which is more informative than a thumbnail and costs nothing
        # to produce -- building thumbnails meant fully loading every file.
        self.FilePathListWidget = QtWidgets.QListWidget(self.centralwidget)
        self.FilePathListWidget.setViewMode(QtWidgets.QListView.ListMode)
        self.FilePathListWidget.setAlternatingRowColors(True)
        # Several rows at a time, so deleting or saving a batch is one action.
        self.FilePathListWidget.setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection)
        self.FilePathListWidget.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.FilePathListWidget.setToolTip(
            "Double-click to open. Right-click to rename, save or remove the "
            "selected datasets; Ctrl/Shift-click selects several.")
        menu_layout.addWidget(self.FilePathListWidget, stretch=1)

        # -- what can be done to the data itself: one button per operation,
        # each producing a new dataset that is the old one combined, cut
        # down, binned or rescaled. Analysis and filters are the tool
        # panel's; these change the data and nothing else. Greyed out until
        # the selection is something they take (see list_actions).
        list_ops = QtWidgets.QHBoxLayout()
        list_ops.setSpacing(4)
        self.ListButtons = {}
        for name, text, tooltip in LIST_OPERATIONS:
            button = QtWidgets.QPushButton(text, self.centralwidget)
            button.setToolTip(tooltip)
            button.setProperty("description", tooltip)
            list_ops.addWidget(button, stretch=1)
            self.ListButtons[name] = button
        menu_layout.addLayout(list_ops)

        # -- what the selected data is: extent, point count and step per
        # axis. It stays in the launcher because it is what you compare
        # between files while browsing; the rest of the information (the
        # metadata table) is a window of its own. Filled on demand
        # (right-click -> Show information), since reading a dataset to
        # describe it is work.
        self.DataInfoBox = QtWidgets.QGroupBox("Data Information", self.centralwidget)
        info_grid = QtWidgets.QGridLayout(self.DataInfoBox)
        info_grid.setContentsMargins(6, 4, 6, 4)
        info_grid.setSpacing(3)
        self.DataInfoHeaders = []
        self.DataInfoCells = {}
        for column, name in enumerate(("X", "Y", "Z", "E")):
            header = QtWidgets.QLabel(name)
            header.setAlignment(QtCore.Qt.AlignCenter)
            header_font = header.font()
            header_font.setBold(True)
            header.setFont(header_font)
            info_grid.addWidget(header, 0, column + 1)
            self.DataInfoHeaders.append(header)
        for row, name in enumerate(("min", "max", "num", "step")):
            info_grid.addWidget(QtWidgets.QLabel(name), row + 1, 0)
            for column in range(4):
                cell = QtWidgets.QLineEdit()
                cell.setReadOnly(True)
                cell.setAlignment(QtCore.Qt.AlignCenter)
                cell.setStyleSheet("background-color: white;")
                info_grid.addWidget(cell, row + 1, column + 1)
                self.DataInfoCells[(name, column)] = cell
        menu_layout.addWidget(self.DataInfoBox)

        # -- which build this is, and where to report a bug or get a newer
        # one. Updates come often; without this a user cannot tell which
        # version they are running, or whether a fix has reached them.
        self.AboutLabel = QtWidgets.QLabel(
            f"<b>{version.title()}</b><br>"
            f"Bug reports, questions and the latest version: "
            f"<a href=\"{version.REPOSITORY}\">{version.REPOSITORY}</a>"
            f" (bugs: <a href=\"{version.REPOSITORY}/issues\">Issues</a>)",
            self.centralwidget)
        self.AboutLabel.setWordWrap(True)
        self.AboutLabel.setOpenExternalLinks(True)
        self.AboutLabel.setTextInteractionFlags(QtCore.Qt.TextBrowserInteraction)
        self.AboutLabel.setToolTip(
            "Please give this version number when reporting a bug.")
        menu_layout.addWidget(self.AboutLabel)

        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        MainWindow.setStatusBar(self.statusbar)


    def _tool_panel(self, parent):
        """The tool panel: the tools for the selected datasets, then the
        memory controls at the bottom.

        One button per tool, so every one of them is in sight -- they used
        to be behind a right-click or on a tab of a window the user had to
        know to open. Grouped by the dimension of data a tool takes, and each
        greyed out until the list's selection is something it can take; its
        tooltip then says why (ARPES_viewer.update_tool_buttons()).
        """
        panel = QtWidgets.QFrame(parent)
        panel.setObjectName("ToolPanel")
        panel.setFixedWidth(228)
        layout = QtWidgets.QVBoxLayout(panel)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(6)

        title = QtWidgets.QLabel("Tools")
        title.setObjectName("ToolPanelTitle")
        title.setAlignment(QtCore.Qt.AlignCenter)
        layout.addWidget(title)

        self.ToolButtons = {}
        self.ToolGroups = {}
        for heading, tools in TOOL_GROUPS:
            group = QtWidgets.QGroupBox(heading, panel)
            box = QtWidgets.QVBoxLayout(group)
            box.setContentsMargins(4, 4, 4, 4)
            box.setSpacing(2)
            for name, text, tooltip in tools:
                button = QtWidgets.QPushButton(text, group)
                button.setToolTip(tooltip)
                button.setProperty("description", tooltip)
                box.addWidget(button)
                self.ToolButtons[name] = button
            layout.addWidget(group)
            self.ToolGroups[heading] = group
        layout.addStretch(1)

        # -- memory: how much this program is holding, and a way to give
        # some of it back. Several spatial scans open at once is where it
        # adds up, so the numbers are always on show rather than in a menu.
        # Refreshed by a timer in ARPES_viewer.py.
        self.MemoryBox = QtWidgets.QGroupBox("Memory", panel)
        memory_layout = QtWidgets.QVBoxLayout(self.MemoryBox)
        memory_layout.setContentsMargins(4, 4, 4, 4)
        memory_layout.setSpacing(3)
        self.MemoryUsageBar = QtWidgets.QProgressBar()
        self.MemoryUsageBar.setRange(0, 100)
        self.MemoryUsageBar.setTextVisible(True)
        self.MemoryUsageBar.setFormat("n/a")
        self.MemoryUsageBar.setToolTip(
            "This program's resident memory, as a share of the machine's RAM.")
        memory_layout.addWidget(self.MemoryUsageBar)
        self.MemoryLabel = QtWidgets.QLabel("")
        self.MemoryLabel.setWordWrap(True)
        self.MemoryLabel.setStyleSheet("font-weight: normal;")
        memory_layout.addWidget(self.MemoryLabel)
        self.FreeMemoryButton = QtWidgets.QPushButton("Free memory")
        self.FreeMemoryButton.setToolTip(
            "Drop cached copies of computed datasets (they are re-read from "
            "their auto-saved files when next opened) and give freed memory "
            "back to the system. Open windows are not touched.")
        self.CloseViewersButton = QtWidgets.QPushButton("Close all viewers")
        self.CloseViewersButton.setToolTip(
            "Close every open viewer window, then free memory. The file list "
            "is kept, so anything can be reopened.")
        # -- the auto-saved working copies, for clearing out by hand.
        self.OpenAutosaveButton = QtWidgets.QPushButton("Open auto-save folder")
        self.OpenAutosaveButton.setToolTip(
            "Open the folder where computed datasets are auto-saved, in the "
            "system's file browser, to delete old working copies by hand. "
            "The folder of the session now running is in use; leave it.")
        for button in (self.FreeMemoryButton, self.CloseViewersButton,
                       self.OpenAutosaveButton):
            memory_layout.addWidget(button)
        layout.addWidget(self.MemoryBox)
        return panel


#: The operations under the main list: (action, button text, tooltip). The
#: action names are list_actions' and ARPES_viewer.run_list_action()'s.
LIST_OPERATIONS = (
    ("combine", "Combine...",
     "Join or stack the selected datasets into one: cuts into a map along an "
     "axis of your own (temperature, gate voltage, photon energy...), or the "
     "pieces of a map or kz map measured in segments back into one. Two or "
     "more datasets of one kind."),
    ("truncate", "Truncate...",
     "Crop the selected datasets to a range on each axis, by value or by "
     "index. Several at once, if they have the same kind and axis lengths."),
    ("compress", "Compress...",
     "Bin the selected datasets: sum every N points of an axis into one, "
     "which keeps the counts and improves the noise."),
    ("normalise", "Self-normalise...",
     "Divide each spectrum (each EDC, each cut...) by its own total or peak, "
     "optionally over a window of your own."),
)

#: The tool panel: (group heading, [(action, button text, tooltip), ...]).
#: Grouped by the dimension of the data a tool takes -- 1D curves (EDC,
#: MDC), 2D cuts, 3D maps (angle maps, kz maps, k-maps).
TOOL_GROUPS = (
    ("1D && 2D \u2014 curves, cuts", (
        ("figure", "Figure composer...",
         "One panel per selected dataset, in the figure composer."),
    )),
    ("2D && 3D \u2014 cuts, maps", (
        ("smooth", "Smooth...",
         "Gaussian, Savitzky-Golay or moving-average smoothing, with spike "
         "removal. A map is smoothed plane by plane, on the plane you pick."),
        ("derivative", "Second derivative...",
         "-d\u00b2I/dE\u00b2 or -d\u00b2I/dk\u00b2 after a pre-smooth, to "
         "sharpen the bands."),
        ("curvature", "Curvature...",
         "The 1D or 2D curvature (Zhang et al.), which follows the band "
         "better than a second derivative."),
        ("symmetry", "Symmetrise...",
         "Mirror a cut about a line; for a map, mirror or rotate its "
         "constant-energy planes."),
        ("degrid", "Remove detector grid...",
         "Remove the detector's grid pattern. From a whole map (or kz map), "
         "using the map itself as the reference -- do this first, before a "
         "k conversion or Fermi-surface correction. From a cut, with a grid "
         "found on a map with the same settings if one is in the list, by "
         "notching the cut's own grid peaks otherwise."),
    )),
    ("2D only \u2014 cuts", (
        ("background", "Background subtraction...",
         "Remove the intensity that is not the bands: Shirley or Tougaard "
         "inelastic background, a floor common to every angle, a smooth "
         "polynomial -- or divide out the Fermi cut-off."),
        ("stack", "Stack plot (EDC / MDC)...",
         "A waterfall of EDCs or MDCs across the cut."),
        ("mdcfit", "MDC fit (k-cut only)...",
         "Fit the momentum distribution curves line by line, then read the "
         "band off the fitted centres: Fermi velocity, effective mass, "
         "self-energy. Only for a cut converted to momentum (k)."),
        ("edcfit", "EDC fit...",
         "Fit the energy distribution curves line by line, optionally with "
         "the Fermi cut-off in the model. Any cut; v_F and m* from the "
         "fitted band need a cut in momentum."),
        ("corelevel", "Core-level check...",
         "Which element and orbital each peak of the angle-integrated EDC is: "
         "the sample's elements against tabulated binding energies, after "
         "checking E_F and the photon energy; unexplained peaks get possible "
         "contaminants."),
        ("arithmetic", "Cut arithmetic (two cuts)...",
         "Combine two cuts: linear or circular dichroism, dividing by a "
         "reference, A \u2212 B, A / B, (A \u2212 B)/(A + B), A + B. The first "
         "selected is A; the window can swap them."),
    )),
    ("3D only \u2014 maps", (
        ("slices3d", "3D slices...",
         "The three orthogonal planes through the cube, each placed where "
         "you choose, each with its own colours."),
        ("notch3d", "Notched cube...",
         "The cube with one corner cut away: a constant-energy map and two "
         "dispersions in one picture."),
        ("projection3d", "Volume projection...",
         "The whole volume projected along the view: maximum intensity, "
         "sum, alpha composite, or a shaded isosurface."),
    )),
)


class InfoWindow(QtWidgets.QWidget):
    """The selected dataset's metadata, in a window of its own.

    Opened from the list's right-click menu ("Show information"). The
    metadata table is long -- a hundred-odd fields -- and under the file
    list it left neither the list nor itself much room; in its own window it
    can be moved beside a viewer, kept open while other files are browsed,
    or closed and forgotten. The short summary of what the data *is* stays
    in the launcher (``Ui_MainWindow.DataInfoBox``), where it is what you
    compare between files.
    """

    def __init__(self, parent=None):
        # Qt.Window, not a plain child: a QWidget given a parent without it
        # is laid out *inside* that parent rather than standing on its own.
        super().__init__(parent, QtCore.Qt.Window)
        self.setWindowTitle("Dataset information")
        self.resize(560, 620)
        layout = QtWidgets.QVBoxLayout(self)

        self.TitleLabel = QtWidgets.QLabel("No dataset")
        self.TitleLabel.setAlignment(QtCore.Qt.AlignCenter)
        font = self.TitleLabel.font()
        font.setBold(True)
        self.TitleLabel.setFont(font)
        layout.addWidget(self.TitleLabel)

        self.KindLabel = QtWidgets.QLabel("No file selected")
        self.KindLabel.setAlignment(QtCore.Qt.AlignCenter)
        layout.addWidget(self.KindLabel)

        # -- the file's own metadata: a different question (photon energy,
        # temperature, slit), so it keeps its filter and its CSV export --
        self.InfoFilterLE = QtWidgets.QLineEdit()
        self.InfoFilterLE.setPlaceholderText("Filter metadata (e.g. 'mono', 'MBS')")
        layout.addWidget(self.InfoFilterLE)

        self.FileInfoTW = QtWidgets.QTableWidget()
        self.FileInfoTW.setColumnCount(2)
        self.FileInfoTW.setHorizontalHeaderLabels(["Field", "Value"])
        self.FileInfoTW.horizontalHeader().setStretchLastSection(True)
        self.FileInfoTW.verticalHeader().setVisible(False)
        self.FileInfoTW.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.FileInfoTW.setAlternatingRowColors(True)
        layout.addWidget(self.FileInfoTW, stretch=1)

        self.ExportMetaButton = QtWidgets.QPushButton("Export metadata (.csv)")
        self.ExportMetaButton.setToolTip("Write the metadata table to CSV.")
        layout.addWidget(self.ExportMetaButton)
