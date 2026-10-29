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

from arpes_viewer import version
from arpes_viewer.ui.common.platform_style import MOD_KEY, ui_font_family


class Ui_MainWindow:
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(840, 940)
        MainWindow.setWindowTitle(version.title())
        # The panel's colours. The disabled rules are not optional: a colour
        # set for every widget is also the colour of a *disabled* widget and
        # of a disabled menu entry, so without them a greyed-out entry of the
        # list's right-click menu looked exactly like an enabled one.
        MainWindow.setStyleSheet(
            "* { background-color: rgb(245, 241, 249); "
            f"font-family: \"{ui_font_family()}\"; "
            "color: rgb(102, 126, 161); }"
            "\n*:disabled { color: rgb(190, 190, 200); }"
            "\nQMenu::item { padding: 3px 24px 3px 20px; }"
            "\nQMenu::item:selected { background-color: rgb(214, 224, 240); "
            "color: rgb(40, 60, 95); }"
            "\nQMenu::item:disabled { color: rgb(190, 190, 200); }"
            "\nQMenu::item:disabled:selected { background-color: rgb(236, 234, 242); "
            "color: rgb(175, 175, 188); }"
            "\nQMenu::separator { height: 1px; background: rgb(210, 214, 226); "
            "margin: 4px 8px; }"
            # The two panels: each a box of its own, on a slightly darker
            # ground so the gap between them reads as a gap.
            "\nQWidget#centralwidget { background-color: rgb(232, 228, 240); }"
            "\nQFrame#ToolPanel, QFrame#MainPanel { "
            "border: 1px solid rgb(200, 204, 222); border-radius: 6px; }"
            "\nQFrame#ToolPanel QGroupBox { font-weight: bold; }"
            # The tool panel uses the main panel's font, size included; only
            # its headings are set apart, in grey italics.
            "\nQLabel#ToolHeading { font-weight: normal; font-style: italic;"
            " color: rgb(120, 130, 160); padding-top: 3px; }"
            "\nQFrame#ToolPanel QPushButton { text-align: left; "
            "padding: 3px 8px; font-weight: normal; }"
            "\nQTabWidget#ToolTabs QTabBar::tab { padding: 3px 10px;"
            " min-width: 60px; }"
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
            "to open a map's cuts. The tools on the left -- and its Datasets "
            "box: information, the experiment log, renaming, saving -- light "
            "up for whatever the selection can take.",
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
            "Double-click to open; right-click for a map's slit and scan-axis "
            f"cuts. {MOD_KEY}/Shift-click selects several.")
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
        # (Datasets -> Experiment info), since reading a dataset to
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
            f" (bugs: <a href=\"{version.REPOSITORY}/issues\">Issues</a>)<br>"
            f"With a bug report, upload the error log (attach the file) or "
            f"paste its text: <a href=\"errors:open\">Open the error logs</a>"
            f" · <a href=\"errors:copy\">Copy the newest error log</a>",
            self.centralwidget)
        self.AboutLabel.setWordWrap(True)
        self.AboutLabel.setOpenExternalLinks(False)
        self.AboutLabel.linkActivated.connect(
            lambda link: about_link(link, MainWindow))
        self.AboutLabel.setTextInteractionFlags(QtCore.Qt.TextBrowserInteraction)
        self.AboutLabel.setToolTip(
            "Please give this version number when reporting a bug. Errors are "
            "saved by themselves in error_logs/ in the program's folder: "
            "upload that file with the report (GitHub Issues: drag it into "
            "the text box), or copy its text and paste it there or into an "
            "e-mail.")
        menu_layout.addWidget(self.AboutLabel)

        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        MainWindow.setStatusBar(self.statusbar)


    def _tool_panel(self, parent):
        """The tool panel: the tools for the selected datasets, then what
        is done to the rows themselves (the Datasets box), then the memory
        controls at the bottom.

        One button per tool, so every one of them is in sight -- they used
        to be behind a right-click or on a tab of a window the user had to
        know to open. On tabs, Basic (shown first), Advanced and Calculation
        (TOOL_SECTIONS), each grouped by the dimension of data a tool
        takes (a third, Calculation, for what is not a dataset), and each
        button
        greyed out until the list's selection is something it can take; its
        tooltip then says why (ARPES_viewer.update_tool_buttons()).
        """
        panel = QtWidgets.QFrame(parent)
        panel.setObjectName("ToolPanel")
        layout = QtWidgets.QVBoxLayout(panel)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)

        title = QtWidgets.QLabel("Tools")
        title.setObjectName("ToolPanelTitle")
        title.setAlignment(QtCore.Qt.AlignCenter)
        layout.addWidget(title)

        # Three tabs, Basic, Advanced and Calculation, with Basic up front: the steps nearly
        # every dataset goes through, without the analysis tools crowding
        # them. In each, a small heading per dimension of data and the
        # buttons under it, in a scroll area so a short screen scrolls the
        # tools rather than squeezing them -- the memory box below stays put.
        self.ToolTabs = QtWidgets.QTabWidget(panel)
        self.ToolTabs.setObjectName("ToolTabs")
        self.ToolButtons = {}
        self.ToolSections = {}
        for section, groups in TOOL_SECTIONS:
            page = QtWidgets.QWidget()
            page.setObjectName("ToolSection")
            column = QtWidgets.QVBoxLayout(page)
            column.setContentsMargins(4, 4, 4, 4)
            column.setSpacing(2)
            for heading, names in groups:
                label = QtWidgets.QLabel(heading.replace("&&", "&"), page)
                label.setObjectName("ToolHeading")
                column.addWidget(label)
                for name in names:
                    text, tooltip = TOOLS[name]
                    button = QtWidgets.QPushButton(text, page)
                    button.setToolTip(tooltip)
                    button.setProperty("description", tooltip)
                    column.addWidget(button)
                    self.ToolButtons[name] = button
            column.addStretch(1)
            scroller = QtWidgets.QScrollArea()
            scroller.setWidgetResizable(True)
            scroller.setFrameShape(QtWidgets.QFrame.NoFrame)
            scroller.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
            scroller.setWidget(page)
            self.ToolTabs.addTab(scroller, section)
            self.ToolSections[section] = page
        self.ToolTabs.setCurrentIndex(0)
        layout.addWidget(self.ToolTabs, stretch=1)
        # 262 px fits the longest button in the Windows font; a wider font
        # (macOS's) widens the panel rather than cutting the labels short.
        # The extra is the scroll bar, the tab frame and the margins.
        widest = max(page.sizeHint().width()
                     for page in self.ToolSections.values())
        scrollbar = panel.style().pixelMetric(
            QtWidgets.QStyle.PM_ScrollBarExtent)
        panel.setFixedWidth(max(262, widest + scrollbar + 24))

        # -- what is done to the rows themselves rather than to their data:
        # describing, logging, naming, saving and removing them. These were
        # the right-click menu's, which kept only the three ways of opening
        # a row; here they are in sight, and grey out like the tools.
        self.DatasetBox = QtWidgets.QGroupBox("Datasets", panel)
        dataset_layout = QtWidgets.QVBoxLayout(self.DatasetBox)
        dataset_layout.setContentsMargins(4, 4, 4, 4)
        dataset_layout.setSpacing(2)
        self.DatasetButtons = {}
        actions = {name: (text, tip) for name, text, tip in DATASET_ACTIONS}
        for row_names in DATASET_ROWS:
            row = QtWidgets.QHBoxLayout()
            row.setSpacing(2)
            for name in row_names:
                text, tooltip = actions[name]
                button = QtWidgets.QPushButton(text, self.DatasetBox)
                button.setToolTip(tooltip)
                button.setProperty("description", tooltip)
                row.addWidget(button, stretch=1)
                self.DatasetButtons[name] = button
            dataset_layout.addLayout(row)
        layout.addWidget(self.DatasetBox)

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
        # Three short buttons on one row: the box sits under the tools and
        # every line it takes is a line of tools scrolled out of sight.
        self.FreeMemoryButton = QtWidgets.QPushButton("Free RAM")
        self.FreeMemoryButton.setToolTip(
            "Free memory: drop cached copies of computed datasets (they are "
            "re-read from their auto-saved files when next opened) and give "
            "freed memory back to the system. Open windows are not touched.")
        self.CloseViewersButton = QtWidgets.QPushButton("Close all")
        self.CloseViewersButton.setToolTip(
            "Close all viewers: close every open viewer window, then free "
            "memory. The file list is kept, so anything can be reopened.")
        # -- the auto-saved working copies, for clearing out by hand.
        self.OpenAutosaveButton = QtWidgets.QPushButton("Autosaves")
        self.OpenAutosaveButton.setToolTip(
            "Open the auto-save folder, where computed datasets are written, "
            "in the system's file browser, to delete old working copies by "
            "hand. The folder of the session now running is in use; leave it.")
        memory_row = QtWidgets.QHBoxLayout()
        memory_row.setSpacing(2)
        for button in (self.FreeMemoryButton, self.CloseViewersButton,
                       self.OpenAutosaveButton):
            button.setStyleSheet("text-align: center; padding: 3px 4px;")
            memory_row.addWidget(button, stretch=1)
        memory_layout.addLayout(memory_row)
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

#: The tool panel's Datasets box, top to bottom: (action, button text,
#: tooltip). The action names are list_actions' and
#: ARPES_viewer.run_list_action()'s.
DATASET_ACTIONS = (
    ("info", "Experiment info",
     "The selected dataset's kind, axes and metadata: the photon energy, "
     "temperature, analyser settings, motor positions... in a window of "
     "its own, with a filter and a CSV export."),
    ("explogformat", "Log format...",
     "Set the experiment log up, once: its folder and file name, its "
     "separator, and its columns -- chosen from a dataset's fields and put "
     "in the order you like. An existing log shows its own format, to "
     "adjust."),
    ("explog", "Auto Exp. log",
     "Add the selected datasets to the experiment log (.csv), one row each, "
     "in its format, after its last record. The first time, the format is "
     "asked for."),
    ("rename", "Rename",
     "Rename: change the name the dataset is listed (and saved) under. The "
     "file on disk keeps its name."),
    ("save", "Save",
     "Save the selected datasets: to a .nxs file of your own, which this "
     "program opens again, or for Igor Pro (.itx text, or one .ibw binary "
     "wave per dataset), axes, units and metadata included."),
    ("remove", "Remove",
     "Remove the selected rows from the list (after asking). Files on disk "
     "are not deleted."),
    ("sumrepeat", "Sum repeated scans...",
     "Add up a measurement taken in several batches -- the same fine cut or "
     "fine map acquired again and again, so a slowly ageing sample can be "
     "seen and its late scans left out. Checks the scans really are "
     "repeats, shows how each one differs from the rest, and can line up a "
     "drifting energy scale before summing."),
    ("log", "Session log...",
     "What has been auto-saved, saved and deleted this session, and where "
     "each of it is -- the file to look at after a crash."),
)

#: How the Datasets box lays them out: one row each, but the log's format
#: and the log itself side by side, and the three short row operations --
#: Rename, Save, Remove -- on one line.
DATASET_ROWS = (("info",), ("explogformat", "explog"),
                ("rename", "save", "remove"), ("sumrepeat",), ("log",))

#: Every tool of the tool panel: action -> (button text, tooltip). The
#: action names are list_actions' and ARPES_viewer.run_list_action()'s.
TOOLS = {
    "figure": ("Figure composer...",
             "A figure laid out like a manuscript page: one panel per "
             "selected dataset, or, with nothing selected, a blank page of "
             "empty panels to fill from the list or with pictures."),
    "fermi": ("Fermi edge fit...",
             "Fit the Fermi edge: E_F and the resolution, then offset the "
             "energy axis to E - E_F. A cut opens with its viewer; a map "
             "with its slit cut, whose slider picks the cut to fit -- one "
             "E_F for the whole map, or one per cut, interpolated between "
             "them."),
    "smooth": ("Smooth...",
             "Gaussian, Savitzky-Golay or moving-average smoothing, with "
             "spike removal. A map is smoothed plane by plane, on the plane "
             "you pick."),
    "derivative": ("Second derivative...",
             "-d\xb2I/dE\xb2 or -d\xb2I/dk\xb2 after a pre-smooth, to "
             "sharpen the bands."),
    "curvature": ("Curvature...",
             "The 1D or 2D curvature (Zhang et al.), which follows the band "
             "better than a second derivative."),
    "symmetry": ("Symmetrise...",
             "Mirror a cut about a line; for a map, mirror or rotate its "
             "constant-energy planes."),
    "degrid": ("Remove detector grid...",
             "Remove the detector's grid pattern. From a whole map (or kz "
             "map), using the map itself as the reference -- do this first, "
             "before a k conversion or Fermi-surface correction. From a cut, "
             "with a grid found on a map with the same settings if one is "
             "in the list, by notching the cut's own grid peaks otherwise."),
    "background": ("Background subtraction...",
             "Remove the intensity that is not the bands: Shirley or "
             "Tougaard inelastic background, a floor common to every angle, "
             "a smooth polynomial -- or divide out the Fermi cut-off."),
    "stack": ("Stack plot (EDC / MDC)...",
             "A waterfall of EDCs or MDCs across the cut."),
    "mdcfit": ("MDC fit (k-cut only)...",
             "Fit the momentum distribution curves line by line, then read "
             "the band off the fitted centres: Fermi velocity, effective "
             "mass, self-energy. Only for a cut converted to momentum (k)."),
    "edcfit": ("EDC fit...",
             "Fit the energy distribution curves line by line, optionally "
             "with the Fermi cut-off in the model. Any cut; v_F and m* from "
             "the fitted band need a cut in momentum."),
    "corelevel": ("Core-level check...",
             "Which element and orbital each peak of the angle-integrated "
             "EDC is: the sample's elements against tabulated binding "
             "energies, after checking E_F and the photon energy; "
             "unexplained peaks get possible contaminants."),
    "arithmetic": ("Cut arithmetic (two cuts)...",
             "Combine two cuts: linear or circular dichroism, dividing by a "
             "reference, A \u2212 B, A / B, (A \u2212 B)/(A + B), A + B. The "
             "first selected is A; the window can swap them."),
    "kzalign": ("Align Fermi levels (kz)...",
             "Fit the Fermi edge of every spectrum of a photon-energy scan, "
             "shift each one to put its own edge at zero, and crop to the "
             "energy range they all still cover. Before the kz conversion."),
    "slices3d": ("3D slices...",
             "The three orthogonal planes through the cube drawn in 3-D, "
             "turned with the mouse, each placed where you choose and with "
             "its own colours. For flat panels on a page, see Slice series "
             "figure."),
    "notch3d": ("Notched cube...",
             "The cube with one corner cut away: a constant-energy map and "
             "two dispersions in one picture."),
    "stack3d": ("Stacked slices...",
             "Parallel slices of the cube stacked in 3-D -- constant-energy "
             "maps one above another up the energy axis, or cuts one behind "
             "another along k_z -- each layer's background transparent, "
             "with a curve marked through the layers (a Dirac point "
             "followed). For figures of how a band evolves."),
    "projection3d": ("Volume projection...",
             "The whole volume projected along the view: maximum intensity, "
             "sum, alpha composite, or a shaded isosurface."),
    "kmap": ("k conversion (map)...",
             "Convert an angle map to momentum (\u00c5\u207b\u00b9). It opens on "
             "its constant-energy contour, where the readout cursor picks "
             "\u0393; the settings start from the last conversion of a map "
             "of the same size."),
    "kcut": ("k conversion (cut)...",
             "Convert a cut to momentum. By default a vertical line on the "
             "cut picks \u0393's slit angle (\u0393 taken to lie on the "
             "cut); Advanced inherits \u0393 from a converted map. Starts "
             "from the last conversion of a cut of the same size."),
    "kzconv": ("kz conversion...",
             "Convert a photon-energy scan to (k_z, k_par), choosing the "
             "inner potential against the lattice's own period. Flatten the "
             "Fermi edge first. Starts from the last conversion of a scan of "
             "the same size."),
    "slicefig": ("Slice series figure...",
                 "A page of 2-D slices through the map as one figure for a "
                 "paper: constant-energy contours at a series of energies, or "
                 "cuts at a series of angles or momenta, with a shared colour "
                 "scale and labels on the outer edges only. Flat panels, not "
                 "a 3-D view (that is 3D slices)."),
    "ridge": ("Ridge detection...",
              "Trace the bands as ridges of the intensity: the Hessian of "
              "the smoothed image, kept where it curves down across the band "
              "and the slope across it changes sign. Sharper than a second "
              "derivative, with no false double lines; the ridge points can "
              "be listed or saved for a fit. A map is done plane by plane."),
    "tbfit": ("Tight-binding fit...",
              "Fit a tight-binding band E(k) = e0 + sum of hoppings x "
              "cos(k.R) to band points -- picked on a k-map, traced by ridge "
              "detection, or read from a file -- shell by shell of "
              "neighbours of a 2-D or 3-D lattice. Step by step; the model's "
              "contours go back onto the map and its band into the list."),
    "detcorr": ("Detector corrections...",
                "Before anything else, on the detector's own pixels: make a "
                "trapezoidal (keystone) image rectangular again, and patch "
                "dead or bad regions of the detector from their "
                "neighbours. A cut, or every slice of a map at once."),
    "lattice": ("Crystal structure...",
                "The crystal in real space, in 3-D: typed in, taken from the "
                "Brillouin-zone tool's last crystal, or read from a .cif "
                "file. Draw a cleavage plane (hkl) and slide it through the "
                "cell; atoms' colours, gloss, transparency and size are "
                "yours to set. The picture goes to the figure composer."),
    "bands": ("Calculated bands...",
              "A band structure from the theory side -- a .bxsf band grid "
              "(Wien2k, QE, VASP, Wannier90), a Wannier90 _hr.dat "
              "Hamiltonian, or bands along a path (_band.dat, "
              ".spaghetti_ene, VASP EIGENVAL) -- along a high-symmetry path "
              "over the measurement assembled from the selected k-maps, and "
              "over a cut or a map where it was really measured (its k_z). "
              "Fits the renormalisation Z and the shift; panels go to the "
              "figure composer. Select k-maps, kz maps in k and cuts in k."),
    "bz": ("Brillouin zone...",
           "Overlay a Brillouin zone (or a moire zone) on a map in momentum "
           "-- a k-map or a converted kz map. Once plotted, the zone stays "
           "with the dataset: its viewer shows it with the Brillouin zone "
           "tick box, and its exports can draw it."),
}

#: The tool panel's layout: three tabs -- Basic, the steps nearly every
#: dataset goes through (E_F, the detector grid, the background, the k
#: conversions, the core levels), Advanced, the analysis that follows, and
#: Calculation, what is worked out rather than measured (the crystal) --
#: each grouped by the dimension of data a tool takes: 1D curves (EDC, MDC),
#: 2D cuts, 3D maps (angle maps, kz maps, k-maps).
TOOL_SECTIONS = (
    ("Basic", (
        ("2D && 3D \u2014 cuts, maps", ("fermi", "degrid", "detcorr")),
        ("2D only \u2014 cuts", ("background", "kcut", "corelevel")),
        ("3D only \u2014 maps", ("kmap", "kzalign", "kzconv", "bz")),
    )),
    ("Advanced", (
        ("1D && 2D \u2014 curves, cuts", ("figure",)),
        ("2D && 3D \u2014 cuts, maps", ("smooth", "derivative", "curvature",
                                       "symmetry", "ridge")),
        ("2D only \u2014 cuts", ("stack", "mdcfit", "edcfit", "arithmetic")),
        ("3D only \u2014 maps", ("slicefig", "slices3d", "notch3d",
                                  "projection3d", "stack3d")),
    )),
    # What is worked out rather than measured: the crystal itself, band
    # models fitted to the data, and calculated band structures.
    ("Calculation", (
        ("Crystal \u2014 real space", ("lattice",)),
        ("Band models \u2014 fitted to the data", ("tbfit",)),
        ("Band structures \u2014 calculated", ("bands",)),
    )),
)


class InfoWindow(QtWidgets.QWidget):
    """The selected dataset's metadata, in a window of its own.

    Opened from the tool panel's Datasets box ("Experiment info"). The
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
        self.setWindowTitle("Experiment info")
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


def about_link(link: str, window=None):
    """A link of the main panel's foot: the error logs (open their folder,
    or copy the newest one's text to paste into a bug report), else a web
    page."""
    from PyQt5.QtCore import QUrl
    from PyQt5.QtGui import QDesktopServices
    from PyQt5.QtWidgets import QApplication

    from arpes_viewer.tools.runtime import errorlog

    def say(text):
        bar = window.statusBar() if window is not None else None
        if bar is not None:
            bar.showMessage(text, 12000)

    if link == "errors:open":
        folder = errorlog.log_folder()
        if folder:
            QDesktopServices.openUrl(QUrl.fromLocalFile(folder))
            say(f"Error logs: {folder}")
        return folder
    if link == "errors:copy":
        path = errorlog.newest()
        if not path:
            say("No error log yet: nothing has gone wrong.")
            return ""
        with open(path, encoding="utf-8", errors="replace") as handle:
            QApplication.clipboard().setText(handle.read())
        say(f"Copied {path} -- paste it into the bug report "
            "(or attach the file itself).")
        return path
    QDesktopServices.openUrl(QUrl(link))
    return link
