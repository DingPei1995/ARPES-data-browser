"""
ui/main_window.py
==============
The launcher window's layout: a file picker, the global colormap controls,
and the metadata table for whichever file is selected. Hand-written PyQt5,
no Qt-Designer ``.ui``/``.qrc`` needed.

Visualisation is NOT here -- double-clicking a file in the list opens a
viewer window from ``ui/windows.py``. Keeping the launcher small
means several files (and several cuts of one file) can be open side by side
and arranged freely, which the single stacked layout this replaced could not
do.

Wiring lives in ``ARPES_viewer.py`` (``SetConnect()``), the same split the
original app used.
"""
from PyQt5 import QtCore, QtWidgets


class Ui_MainWindow:
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(560, 900)
        MainWindow.setWindowTitle("ARPES viewer")
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
        )

        self.centralwidget = QtWidgets.QWidget(MainWindow)
        MainWindow.setCentralWidget(self.centralwidget)
        menu_layout = QtWidgets.QVBoxLayout(self.centralwidget)

        # -- folder + explicit file selection --
        folder_row = QtWidgets.QHBoxLayout()
        self.FilePathLE = QtWidgets.QLineEdit(self.centralwidget)
        self.FilePathLE.setPlaceholderText("Folder containing .nxs files")
        self.BrowseFolderButton = QtWidgets.QPushButton("Folder...", self.centralwidget)
        folder_row.addWidget(self.FilePathLE, stretch=1)
        folder_row.addWidget(self.BrowseFolderButton)
        menu_layout.addLayout(folder_row)

        pick_row = QtWidgets.QHBoxLayout()
        # Loading has its own window (ui.loader_dialog): which beamline wrote
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
            "<b>Double-click a file to open it</b> in its own window; right-click a "
            "single file for \"Show information\". Map files open on the constant-energy "
            "contour, with buttons there for the two cuts.",
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

        button_row = QtWidgets.QHBoxLayout()
        self.OpenSelectedButton = QtWidgets.QPushButton("Open selected", self.centralwidget)
        self.OpenSelectedButton.setToolTip("Same as double-clicking the highlighted file.")
        self.DataOpsButton = QtWidgets.QPushButton("Data operations...", self.centralwidget)
        self.DataOpsButton.setToolTip(
            "Truncate, self-normalise or compress the selected datasets. Several at "
            "once, as long as they have the same kind and the same axis lengths.")
        self.ProcessButton = QtWidgets.QPushButton("Process...", self.centralwidget)
        self.ProcessButton.setToolTip(
            "Smooth, differentiate, take the curvature, symmetrise, subtract a "
            "background, normalise or fit peaks. Every result is added to this "
            "list as a new dataset carrying the steps that made it.")
        button_row.addWidget(self.OpenSelectedButton, stretch=1)
        button_row.addWidget(self.DataOpsButton, stretch=1)
        button_row.addWidget(self.ProcessButton, stretch=1)
        menu_layout.addLayout(button_row)

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

        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        MainWindow.setStatusBar(self.statusbar)


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
