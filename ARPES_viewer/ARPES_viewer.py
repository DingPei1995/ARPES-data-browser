#!/usr/bin/env python3
"""
ARPES_viewer.py
============
Launcher for the nano-ARPES ``*.nxs`` viewer -- the ``*.nxs`` counterpart of
the original ``h5Loader.py``.

The main window is deliberately just a browser: pick a folder, pick the
file(s) you care about, click one to read its metadata, and **double-click**
it to open it in its own viewer window (``ui/windows.py``):

- a real-space scan opens the spatial map together with the E-vs-k spectrum
  at the cursor;
- a single cut opens that spectrum;
- a deflector map opens on its constant-energy contour, whose two buttons
  open the orthogonal cuts, each in a further window.

One window per thing you are looking at means several files -- or a contour
and both its cuts -- can be on screen at once and arranged freely, which the
single stacked layout this replaced could not do. The colormap chosen here
applies to every open window; metadata stays here because it describes the
file rather than any one view, and each viewer exports what it itself shows.

Differences from ``h5Loader.py`` worth knowing about if you maintain both:
no dependency on ``pyNanoScanSystem``/``tools_packages`` or ``win32*`` (see
``tools/system.py`` and ``ui/notify.py`` for the stand-ins), and no "Sync with running scan" button,
which talked to the live acquisition software rather than to the file format.
"""
import csv
import os
import sys
import time
import traceback

from PyQt5 import QtCore
from PyQt5.QtWidgets import (QApplication, QMainWindow, QTableWidgetItem,
                              QListWidgetItem, QFileDialog, QMenu, QMessageBox,
                              QInputDialog, QDialog)

os.environ["QT_ENABLE_HIGHDPI_SCALING"] = "1"
QApplication.setHighDpiScaleFactorRoundingPolicy(QtCore.Qt.HighDpiScaleFactorRoundingPolicy.PassThrough)

from tools import system
from ui import notify
import numpy as np

from loader.nxs_file import list_datasets, save_dataset, CUBE_KINDS
from tools.dataops import axis_summary, same_format, ARRAY_AXES
from ui.widgets import NxsData, MemoryData, DEFAULT_COLORMAP
# Imported under short names, not as `import ui.jobs`: the launcher's own
# global `ui` (the built main window, used on nearly every line below) would
# otherwise shadow the `ui` package the moment it is assigned, and
# `ui.main_window` would start looking for an attribute on a window.
from ui import jobs as nxs_jobs
from ui import main_window as main_window_ui
from ui import windows as viewer_windows_module
from loader import registry as nxs_loaders
from loader import session as nxs_session

CONFIG_PATH = os.path.join(os.path.expanduser("~"), ".arpes_viewer", "config.json")

win = None
ui = None
current_data: NxsData = None     # the file whose metadata is on show
current_info_rows = []           # full metadata, before the filter box
#: (abspath, entry) -> {path, name, entry, kind, title, start_time, memory}.
#: One record per dataset, because a file can hold several independent ones,
#: and because a converted k-map is a dataset with no file behind it.
loaded_items = {}
current_key = None               # which dataset's metadata is on show
viewer_windows = []                 # open viewer windows, kept alive here
data_ops_dialog = None              # the Data operations dialog, while open
process_dialog = None               # the Process dialog (2D or 3D), while open
info_window = None                  # the Dataset information window, once asked for
default_colormap = DEFAULT_COLORMAP   # what a newly opened window starts with
default_flip = False

#: This run's working folder for computed datasets. Every derived dataset is
#: written here the moment it exists (see :mod:`loader.session`), so that losing
#: the program stops meaning losing the session's work, and so that a listed
#: dataset can be a path rather than a resident cube.
session = nxs_session.SessionStore()

#: Kinds that are *not* auto-saved: a 4-D spatial scan is the measurement
#: itself, gigabytes of it, and writing a copy after every operation would
#: cost more than it protects. Everything smaller is written.
NO_AUTOSAVE_KINDS = ("spem_4d", "spem_1d")

#: What a kind is called in the list and in the Data Information panel.
#: The internal names are terse on purpose; these are what a person reads.
KIND_LABELS = {
    "cut": "Cut",
    "map": "Map",
    "k_map": "k-map",
    "kz_map": "kz map",
    "kz_map_k": "kz map (k)",
    "spem_4d": "SPEM",
    "spem_1d": "SPEM",
}


# --------------------------------------------------------------------------
# File selection
# --------------------------------------------------------------------------
def BrowseFolder():
    start = ui.FilePathLE.displayText() or os.path.expanduser("~")
    folder = QFileDialog.getExistingDirectory(win, "Choose folder with .nxs files", start)
    if folder:
        ui.FilePathLE.setText(folder)


def PickFiles():
    """Explicitly choose which data files to list, rather than scanning and
    loading a whole folder (slow when it holds multi-GB scans)."""
    from ui.loader_dialog import LoaderDialog

    start = ui.FilePathLE.displayText() or os.path.expanduser("~")
    if not os.path.isdir(start):
        start = os.path.expanduser("~")
    # The same filter the Loader window builds, from what the loaders say
    # they read -- so a beamline whose files are not .nxs is visible here too.
    paths, _ = QFileDialog.getOpenFileNames(win, "Select data file(s)", start,
                                            LoaderDialog._file_filter())
    if not paths:
        return
    ui.FilePathLE.setText(os.path.dirname(paths[0]))
    add_files(paths)


def open_loader():
    """The Loader window: files, which beamline's reader, and the axis
    options, all before anything reaches the list."""
    from ui.loader_dialog import LoaderDialog

    dialog = LoaderDialog(win, ui.FilePathLE.displayText())
    if dialog.exec_() != QDialog.Accepted:
        return None
    paths, options = dialog.chosen()
    if not paths:
        return None
    ui.FilePathLE.setText(os.path.dirname(paths[0]))
    add_files(paths, options=options)
    return paths


def add_files(paths, options=None, background: bool = True):
    """List the chosen files, one row per *dataset* they contain.

    A .nxs file routinely holds several complete measurements -- up to six
    Cuts or Maps recorded one after another, each its own NeXus entry. They
    are separate measurements sharing a file, so each gets its own row and
    opening one opens only that entry. Files with a single dataset still
    show a single row, labelled just by their kind.

    The listing reads only structure and dataset shapes
    (loader.nxs_file.list_datasets), so adding a folder of multi-GB scans is
    instant.
    """
    if background and len(paths) > 1 and not nxs_jobs.busy():
        def work(report):
            scanned = []
            for index, path in enumerate(paths):
                report(index / len(paths),
                       f"Reading {os.path.basename(path)} "
                       f"({index + 1} of {len(paths)})")
                scanned.append((path, nxs_loaders.list_entries(path, options)))
            return scanned

        started = nxs_jobs.run_job(
            win, "Loading data", work,
            on_done=lambda scanned: _add_scanned(scanned, options),
            on_error=lambda exc: QMessageBox.warning(
                win, "Load", f"Could not read the files:\n{exc}"),
            on_cancel=lambda: ui.statusbar.showMessage("Loading cancelled"))
        if started:
            return

    _add_scanned([(path, nxs_loaders.list_entries(path, options))
                  for path in paths], options)


def _add_scanned(scanned, options=None):
    """Put the results of scanning some files into the list.

    Split from :func:`add_files` because the scanning is file I/O that
    belongs on the worker thread and this is widget work that has to be on
    the GUI thread -- which is the whole shape of moving a long operation
    off the event loop.
    """
    added = 0
    for path, datasets in scanned:
        name = os.path.basename(path)
        if not datasets:
            # Unreadable, or a layout with nothing openable in it. Still
            # listed, so the file does not silently vanish from the browser.
            datasets = [{"entry": None, "kind": "unknown", "title": None,
                         "start_time": None}]

        multiple = len(datasets) > 1
        for info in datasets:
            # An entry may name the file it should really be read from. A
            # CASSIOPEE folder uses this: every member of a numbered series
            # offers the same assembled map, and each names the folder's
            # first member as its path. So the key is the same whether one
            # file was selected or all sixty, and the folder is listed once
            # instead of sixty times over.
            source = info.get("path") or path
            key = (os.path.abspath(source), info["entry"])
            if key in loaded_items:
                continue
            # A loader may know what the dataset is actually called -- the
            # native one does, because the name was saved with it. Only fall
            # back to the filename when nothing better is offered.
            row_name = info.get("name") or name
            loaded_items[key] = {"memory": None,
                                 "backing": None, "exported": True,
                                 "options": options, **info,
                                 "path": source, "name": row_name}

            # Only name the entry when the file holds more than one, so
            # ordinary single-measurement files stay uncluttered.
            label = f"{row_name}   [{info['kind']}]"
            if multiple and info.get("entry") and info["entry"] != row_name:
                label = f"{row_name}  \u00b7  {info['entry']}   [{info['kind']}]"
            item = QListWidgetItem(label)
            # The display text carries the labels, so the identity of the row
            # travels in its data rather than being parsed back out of it.
            item.setData(QtCore.Qt.UserRole, key)
            tip = [path]
            if info["entry"]:
                tip.append(f"entry: {info['entry']}")
            if info.get("title"):
                tip.append(f"title: {info['title']}")
            if info.get("start_time"):
                tip.append(f"recorded: {info['start_time']}")
            tip.append(info["kind"])
            item.setToolTip("\n".join(tip))
            ui.FilePathListWidget.addItem(item)
            added += 1

    ui.statusbar.showMessage(
        f"Added {added} dataset(s); {len(loaded_items)} in the list")


def listed_names():
    """Every name currently in the browser, so a new one can avoid them."""
    return [record["name"] for record in loaded_items.values()]


def unique_name(base: str) -> str:
    """``base``, or ``base 2`` / ``base 3`` ... if it is already listed."""
    taken = set(listed_names())
    if base not in taken:
        return base
    n = 2
    while f"{base} {n}" in taken:
        n += 1
    return f"{base} {n}"


def autosave(data, name: str):
    """Write a just-computed dataset to the session folder.

    Called for every derived dataset the moment it exists, which is the
    whole answer to "the program died and I lost the morning's work": there
    is no window in which a result exists only as a numpy array. Returns the
    file it wrote, or None if this kind is not auto-saved or the write
    failed.

    A failure here is reported and then ignored. Auto-saving is a safety
    net; if the disk is full or the folder is read-only, the right outcome
    is a working program that says so, not a refused k conversion.
    """
    if getattr(data, "kind", None) in NO_AUTOSAVE_KINDS:
        return None
    try:
        return session.store_data(data, name)
    except Exception as exc:                               # noqa: BLE001
        traceback.print_exc()
        ui.statusbar.showMessage(f"Could not auto-save '{name}': {exc}")
        return None


def add_memory_dataset(data, label: str):
    """List a dataset that was computed rather than read -- a k-space map,
    an arbitrary cut, a processed result.

    It gets a row like any other, keyed by its own identity rather than by a
    file and entry. It is also written to the session folder straight away,
    so the row has a file behind it exactly as a measured dataset does; the
    in-memory copy then becomes a cache that the memory budget may drop and
    re-read. Converting the same map twice gives two rows, since the
    settings (and so the result) differ.
    """
    # Show the same word a file's row shows for that kind, so a computed map
    # and a measured one read alike in the list.
    kind = KIND_LABELS.get(getattr(data, "kind", ""),
                           getattr(data, "kind", "unknown"))
    label = unique_name(label)
    key = ("<converted>", f"{label}#{len(loaded_items)}")
    backing = autosave(data, label)
    loaded_items[key] = {"path": getattr(data, "path", ""), "name": label,
                         "entry": None, "kind": kind, "title": None,
                         "start_time": None, "memory": data,
                         "backing": backing, "exported": False}
    if backing is not None:
        session.budget.put(backing, data)
    item = QListWidgetItem(f"{label}   [{kind}]")
    item.setData(QtCore.Qt.UserRole, key)
    where = (f"auto-saved to {os.path.basename(backing)}" if backing
             else "computed in this session, not auto-saved")
    item.setToolTip(f"{label}\n{where}\n{kind}")
    ui.FilePathListWidget.addItem(item)
    # Make the new row the *only* selection, so a right-click straight after
    # a conversion acts on it alone rather than on it plus whatever was
    # selected before.
    ui.FilePathListWidget.clearSelection()
    ui.FilePathListWidget.setCurrentItem(item)
    item.setSelected(True)
    ui.statusbar.showMessage(
        f"Added {label} [{kind}] to the list"
        + (" and auto-saved it" if backing else ""))
    return item


def selected_key():
    """Identity of the highlighted row: ``(abspath, entry)``, or
    ``("<converted>", tag)`` for a computed one. The row's label carries the
    kind and entry name for reading, which are not part of that identity."""
    items = ui.FilePathListWidget.selectedItems()
    if not items:
        return None
    return items[0].data(QtCore.Qt.UserRole)


def label_for(key) -> str:
    """How to name one dataset in messages: the filename, plus the entry
    when its file holds more than one."""
    record = loaded_items.get(key)
    if record is None:
        return "(unknown)"
    if record.get("memory") is not None:
        return record["name"]
    siblings = [k for k in loaded_items if k[0] == key[0]]
    if len(siblings) > 1 and record["entry"]:
        return f"{record['name']} \u00b7 {record['entry']}"
    return record["name"]


def current_label() -> str:
    """How to name the selected dataset in messages and window titles: the
    filename, plus the entry when its file holds more than one."""
    record = loaded_items.get(current_key)
    if record is None:
        return "(nothing selected)"
    if record.get("memory") is not None:
        return record["name"]
    siblings = [k for k in loaded_items if k[0] == current_key[0]]
    if len(siblings) > 1 and record["entry"]:
        return f"{record['name']} \u00b7 {record['entry']}"
    return record["name"]


def ClearList():
    loaded_items.clear()
    ui.FilePathListWidget.clear()
    ui.statusbar.showMessage("List cleared")


# --------------------------------------------------------------------------
# Right-click: rename, save, remove
# --------------------------------------------------------------------------
def selected_keys():
    """Identities of every highlighted row, in list order."""
    return [item.data(QtCore.Qt.UserRole)
            for item in ui.FilePathListWidget.selectedItems()]


def _row_for_key(key):
    for i in range(ui.FilePathListWidget.count()):
        item = ui.FilePathListWidget.item(i)
        if item.data(QtCore.Qt.UserRole) == key:
            return item
    return None


def _row_label(key) -> str:
    """The text a row should carry, from its record. Kept in one place so
    renaming and listing cannot drift apart."""
    record = loaded_items[key]
    label = f"{record['name']}   [{record['kind']}]"
    if record.get("memory") is None and record.get("entry"):
        siblings = [k for k in loaded_items if k[0] == key[0]]
        if len(siblings) > 1:
            label = (f"{record['name']}  \u00b7  {record['entry']}   "
                     f"[{record['kind']}]")
    return label


def show_list_menu(position):
    """The list's right-click menu. Acts on the whole selection, so a batch
    of datasets can be saved or removed in one go."""
    keys = selected_keys()
    menu = QMenu(ui.FilePathListWidget)

    open_action = menu.addAction("Open")
    open_action.setEnabled(len(keys) == 1)
    info_action = menu.addAction("Show information")
    info_action.setEnabled(len(keys) == 1)
    rename_action = menu.addAction("Rename...")
    rename_action.setEnabled(len(keys) == 1)
    menu.addSeparator()
    figure_action = menu.addAction(
        "Plot as a figure..." if len(keys) <= 1
        else f"Plot {len(keys)} as one figure...")
    figure_action.setEnabled(bool(keys))
    figure_action.setToolTip(
        "One panel per selected dataset, in the figure composer.")
    process_action = menu.addAction("Process...")
    process_action.setEnabled(bool(keys))
    process_action.setToolTip(
        "Smooth, differentiate, curvature, symmetrise, subtract a "
        "background, despike.")
    stack_action = menu.addAction("Stack plot (EDC / MDC)...")
    stack_action.setEnabled(len(keys) == 1)
    fit_action = menu.addAction("MDC / EDC fit...")
    fit_action.setEnabled(len(keys) == 1)
    fit_action.setToolTip(
        "Fit peaks line by line, then read the band off the fitted centres: "
        "Fermi velocity, effective mass, self-energy. Converted cuts only.")
    view3d_action = menu.addAction("3D view...")
    view3d_action.setEnabled(len(keys) == 1)
    view3d_action.setToolTip(
        "Orthogonal slices, the notched cube, an isosurface and the volume "
        "projections. Maps only.")
    # Was "Compare the two...", which already computed A - B, A / B and
    # (A - B)/(A + B) and put the result in the list -- arithmetic under a
    # name that promised a look. It is now the same window as a cut viewer's
    # "Cut arithmetic..." button, named for what it does.
    arithmetic_action = menu.addAction("Cut arithmetic on the two...")
    arithmetic_action.setEnabled(len(keys) == 2)
    arithmetic_action.setToolTip(
        "Combine two cuts: linear or circular dichroism, dividing by a "
        "reference, A \u2212 B, A / B, (A \u2212 B)/(A + B), A + B. The first "
        "selected is A; the window can swap them.")
    menu.addSeparator()
    save_action = menu.addAction(
        "Save dataset..." if len(keys) <= 1 else f"Save {len(keys)} datasets...")
    save_action.setEnabled(bool(keys))
    menu.addSeparator()
    remove_action = menu.addAction(
        "Remove from list" if len(keys) <= 1 else f"Remove {len(keys)} from list")
    remove_action.setEnabled(bool(keys))
    menu.addSeparator()
    log_action = menu.addAction("Session log...")
    log_action.setToolTip(
        "What has been auto-saved, saved and deleted this session, and "
        "where each of it is -- the file to look at after a crash.")

    chosen = menu.exec_(ui.FilePathListWidget.mapToGlobal(position))
    if chosen is None:
        return
    if chosen is open_action:
        open_selected()
    elif chosen is info_action:
        show_information()
    elif chosen is rename_action:
        rename_selected()
    elif chosen is figure_action:
        plot_selected_as_figure()
    elif chosen is process_action:
        open_processing()
    elif chosen is stack_action:
        open_stack_plot()
    elif chosen is fit_action:
        open_curve_fit()
    elif chosen is view3d_action:
        open_volume_view()
    elif chosen is arithmetic_action:
        open_cut_arithmetic()
    elif chosen is save_action:
        save_selected()
    elif chosen is log_action:
        show_operations_log()
    elif chosen is remove_action:
        remove_selected()


def plot_selected_as_figure():
    """Every selected dataset as one panel of a figure.

    The MATLAB tool's ``massplot_2D``, minus its guesswork: each panel keeps
    its own axes and labels, and anything that is not a plain 2-D dataset is
    reported rather than silently squeezed into one.
    """
    from ui.figure import FigureWindow, panel_from_dataset

    keys = selected_keys()
    if not keys:
        return None
    panels, problems = [], []
    QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
    try:
        for key in keys:
            record = loaded_items[key]
            data = record.get("memory")
            opened = data is None
            if opened:
                try:
                    data = NxsData.acquire(record["path"], record["entry"])
                except Exception as exc:
                    problems.append(f"{record['name']}: {exc}")
                    continue
            try:
                panels.append(panel_from_dataset(data, record["name"]))
            except Exception as exc:
                problems.append(f"{record['name']}: {exc}")
            finally:
                if opened:
                    data.close()
    finally:
        QApplication.restoreOverrideCursor()

    if not panels:
        QMessageBox.warning(win, "Plot as a figure",
                            "Nothing could be plotted:\n" + "\n".join(problems))
        return None
    window = figure_window(panels)
    if problems:
        ui.statusbar.showMessage(
            f"{len(panels)} panel(s) plotted, {len(problems)} skipped")
    else:
        ui.statusbar.showMessage(f"{len(panels)} panel(s) plotted")
    return window


#: figure composers opened from the launcher, kept alive
figure_windows = []


def figure_window(panels):
    from ui.figure import FigureWindow

    window = FigureWindow(panels, "Figure",
                          dataset_source={"entries": listed_datasets,
                                          "loader": load_dataset})
    window.closed.connect(lambda w: figure_windows.remove(w)
                          if w in figure_windows else None)
    figure_windows.append(window)
    window.show()
    return window


def rename_selected():
    """Rename one dataset. This is the browser's own label for it -- the
    file on disk is untouched, and a converted k-map has no file anyway."""
    keys = selected_keys()
    if len(keys) != 1:
        return
    key = keys[0]
    record = loaded_items[key]
    new_name, ok = QInputDialog.getText(
        win, "Rename dataset", "Name:", text=record["name"])
    new_name = (new_name or "").strip()
    if not ok or not new_name or new_name == record["name"]:
        return
    if new_name in listed_names():
        QMessageBox.warning(win, "Rename",
                            f"'{new_name}' is already in the list.")
        return

    record["name"] = new_name
    memory = record.get("memory")
    if memory is not None and hasattr(memory, "source_label"):
        memory.source_label = new_name      # keep the dataset's own idea in step
    item = _row_for_key(key)
    if item is not None:
        item.setText(_row_label(key))
    ui.statusbar.showMessage(f"Renamed to {new_name}")


def remove_selected():
    """Take datasets out of the browser.

    Only the list is touched: a dataset read from a file is simply no longer
    listed, and the file stays where it is. A converted k-map exists only in
    this session, so removing it does discard it -- which is why the
    confirmation says so.
    """
    keys = selected_keys()
    if not keys:
        return
    # A computed dataset used to exist only in memory, so removing its row
    # destroyed it. It now has an auto-saved copy behind it, so the warning
    # is about *that* copy going, and only for the ones that really have no
    # other home.
    computed = [loaded_items[k]["name"] for k in keys
                if loaded_items[k].get("backing")
                and not loaded_items[k].get("exported")]
    unsaved = [loaded_items[k]["name"] for k in keys
               if loaded_items[k].get("memory") is not None
               and not loaded_items[k].get("backing")]
    detail = ""
    if computed:
        detail += (f"\n\n{len(computed)} of them were computed in this "
                   f"session. Their auto-saved copies in the session folder "
                   f"go too, and they have not been saved anywhere else.")
    if unsaved:
        detail += (f"\n\n{len(unsaved)} of them exist only in memory "
                   f"(too large to auto-save); removing them discards them.")
    answer = QMessageBox.question(
        win, "Remove from list",
        f"Remove {len(keys)} dataset(s) from the list?"
        f"\n\nFiles on disk are not deleted.{detail}",
        QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
    if answer != QMessageBox.Yes:
        return

    global current_data, current_key
    for key in keys:
        record = loaded_items.pop(key, None)
        # Taking the row away takes its auto-saved copy with it; leaving the
        # file behind would fill the session folder with results nothing
        # refers to any more.
        if record is not None and record.get("backing"):
            session.discard(record["backing"])
        item = _row_for_key(key)
        if item is not None:
            ui.FilePathListWidget.takeItem(ui.FilePathListWidget.row(item))
        if key == current_key:
            if current_data is not None and record is not None \
                    and record.get("memory") is None:
                current_data.close()
            current_data = None
            current_key = None
            fill_data_info(None)
            if info_window is not None:
                info_window.TitleLabel.setText("No dataset")
                info_window.KindLabel.setText("No file selected")
                fill_info_table({}, {})
    ui.statusbar.showMessage(f"Removed {len(keys)} dataset(s) from the list")


def _dataset_for_saving(key):
    """Gather one row into the plain form save_dataset writes.

    Loads the dataset if it is not the one on show; a k-map is already in
    memory. Returns None (with a message) for anything that has no simple
    axes-and-cube form.
    """
    record = loaded_items[key]
    if record.get("memory") is not None:
        data = record["memory"]
    else:
        try:
            data = NxsData.acquire(record["path"], record["entry"])
        except Exception as exc:
            traceback.print_exc()
            return None, f"{record['name']}: could not be read ({exc})"
    try:
        # nxs_session.scan_to_dict works from loader.nxs_file.AXIS_SLOTS, the
        # one table that says which axes a kind has -- this used to be a
        # fifth hand-written copy of it, and the one that quietly wrote a
        # 4-D spatial scan with only three of its four axes.
        item = nxs_session.scan_to_dict(data.kind, data.scan, record["name"])
        # Saving means writing every point, lazy array or not.
        item["value"] = np.asarray(item["value"])
        return item, None
    except ValueError as exc:
        return None, f"{record['name']}: {exc}"
    finally:
        if record.get("memory") is None:
            data.close()


def _rebind_saved_rows(keys, path: str, entries):
    """Point the rows that were just saved at the file they were saved to,
    and delete their working copies.

    Two things happen here, and the order matters. The row's identity moves
    to ``(the file the user chose, the entry inside it)``, so re-opening it
    later reads *that* file -- and only then is the session copy removed,
    because deleting it first would leave a row briefly pointing at nothing.

    Without this the session folder accumulated a duplicate of everything
    already saved, and the copy, not the saved file, stayed the thing the
    program would re-read.
    """
    for key, entry in zip(keys, entries or []):
        record = loaded_items.get(key)
        if record is None:
            continue
        record["exported"] = True
        backing = record.get("backing")
        if not backing:
            continue                       # read from a file to begin with
        record["path"] = path
        record["entry"] = entry
        record["backing"] = None
        record["memory"] = None            # re-read from the saved file now
        session.released(backing, f"{os.path.basename(path)}::{entry}")
        item = _row_for_key(key)
        if item is not None:
            item.setToolTip(f"{record['name']}\nsaved to {path}\n"
                            f"entry: {entry}\n{record['kind']}")


def save_selected():
    """Write the selected datasets to one .nxs file this program can reopen.

    Several datasets go into one file as several entries, so saving a
    selection and opening the result gives the same rows back. Each keeps
    the name it has in the list, which is what "Rename..." is for.
    """
    keys = selected_keys()
    if not keys:
        return
    folder = ui.FilePathLE.displayText() or os.path.expanduser("~")
    default = loaded_items[keys[0]]["name"] if len(keys) == 1 else "datasets"
    default = "".join(ch for ch in default if ch.isalnum() or ch in "._- ")
    path, _ = QFileDialog.getSaveFileName(
        win, "Save dataset(s)", os.path.join(folder, f"{default}.nxs"),
        "NeXus files (*.nxs)")
    if not path:
        return
    if not path.lower().endswith(".nxs"):
        path += ".nxs"

    QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
    try:
        payload, problems, saved_keys = [], [], []
        for key in keys:
            item, problem = _dataset_for_saving(key)
            if item is not None:
                payload.append(item)
                # In step with `payload`, so the entry names save_dataset
                # reports back can be matched to the rows that produced them.
                saved_keys.append(key)
            if problem:
                problems.append(problem)
        if not payload:
            QMessageBox.warning(win, "Save",
                                "Nothing could be saved:\n" + "\n".join(problems))
            return
        entries = save_dataset(path, payload)
    except Exception as exc:
        traceback.print_exc()
        QMessageBox.warning(win, "Save", f"Could not write the file:\n{exc}")
        return
    finally:
        QApplication.restoreOverrideCursor()

    _rebind_saved_rows(saved_keys, path, entries)

    message = f"Saved {len(payload)} dataset(s) to {os.path.basename(path)}"
    if problems:
        message += f" ({len(problems)} skipped)"
        QMessageBox.information(win, "Save", message + ":\n" + "\n".join(problems))
    ui.statusbar.showMessage(message)


# --------------------------------------------------------------------------
# Metadata (single click)
# --------------------------------------------------------------------------
def load_selected(item=None):
    """Load the highlighted dataset and show its metadata. Returns the
    NxsData (or KMapData), or None if it could not be read.

    Metadata is per entry, not per file: five Maps in one file were recorded
    at different times with different settings, so selecting a row reads
    that entry's own metadata.
    """
    global current_data, current_key
    key = selected_key()
    if key is None:
        return None
    record = loaded_items.get(key)
    if record is None:
        return None
    if (record.get("memory") is None and not record.get("backing")
            and not os.path.isfile(record["path"])):
        ui.statusbar.showMessage(f"File no longer available: {key[0]}")
        return None

    current_key = key
    if current_data is not None:
        current_data.close()      # release the previous dataset's claim
        current_data = None
    try:
        current_data = load_dataset(key)
    except Exception as exc:
        traceback.print_exc()
        current_data = None
        ensure_info_window().KindLabel.setText("Load error")
        fill_data_info(None)
        fill_info_table({}, {})
        ui.statusbar.showMessage(f"Failed to load {current_label()}: {exc}")
        return None

    info = ensure_info_window()
    info.TitleLabel.setText(current_label())
    info.KindLabel.setText(f"kind: {current_data.kind}")
    fill_data_info(current_data)
    fill_info_table(current_data.scan.info, current_data.scan.fourd_info)
    ui.statusbar.showMessage(f"{current_label()}: {current_data.kind} "
                              f"-- double-click to open it")
    return current_data


#: Which axes each kind puts in the X / Y / Z / E columns of the Data
#: Information panel, as attribute names on the scan, and what to call them.
#: A cut is (angle, energy); a map is its two in-plane axes (angles, or
#: momenta once converted) plus energy; a spatial scan is the two sample
#: coordinates, then the analyser angle and the energy -- which is why it is
#: the one kind that needs the fourth column.
DATA_INFO_AXES = {
    "cut": (("x", "angle"), ("y", "energy")),
    "map": (("x", "angle/k"), ("k", "angle/k"), ("z", "energy")),
    "k_map": (("x", "k"), ("k", "k"), ("z", "energy")),
    "kz_map": (("x", "photon energy"), ("k", "angle/k"), ("z", "energy")),
    "kz_map_k": (("x", "kz"), ("k", "k"), ("z", "energy")),
    "spem_4d": (("x", "spatial"), ("y", "spatial"), ("k", "angle"), ("z", "energy")),
    "spem_1d": (("x", "spatial"), ("k", "angle"), ("z", "energy")),
}


def ensure_info_window():
    """The information window, made the first time it is asked for.

    It is not part of the launcher's layout: information is read on
    request, and a window of its own can sit beside a viewer, stay open
    while other files are browsed, or be closed and forgotten.
    """
    global info_window
    if info_window is None:
        info_window = main_window_ui.InfoWindow(win)
        info_window.InfoFilterLE.textChanged.connect(apply_info_filter)
        info_window.ExportMetaButton.clicked.connect(ExportMetadataCsv)
    return info_window


def show_information():
    """Right-click -> Show information: read the highlighted dataset and put
    what it is in front of the user."""
    window = ensure_info_window()
    if load_selected() is None:
        return None
    window.show()
    window.raise_()
    window.activateWindow()
    return window


def fill_data_info(data):
    """Describe the dataset's extent in the Data Information panel: the
    first and last value, the number of points and the step, per axis.

    The column headers carry the axis's own label, so "Z" never has to be
    guessed at -- on a map it is the energy, on a spatial scan the analyser
    angle, and the panel says so.
    """
    for cell in ui.DataInfoCells.values():
        cell.setText("")
    ui.DataInfoBox.setTitle("Data Information" if data is None
                            else f"Data Information \u2014 {current_label()}")
    columns = DATA_INFO_AXES.get(getattr(data, "kind", ""), ()) if data is not None else ()
    labels = getattr(getattr(data, "scan", None), "labels", {}) or {}
    for index, header in enumerate(ui.DataInfoHeaders):
        name = ("X", "Y", "Z", "E")[index]
        if index < len(columns):
            attribute, role = columns[index]
            label = labels.get(attribute) or role
            # A spatial axis is already called "X (mm)"; repeating the column
            # letter above it would just say X twice.
            header.setText(label if label.startswith(name + " ")
                           else f"{name}\n{label}")
            header.setVisible(True)
        else:
            header.setText(name)
            header.setVisible(index < 3)      # X/Y/Z always, E only for 4D
        for row in ("min", "max", "num", "step"):
            ui.DataInfoCells[(row, index)].setVisible(
                index < max(3, len(columns)))
    if data is None:
        return
    for index, (attribute, _role) in enumerate(columns):
        summary = axis_summary(getattr(data.scan, attribute, None))
        if summary is None:
            continue
        lo, hi, num, step = summary
        for row, value in (("min", f"{lo:.6g}"), ("max", f"{hi:.6g}"),
                           ("num", str(num)), ("step", f"{step:.6g}")):
            ui.DataInfoCells[(row, index)].setText(value)


def fill_info_table(info: dict, fourd_info: dict):
    global current_info_rows
    current_info_rows = ([(k, v) for k, v in info.items()]
                         + [(f"Motor.{k}", v) for k, v in fourd_info.items()])
    apply_info_filter(ensure_info_window().InfoFilterLE.displayText())


def apply_info_filter(pattern: str):
    pattern = (pattern or "").strip().lower()
    rows = [(k, v) for k, v in current_info_rows
            if not pattern or pattern in str(k).lower() or pattern in str(v).lower()]
    table = ensure_info_window().FileInfoTW
    table.setRowCount(len(rows))
    for i, (key, value) in enumerate(rows):
        table.setItem(i, 0, QTableWidgetItem(str(key)))
        table.setItem(i, 1, QTableWidgetItem(str(value)))
    table.resizeColumnsToContents()


# --------------------------------------------------------------------------
# Opening viewers (double click)
# --------------------------------------------------------------------------
def open_selected(item=None):
    """Open the highlighted dataset in its own viewer window -- that entry
    only, even when its file holds five more beside it."""
    key = selected_key()
    if key is None:
        return None
    if current_data is None or key != current_key:
        if load_selected() is None:
            return None
    record = loaded_items[key]

    # A reference of its own, so this window's lifetime is independent of
    # the browser selection and of any other window on the same file. A
    # dataset computed here is handed over as it stands, or re-read from the
    # copy it was auto-saved to if the budget has since dropped it.
    window_data = load_dataset(key)
    if window_data is None:
        return None

    window = viewer_windows_module.open_viewer(
        window_data, current_label(), default_colormap, default_flip)
    if window is None:
        window_data.close()
        ui.statusbar.showMessage(
            f"'{current_label()}' has no view: unsupported layout "
            f"(case {current_data.scan.info.get('_case')}, "
            f"group {current_data.scan.info.get('_group')})")
        return None

    _adopt_viewer(window, key)
    window.closed.connect(_forget_viewer)
    viewer_windows.append(window)
    window.show()
    ui.statusbar.showMessage(
        f"Opened {current_label()} ({len(viewer_windows)} window(s) open)")
    return window


def _adopt_viewer(window, key=None):
    """Hook a new viewer up to the list it came from.

    A viewer can compute datasets of its own -- a k-space map, an
    arbitrary-direction cut, a saved slice, a Fermi-corrected map. They land
    in this list under whatever name their dialog was given, the dialog asks
    which names are taken so its default does not collide, and a viewer that
    wants to show its result straight away opens it through here, so the new
    window is registered exactly like a double-clicked one.
    """
    window._name_source = listed_names
    window._dataset_opener = open_dataset
    # So a viewer can work from *another* listed dataset -- the Fermi-surface
    # correction fits a gold reference, which is a different measurement.
    window._dataset_entries = listed_datasets
    window._dataset_loader = load_dataset
    # So a cut viewer can say "now click the other cut in the list" -- and
    # recognise its own row when it is clicked. The row, not the object: a
    # cut read from a file is re-read on every load, so the same row hands
    # back a different object each time.
    window._await_selection = await_list_selection
    window._list_key = key
    if hasattr(window, "kMapCreated"):
        window.kMapCreated.connect(
            lambda data: add_memory_dataset(data, getattr(data, "source_label", "k map")))
    if hasattr(window, "datasetCreated"):
        window.datasetCreated.connect(
            lambda data: add_memory_dataset(data,
                                            getattr(data, "source_label", "computed")))


def listed_datasets():
    """(label, key) for every row, for a viewer that needs to work from
    another dataset."""
    return [(f"{label_for(key)}   [{loaded_items[key]['kind']}]", key)
            for key in loaded_items]


def load_dataset(key):
    """The dataset behind one row, loaded if it is not already in memory.

    Three cases, in this order:

    * still resident -- hand it back;
    * computed earlier and auto-saved -- re-read it from the session folder,
      which is what lets the memory budget drop a cube it is not using and
      what makes a row cost a path rather than a cube;
    * read from a file -- open the file, as always.
    """
    record = loaded_items[key]
    if record.get("memory") is not None:
        return record["memory"]

    backing = record.get("backing")
    if backing:
        cached = session.budget.get(backing)
        if cached is not None:
            return cached
        data = _reload_backing(record, backing)
        if data is not None:
            return data

    options = record.get("options")
    return _acquire(record["path"], record["entry"], options)


def _acquire(path, entry, options):
    """Open a file, on the worker thread when the reader says it will be slow.

    Most readers open a file with one call into HDF5 and there is nothing to
    report or to wait for. A few genuinely take a while -- assembling a
    CASSIOPEE folder means parsing a hundred text files -- and those are the
    ones that declare a ``progress`` argument. Those go through
    :func:`ui.jobs.run_blocking`, so the window keeps painting and there is a
    bar and a Cancel instead of a frozen panel.

    The distinction is the reader's own: a loader that can report progress is
    one that expected to need to. Everything else keeps the direct path it
    has always had, which is both faster and a smaller thing to get wrong.
    """
    loader = None
    try:
        if options is not None and getattr(options, "loader", None):
            loader = nxs_loaders.get_loader(options.loader)
        else:
            loader = nxs_loaders.detect(path)
    except Exception:                                      # noqa: BLE001
        loader = None

    slow = loader is not None and nxs_loaders._accepts_progress(loader.load)
    if not slow or nxs_jobs.busy():
        return NxsData.acquire(path, entry, options=options)

    title = f"Reading {os.path.basename(path)}"

    def work(report):
        def progress(done, total, label):
            report(done / max(1, total), f"{label}  ({done + 1} of {total})")
        return NxsData.acquire(path, entry, options=options, progress=progress)

    try:
        return nxs_jobs.run_blocking(win, title, work)
    except nxs_jobs.JobCancelled:
        # Raised rather than returned as None: every caller of load_dataset
        # goes straight on to use the dataset, and all but one of them guard
        # with try/except rather than a None check. An exception is the
        # contract they already handle; a None would be an AttributeError
        # one line later.
        ui.statusbar.showMessage("Loading cancelled")
        raise RuntimeError("loading was cancelled") from None


def _reload_backing(record, backing):
    """Re-read a derived dataset from the file it was auto-saved to.

    Returns a :class:`NxsData` over the session file rather than the
    original :class:`MemoryData`: the two present the same interface to
    every viewer, and this one keeps the array on disk until it is used.
    """
    if not os.path.isfile(backing):
        return None
    try:
        data = NxsData.acquire(backing, None)
    except Exception as exc:                               # noqa: BLE001
        traceback.print_exc()
        ui.statusbar.showMessage(
            f"Could not re-read '{record['name']}' from its auto-saved copy: {exc}")
        return None
    session.budget.put(backing, data)
    return data


def open_dataset(data):
    """Open a dataset object that is already in the list (a viewer has just
    computed it) in its own viewer window."""
    key = None
    for candidate, record in loaded_items.items():
        if record.get("memory") is data:
            key = candidate
            break
    if key is None:
        add_memory_dataset(data, getattr(data, "source_label", "computed"))
        key = selected_key()
    label = loaded_items[key]["name"] if key in loaded_items else "computed"
    window = viewer_windows_module.open_viewer(data, label, default_colormap, default_flip)
    if window is None:
        return None
    _adopt_viewer(window, key)
    window.closed.connect(_forget_viewer)
    viewer_windows.append(window)
    window.show()
    return window


def open_data_operations():
    """Truncate / self-normalise / compress the selected datasets.

    All of them at once, which is the point -- but only if they really are
    the same shape, since one set of bounds or factors cannot describe two
    different axes. Datasets that disagree are refused here with a list of
    what they are, rather than producing a batch of quietly wrong results.
    """
    keys = selected_keys()
    if not keys:
        ui.statusbar.showMessage("Select one or more datasets first")
        QMessageBox.information(win, "Data operations",
                                 "Select the dataset(s) to work on first.")
        return None

    datasets, problems = [], []
    for key in keys:
        record = loaded_items[key]
        try:
            data = record.get("memory")
            if data is None:
                data = NxsData.acquire(record["path"], record["entry"])
        except Exception as exc:
            traceback.print_exc()
            problems.append(f"{label_for(key)}: could not be read ({exc})")
            continue
        if data.kind not in ARRAY_AXES:
            problems.append(f"{label_for(key)}: {data.kind} has no axes to work on")
            continue
        datasets.append((label_for(key), data))

    if problems or not datasets:
        QMessageBox.warning(win, "Data operations",
                             "These cannot be operated on:\n" + "\n".join(problems or
                             ["nothing usable was selected"]))
        return None

    shapes = []
    for name, data in datasets:
        axes = ARRAY_AXES[data.kind]
        shapes.append((data.kind, tuple(np.asarray(getattr(data.scan, a)).size
                                        for a in axes)))
    if not same_format(shapes):
        listing = "\n".join(f"  {name}: {kind} {' x '.join(str(n) for n in shape)}"
                             for (name, _), (kind, shape) in zip(datasets, shapes))
        QMessageBox.warning(
            win, "Data operations",
            "The selected datasets are not the same format, so one set of "
            "bounds or factors cannot describe them all. Select datasets with "
            "the same kind and the same number of points on every axis:\n\n"
            + listing)
        return None

    dialog = viewer_windows_module.DataOperationsDialog(datasets, win)
    dialog.datasetsCreated.connect(_list_computed)
    dialog.show()
    dialog.raise_()
    global data_ops_dialog
    data_ops_dialog = dialog
    return dialog


def open_processing():
    """Smooth / differentiate / curvature / symmetrise / fit the selection.

    Two panels rather than one, and which opens is decided by what is
    selected: a 2-D dataset goes to the plane panel, a 3-D one to the volume
    panel. A cube's first question is always which plane an operation is
    defined on -- the curvature of a constant-energy map and the curvature
    of a dispersion are different quantities that share a formula -- and
    that question has no place on a 2-D panel.
    """
    from ui.process import ProcessDialog
    from ui.volume import VolumeProcessDialog

    keys = selected_keys()
    if not keys:
        ui.statusbar.showMessage("Select one or more datasets first")
        QMessageBox.information(win, "Data processing",
                                "Select the dataset(s) to work on first.")
        return None

    datasets, problems = [], []
    for key in keys:
        try:
            data = load_dataset(key)
        except Exception as exc:
            traceback.print_exc()
            problems.append(f"{label_for(key)}: could not be read ({exc})")
            continue
        if data.kind not in ("cut",) + CUBE_KINDS:
            problems.append(f"{label_for(key)}: {data.kind} is not a cut or a map")
            continue
        datasets.append((label_for(key), data))

    if not datasets:
        QMessageBox.warning(
            win, "Data processing",
            "These cannot be processed:\n" + "\n".join(problems or
            ["nothing usable was selected"]) +
            "\n\nOpen a spatial scan first and process a cut taken from it.")
        return None

    kinds = {data.kind for _, data in datasets}
    if len(kinds) > 1:
        QMessageBox.warning(
            win, "Data processing",
            "Cuts and maps need different panels, so they cannot be processed "
            "in one go. Select one kind at a time.")
        return None

    if kinds == {"cut"}:
        dialog = ProcessDialog(datasets, win, default_colormap, default_flip)
    else:
        dialog = VolumeProcessDialog(datasets, win, default_colormap, default_flip)
    dialog.datasetsCreated.connect(_list_computed)
    dialog.show()
    dialog.raise_()
    global process_dialog
    process_dialog = dialog
    if problems:
        ui.statusbar.showMessage("; ".join(problems))
    return dialog


def open_stack_plot():
    """A waterfall of EDCs or MDCs across the selected dataset."""
    from ui.process import StackWindow

    key = selected_key()
    if key is None:
        return None
    try:
        data = load_dataset(key)
    except Exception as exc:
        QMessageBox.warning(win, "Stack plot", str(exc))
        return None
    if data.kind != "cut":
        QMessageBox.information(
            win, "Stack plot",
            "A stack is made of curves across a 2-D dataset. Open a map and "
            "take a slice of it first, or convert it to a cut.")
        return None
    scan = data.scan
    window = StackWindow(scan.value, (scan.x, scan.y),
                         (scan.labels.get("x", "x"), scan.labels.get("y", "y")),
                         label_for(key), win, default_colormap)
    window.datasetsCreated.connect(_list_computed)
    window.show()
    viewer_windows.append(window)
    return window


def open_curve_fit():
    """Fit the MDCs or EDCs of the selected cut and read the band off them.

    Only a cut whose in-plane axis is already in A^-1 gets in. Nothing in
    the fit itself would object to degrees -- the peaks would come out
    perfectly well -- but every quantity built on the fitted band would be
    wrong: v_F would be in eV.degree, m* would be meaningless, and neither
    number carries its units around to say so. Converting first is one
    button, and it is the honest order of operations.
    """
    from ui.fit import FitPanel, momentum_cut_reason

    key = selected_key()
    if key is None:
        return None
    try:
        data = load_dataset(key)
    except Exception as exc:
        QMessageBox.warning(win, "MDC / EDC fit", str(exc))
        return None
    reason = momentum_cut_reason(data)
    if reason:
        QMessageBox.information(win, "MDC / EDC fit", reason + ".")
        return None

    window = FitPanel(data, label_for(key), win, default_colormap, default_flip)
    window.datasetsCreated.connect(_list_computed)
    window.closed.connect(_forget_viewer)
    window.show()
    viewer_windows.append(window)
    return window


def open_volume_view():
    """The 3-D views of the selected map, over a chosen region."""
    from ui.volume import VolumeWindow, VolumeRangeDialog

    key = selected_key()
    if key is None:
        return None
    try:
        data = load_dataset(key)
    except Exception as exc:
        QMessageBox.warning(win, "3D view", str(exc))
        return None
    if data.kind not in CUBE_KINDS:
        QMessageBox.information(
            win, "3D view",
            "The 3-D views need a cube: a Map or a converted k-map.")
        return None
    x, k, z, cube = data.angle_cube
    scan = data.scan
    labels = (scan.labels.get("x", "x"), scan.labels.get("k", "y"),
              scan.labels.get("z", "z"))

    # Ask what to show before anything is resampled. The 3-D views rebuild
    # from the whole array on every frame, so the window and the sampling
    # are what decide whether they keep up with the mouse.
    chooser = VolumeRangeDialog(np.asarray(cube, dtype=float), (x, k, z),
                                labels, win)
    if chooser.exec_() != chooser.Accepted:
        return None
    QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
    try:
        values, axes = chooser.result()
    except Exception as exc:
        QApplication.restoreOverrideCursor()
        QMessageBox.warning(win, "3D view", str(exc))
        return None
    finally:
        QApplication.restoreOverrideCursor()

    window = VolumeWindow(values, axes, labels,
                          label_for(key), win, default_colormap, default_flip)
    window.closed.connect(_forget_viewer)
    window.show()
    viewer_windows.append(window)
    return window


def open_cut_arithmetic():
    """Two cuts from the list, combined: dichroism, ratios, differences.

    The same window a cut viewer opens from its own "Cut arithmetic..."
    button (see :mod:`ui.cutops`); this is the route for when both cuts are
    already picked in the list.
    """
    from ui.cutops import CutArithmeticDialog

    keys = selected_keys()[:2]
    pair = []
    for key in keys:
        if str(loaded_items.get(key, {}).get("kind", "")).lower() != "cut":
            QMessageBox.information(
                win, "Cut arithmetic",
                f"{label_for(key)} is not a cut. Cut arithmetic works on two "
                f"cuts; take a slice of a map first.")
            return None
        try:
            data = load_dataset(key)
        except Exception as exc:
            QMessageBox.warning(win, "Cut arithmetic", str(exc))
            return None
        pair.append((label_for(key), data))
    try:
        dialog = CutArithmeticDialog(pair[0], pair[1], win, default_colormap,
                                     default_flip, existing_names=listed_names)
    except ValueError as exc:
        QMessageBox.information(win, "Cut arithmetic", str(exc))
        return None
    dialog.datasetsCreated.connect(_list_computed)
    dialog.show()
    global process_dialog
    process_dialog = dialog
    return dialog


def await_list_selection(on_chosen, *, kind: str = "cut", exclude=None,
                         exclude_key=None, on_rejected=None):
    """Call ``on_chosen(label, data)`` with the next dataset clicked in the
    list that is a ``kind``.

    How a viewer asks for "another" dataset without opening a picker of its
    own: the list is where the datasets are, so that is where the choice is
    made. Clicks on the wrong kind, or on ``exclude`` itself, are reported
    through ``on_rejected(message)`` and the wait goes on. ``exclude`` is
    matched by identity and ``exclude_key`` by list row; a viewer passes
    both, since a cut read from a file is a new object on every load. The kind is
    checked from the list's own record *before* anything is loaded, so a
    stray click on a sixty-spectrum map costs nothing.

    Returns a function that cancels the wait.
    """
    widget = ui.FilePathListWidget

    def changed():
        key = selected_key()
        record = loaded_items.get(key) if key is not None else None
        if record is None:
            return
        label = label_for(key)
        if exclude_key is not None and key == exclude_key:
            if on_rejected:
                on_rejected("That is the cut this was opened from. Pick "
                            "another one.")
            return
        if str(record.get("kind", "")).lower() != kind:
            if on_rejected:
                on_rejected(f"\u201c{label}\u201d is a {record.get('kind')}, "
                            f"not a {kind}. Pick a {kind}.")
            return
        try:
            data = load_dataset(key)
        except Exception as exc:                            # noqa: BLE001
            if on_rejected:
                on_rejected(f"\u201c{label}\u201d could not be read: {exc}")
            return
        if exclude is not None and data is exclude:
            if on_rejected:
                on_rejected("That is the cut this was opened from. Pick "
                            "another one.")
            return
        cancel()
        on_chosen(label, data)

    def cancel():
        try:
            widget.itemSelectionChanged.disconnect(changed)
        except TypeError:
            pass                    # already disconnected

    widget.itemSelectionChanged.connect(changed)
    return cancel


def _list_computed(datasets):
    for data in datasets:
        add_memory_dataset(data, getattr(data, "source_label", "computed"))
    ui.statusbar.showMessage(f"Added {len(datasets)} computed dataset(s) to the list")


def _forget_viewer(window):
    if window in viewer_windows:
        viewer_windows.remove(window)


# --------------------------------------------------------------------------
# Metadata export (the viewers export their own data)
# --------------------------------------------------------------------------
def ExportMetadataCsv():
    if not current_info_rows:
        notify.warning("No metadata loaded yet.")
        return
    folder = ui.FilePathLE.displayText() or os.path.expanduser("~")
    # Name the CSV after the dataset, so five maps from one file do not all
    # want to be called the same thing.
    record = loaded_items.get(current_key)
    stem = os.path.splitext(record["name"])[0] if record else "metadata"
    if record and record.get("entry") and len([k for k in loaded_items
                                               if k[0] == current_key[0]]) > 1:
        stem = f"{stem}_{record['entry']}"
    default = os.path.join(folder, stem + "_metadata.csv")
    path, _ = QFileDialog.getSaveFileName(win, "Export metadata as CSV", default,
                                           "CSV (*.csv)")
    if not path:
        return
    if not path.lower().endswith(".csv"):
        path += ".csv"
    try:
        with open(path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["field", "value"])
            for key, value in current_info_rows:
                writer.writerow([key, value])
    except (OSError, PermissionError) as exc:
        notify.warning(f"Could not write the metadata:\n{exc}")
        return
    ui.statusbar.showMessage(f"Exported {len(current_info_rows)} metadata rows to {path}")


# --------------------------------------------------------------------------
# Wiring
# --------------------------------------------------------------------------
def offer_recovery():
    """If an earlier session left datasets behind, offer to open them.

    A session that ended properly leaves nothing: saving a dataset removes
    its working copy and so does removing its row. So anything still here is
    work from a run that stopped without being tidied up -- which is exactly
    the case this whole mechanism exists for, and the moment to say so is
    when the program next starts, not when the user eventually goes digging
    through a folder named after a timestamp and a process id.
    """
    leftovers = session.leftovers()
    if not leftovers:
        return 0

    total_files = sum(len(entry["files"]) for entry in leftovers)
    total_bytes = sum(entry["bytes"] for entry in leftovers)
    newest = leftovers[0]
    when = time.strftime("%Y-%m-%d %H:%M",
                         time.localtime(newest["modified"]))
    answer = QMessageBox.question(
        win, "Recover unsaved work",
        f"{total_files} dataset(s) from {len(leftovers)} earlier session(s) "
        f"were never saved to a file of your own.\n\n"
        f"Most recent: {when} ({newest['bytes'] / 1e6:.0f} MB)\n"
        f"Total: {total_bytes / 1e6:.0f} MB\n\n"
        f"Add them to the list?\n\n"
        f"They are in {session.root}\nand are deleted after "
        f"{nxs_session.KEEP_DAYS} days. Every one of them is in the "
        f"operations log beside them.",
        QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
    if answer != QMessageBox.Yes:
        ui.statusbar.showMessage(
            f"{total_files} unsaved dataset(s) left in {session.root}")
        return 0

    paths = [os.path.join(entry["folder"], filename)
             for entry in leftovers for filename in entry["files"]]
    add_files(paths, background=False)
    # Recovered rows are somebody's unsaved work, not a file they chose, so
    # they count as unsaved until they are saved -- the close prompt has to
    # keep warning about them.
    for key, record in loaded_items.items():
        if os.path.dirname(record["path"]).startswith(session.root):
            record["exported"] = False
            record["backing"] = record["path"]
    session.log("RECOVER", f"{len(paths)} dataset(s)", session.root,
                "listed at startup")
    ui.statusbar.showMessage(f"Recovered {len(paths)} unsaved dataset(s)")
    return len(paths)


def show_operations_log():
    """The operations log, in a window. What was written, when, and where."""
    rows = nxs_session.read_log(session.root)
    if not rows:
        QMessageBox.information(
            win, "Operations log",
            f"Nothing logged yet.\n\nThe log is written to\n{session.log_path}")
        return None
    lines = [" | ".join(row[:6]) for row in rows]
    box = QMessageBox(win)
    box.setWindowTitle("Operations log")
    box.setText(f"{len(rows)} most recent operations")
    box.setInformativeText(f"The full log is {session.log_path}")
    box.setDetailedText("\n".join(lines))
    box.exec_()
    return rows


def unsaved_datasets():
    """Datasets that exist only because this session computed them and that
    the user has not put anywhere of their own.

    The auto-saved copy in the session folder is not "somewhere of their
    own": it is a working copy that gets pruned after a week, so closing
    without being asked would still lose the work in the sense that matters.
    """
    return [record["name"] for record in loaded_items.values()
            if (record.get("backing") or record.get("memory") is not None)
            and not record.get("exported")]


def confirm_close() -> bool:
    """Ask before closing if there is unsaved work. True means go ahead.

    The program used to close on the spot, which was fine while a dataset
    was either a file on disk or nothing -- and painful once a session's
    worth of processing lived only in the list.
    """
    pending = unsaved_datasets()
    if not pending:
        return True

    shown = "\n".join(f"  • {name}" for name in pending[:8])
    if len(pending) > 8:
        shown += f"\n  • ... and {len(pending) - 8} more"
    box = QMessageBox(win)
    box.setWindowTitle("Close ARPES viewer")
    box.setIcon(QMessageBox.Question)
    box.setText(f"{len(pending)} dataset(s) computed in this session have "
                f"not been saved to a file of your own:")
    box.setInformativeText(
        f"{shown}\n\nThey have been auto-saved to\n{session.folder}\n"
        f"which is cleared out after {nxs_session.KEEP_DAYS} days.")
    save = box.addButton("Save to a file...", QMessageBox.AcceptRole)
    discard = box.addButton("Close anyway", QMessageBox.DestructiveRole)
    box.addButton("Cancel", QMessageBox.RejectRole)
    box.setDefaultButton(save)
    box.exec_()

    clicked = box.clickedButton()
    if clicked is discard:
        return True
    if clicked is save:
        # Select everything unsaved, then run the ordinary Save, so there is
        # one save path rather than two that can disagree.
        ui.FilePathListWidget.clearSelection()
        for key, record in loaded_items.items():
            if record["name"] in pending:
                item = _row_for_key(key)
                if item is not None:
                    item.setSelected(True)
        save_selected()
        return not unsaved_datasets()
    return False


class _MainWindow(QMainWindow):
    """The launcher window, with a close that asks first."""

    def closeEvent(self, event):
        if nxs_jobs.busy():
            answer = QMessageBox.question(
                self, "Close ARPES viewer",
                "Something is still running. Close anyway?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if answer != QMessageBox.Yes:
                event.ignore()
                return
        if confirm_close():
            event.accept()
        else:
            event.ignore()


def SetConnect():
    ui.BrowseFolderButton.clicked.connect(BrowseFolder)
    ui.LoadDataButton.clicked.connect(open_loader)
    ui.PickFilesButton.clicked.connect(PickFiles)
    ui.ClearListButton.clicked.connect(ClearList)

    # A click only selects. Reading a dataset to describe it is real work --
    # for a GB-sized spatial scan it is seconds -- so the information panel
    # is filled on request (right-click -> Show information), not on every
    # click. Double-click still opens a viewer.
    ui.FilePathListWidget.itemDoubleClicked.connect(open_selected)
    ui.FilePathListWidget.customContextMenuRequested.connect(show_list_menu)
    ui.OpenSelectedButton.clicked.connect(open_selected)

    # The metadata filter and its CSV export live in the information
    # window now, and are wired when that window is first made.
    ui.DataOpsButton.clicked.connect(open_data_operations)
    ui.ProcessButton.clicked.connect(open_processing)


def InitializeUIWidgets():
    values = system.load_config(CONFIG_PATH, defaults={
        "FilePathLE": os.path.expanduser("~"),
        "ColorMap": DEFAULT_COLORMAP,
        "FlipColorMap": False,
    })
    global default_colormap, default_flip
    ui.FilePathLE.setText(values["FilePathLE"])
    # Each viewer window owns its colormap; these are only the values a newly
    # opened window starts from.
    default_colormap = values.get("ColorMap", DEFAULT_COLORMAP)
    default_flip = bool(values.get("FlipColorMap", False))
    return values


def DumpUIWidgets(values):
    values["FilePathLE"] = ui.FilePathLE.displayText()
    # Remember the most recently opened window's choice as the next default.
    if viewer_windows:
        values["ColorMap"] = viewer_windows[-1].colormap
        values["FlipColorMap"] = viewer_windows[-1].flip
    else:
        values["ColorMap"] = default_colormap
        values["FlipColorMap"] = default_flip
    system.dump_config(CONFIG_PATH, values)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = _MainWindow()
    ui = main_window_ui.Ui_MainWindow()
    ui.setupUi(win)

    config = InitializeUIWidgets()
    SetConnect()

    session.log("START", "", session.folder, "session opened")

    # Order matters: prune first so the recovery offer is only ever about
    # folders recent enough to still be there, then offer what is left.
    try:
        gone = session.prune_old()
        if gone:
            ui.statusbar.showMessage(f"Cleared {gone} old session folder(s)")
    except Exception:                                      # noqa: BLE001
        traceback.print_exc()

    win.show()
    try:
        offer_recovery()
    except Exception:                                      # noqa: BLE001
        traceback.print_exc()
    error = app.exec_()
    session.log("END", "", session.folder, "session closed")
    DumpUIWidgets(config)
    sys.exit(error)
