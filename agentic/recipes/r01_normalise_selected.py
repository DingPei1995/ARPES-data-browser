# Recipe: normalise every row whose list text starts with PREFIX by its total intensity,
# using the viewer's own "Self-normalise" dialog (result rows get the suffix _s_nor).
# Drop into agentic/viewer_cmd/ (live viewer). Adapt PREFIX / the axes ticked, test on one row first.
from PyQt5 import QtWidgets
PREFIX = "Cut_"                      # rows to normalise, matched against the start of the list text
L = A.ui.FilePathListWidget

def clear():                         # close stray dialogs/windows before the next step
    for w in list(QtWidgets.QApplication.topLevelWidgets()):
        if w.isVisible() and w is not A.win: w.close()
    QtWidgets.QApplication.processEvents()

clear(); L.clearSelection(); n = 0
for i in range(L.count()):
    if L.item(i).text().startswith(PREFIX):
        L.item(i).setSelected(True); n += 1
dlg = A.open_data_operations("normalise")
dlg.norm_preset.setCurrentIndex(dlg.norm_preset.count() - 1)      # custom
for c in dlg.norm_checks: c.setChecked(True)                       # sum over both axes = whole cut
QtWidgets.QApplication.processEvents()
print("selected", n, "| axes:", [c.text() for c in dlg.norm_checks], "|", dlg.norm_explain.text()[:140], flush=True)
res = dlg.apply(); print("->", [r.source_label for r in res][:5], "...", flush=True)
dlg.grab().save(HERE + r"\shots\normalise_dialog.png")
dlg.close(); clear()
