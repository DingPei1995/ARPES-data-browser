# Recipe: for every original cut in the list: k conversion (cut) -> Fermi-edge fit -> offset energy axis.
# Uses the viewer's own dialogs. Parameters below are examples from one beamtime: ASK the user and
# change them (slit Gamma, deflector angles, Fermi window in kinetic energy, temperature in the name).
import re, traceback
from PyQt5 import QtWidgets
L = A.ui.FilePathListWidget
GAMMA_SLIT, CUT_DEFLECTOR, GAMMA_DEFLECTOR = 0.8, 23.5, 15.0     # deg
EF_WINDOW = (74.05, 74.45)                                        # eV kinetic, around the Fermi edge
TAG = " [k G15]"                                                  # suffix of the converted row
TEMP_FROM_NAME = r"_(\d+)K"                                       # temperature held in the Fermi fit

def clear():
    for w in list(QtWidgets.QApplication.topLevelWidgets()):
        if w.isVisible() and w is not A.win: w.close()
    QtWidgets.QApplication.processEvents()
def select(prefix):
    for i in range(L.count()):
        if L.item(i).text().startswith(prefix + "   ["):
            L.clearSelection(); L.setCurrentItem(L.item(i)); L.item(i).setSelected(True); return True
    return False

clear()
names = sorted(r["name"] for r in A.loaded_items.values()
               if r["name"].startswith("Cut_") and r.get("memory") is None and "[" not in r["name"])
print(len(names), "cuts", flush=True)
for n in names:
    try:
        T = float(re.search(TEMP_FROM_NAME, n).group(1))
        select(n); d = A.open_k_conversion("cut"); QtWidgets.QApplication.processEvents()
        # Gamma deflector angle is only used when "Advanced" is ticked -- always tick it and verify
        d.advanced.setChecked(True)
        d.gamma_slit.setValue(GAMMA_SLIT); d.deflector_box.setValue(CUT_DEFLECTOR); d.gamma_defl.setValue(GAMMA_DEFLECTOR)
        d.name_box.setText(n + TAG); QtWidgets.QApplication.processEvents()
        if n == names[0]:
            print("dialog:", d.gamma_slit.value(), d.deflector_box.value(), d.gamma_defl.value(), d.advanced.isChecked(), flush=True)
            d.grab().save(HERE + r"\shots\k_conversion_cut_dialog.png")
        [b for b in d.findChildren(QtWidgets.QPushButton) if b.text() == "Convert"][0].click()
        QtWidgets.QApplication.processEvents(); clear()
        select(n + TAG); fd = A.open_fermi_fit()
        fd.value_boxes["temperature"].setValue(T); fd.energy_from.setValue(EF_WINDOW[0]); fd.energy_to.setValue(EF_WINDOW[1])
        fd.estimate(); fd.run_fit(); QtWidgets.QApplication.processEvents(); r = fd.fit_result
        print(f"{n:16s} T={T:5.0f} EF={float(r.values['ef']):.5f} +-{float(r.errors['ef']):.5f} "
              f"res={float(r.values['resolution']):.4f} ok={r.success}", flush=True)
        fd.offset_energy_axis(); QtWidgets.QApplication.processEvents(); clear()   # E - E_F axis
    except Exception:
        traceback.print_exc(); clear()
print("rows:", len(A.loaded_items))
