"""Live ARPES viewer that an agent can drive from files (no GUI control needed).

Start with launch_viewer.py. Reads agentic/config.json (copy config.example.json):
  data_dir  folder with the measurement files / scan folders
  patterns  glob patterns relative to data_dir (default matches CASSIOPEE .txt cuts
            and <scan folder>/*_roi1_.txt maps)
  screen    1-based monitor number for the window (default: 2)
Every *.py dropped into agentic/viewer_cmd/ is executed inside this running GUI process
(globals: A = arpes_viewer.app, HERE, DATA, QtCore, QtWidgets); its printed output goes to
agentic/out/cmd_<name>.log and the file is moved to agentic/viewer_done/.
Nothing here uses the network.
"""
import os, sys, glob, json, time, io, contextlib, traceback
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                       # the repository
CFG = json.load(open(os.path.join(HERE, "config.json"), encoding="utf-8"))
DATA = CFG["data_dir"]; PATTERNS = CFG.get("patterns", ["*.txt", "*/*_roi1_.txt"])
sys.path.insert(0, ROOT)
os.environ["QT_QPA_PLATFORM"] = "windows"
import arpes_viewer.app as A
from PyQt5 import QtCore, QtWidgets

paths = []
for pat in PATTERNS:
    paths += sorted(glob.glob(os.path.join(DATA, pat)))

def label_of(src):                                  # row name = file name, or folder name for scans
    return os.path.splitext(os.path.basename(src))[0] if os.path.dirname(src) == DATA \
        else os.path.basename(os.path.dirname(src))

def setup():
    try:
        w = A.win
        scr = QtWidgets.QApplication.screens()
        want = int(CFG.get("screen", 2)) - 1
        if len(scr) > 1 and 0 <= want < len(scr):
            g = scr[want].geometry(); w.move(g.x() + 40, g.y() + 40)
        A.add_files(paths, background=False)
        for key, rec in A.loaded_items.items():
            rec["name"] = label_of(rec["path"])
            item = A._row_for_key(key)
            if item is not None: item.setText(A._row_label(key))
        w.show(); w.raise_(); QtWidgets.QApplication.processEvents()
        time.sleep(0.5); QtWidgets.QApplication.processEvents()
        w.grab().save(os.path.join(HERE, "shots", "step01_main_panel_loaded.png"))
        print("listed", len(A.loaded_items), "datasets", flush=True)
    except Exception:
        traceback.print_exc(); sys.stdout.flush()

def poll():
    qd, dd = os.path.join(HERE, "viewer_cmd"), os.path.join(HERE, "viewer_done")
    for n in sorted(os.listdir(qd)):
        if not n.endswith(".py"): continue
        src = os.path.join(qd, n); time.sleep(0.2)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            try:
                exec(compile(open(src, encoding="utf-8").read(), src, "exec"),
                     {"A": A, "HERE": HERE, "DATA": DATA, "QtCore": QtCore, "QtWidgets": QtWidgets})
            except Exception:
                traceback.print_exc()
        open(os.path.join(HERE, "out", "cmd_" + n[:-3] + ".log"), "w", encoding="utf-8").write(buf.getvalue())
        os.replace(src, os.path.join(dd, n))

_orig = A.offer_recovery
def patched():
    QtCore.QTimer.singleShot(300, setup)
    t = QtCore.QTimer(A.win); t.timeout.connect(poll); t.start(700); A._poll_timer = t
A.offer_recovery = patched
sys.exit(A.main())
