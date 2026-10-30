# Starts viewer_live.py detached (survives this launcher). Run: python launch_viewer.py
import os, sys, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
py = sys.executable.replace("pythonw.exe", "python.exe")
env = dict(os.environ, QT_QPA_PLATFORM="windows", PYTHONIOENCODING="utf-8")
lf = open(os.path.join(HERE, "out", "viewer_live.log"), "w", encoding="utf-8")
subprocess.Popen([py, os.path.join(HERE, "viewer_live.py")], stdout=lf, stderr=subprocess.STDOUT, env=env,
                 cwd=HERE, creationflags=0x00000008 | 0x00000200 | 0x08000000)   # detached, new group, no console
print("launched")
