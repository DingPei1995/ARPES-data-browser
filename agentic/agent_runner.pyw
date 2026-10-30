# Local job runner: executes *.py jobs dropped into ./queue with THIS machine's Python.
# Standard library only. No network access, no uploads. Stop by creating a file named STOP here.
# Start by double-clicking (associate .pyw with your Python once).
import os, sys, time, subprocess, shutil, datetime
HERE = os.path.dirname(os.path.abspath(__file__))
Q, OUT, DONE = (os.path.join(HERE, d) for d in ("queue", "out", "done"))
py = sys.executable
cand = os.path.join(os.path.dirname(py), "python.exe")
if py.lower().endswith("pythonw.exe") and os.path.exists(cand):
    py = cand
def beat(msg="alive"):
    with open(os.path.join(HERE, "heartbeat.txt"), "w") as f:
        f.write(f"{datetime.datetime.now().isoformat()} {msg} | {py}\n")
stop = os.path.join(HERE, "STOP")
if os.path.exists(stop): os.remove(stop)
beat("started")
while not os.path.exists(stop):
    beat()
    for n in sorted(os.listdir(Q)):
        if not n.endswith(".py"): continue
        src = os.path.join(Q, n)
        time.sleep(0.3)
        log = os.path.join(OUT, n[:-3] + ".log")
        env = dict(os.environ, QT_QPA_PLATFORM=os.environ.get("QT_QPA_PLATFORM", "offscreen"),
                   PYTHONIOENCODING="utf-8", HF_HUB_OFFLINE="1")
        with open(log, "w", encoding="utf-8") as lf:
            t0 = time.time()
            try:
                r = subprocess.run([py, src], stdout=lf, stderr=subprocess.STDOUT, cwd=HERE,
                                   env=env, creationflags=0x08000000, timeout=3600)
                lf.write(f"\n[exit {r.returncode}, {time.time()-t0:.1f}s]\n")
            except Exception as e:
                lf.write(f"\n[runner error] {e!r}\n")
        shutil.move(src, os.path.join(DONE, n))
    time.sleep(1)
beat("stopped")
