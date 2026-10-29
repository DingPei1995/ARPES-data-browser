"""
ARPES viewer.pyw -- start the program with a double-click on Windows,
without the console window.

Windows runs a ``.pyw`` file with ``pythonw.exe``, which opens no console:
there is no black window to close by mistake (closing it closed the
program). What the program prints goes to ``~/.arpes_viewer/launcher.log``
instead, as with ``ARPES viewer.app`` on macOS; if it stops with an error, a
message box shows the end of that log.

To start it from a command prompt, with its output there, run
``python ARPES_viewer.py`` as before.
"""
import os
import runpy
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(os.path.expanduser("~"), ".arpes_viewer")
LOG = os.path.join(CONFIG, "launcher.log")


def _say(text: str):
    """A message box with no Qt (Qt may be what failed to load)."""
    try:
        import ctypes

        ctypes.windll.user32.MessageBoxW(None, text, "ARPES viewer", 0x10)
    except Exception:                                      # noqa: BLE001
        pass


def main() -> int:
    os.makedirs(CONFIG, exist_ok=True)
    log = open(LOG, "w", encoding="utf-8", buffering=1)
    # pythonw has no stdout / stderr: give the program the log for both
    sys.stdout = sys.stderr = log
    print(f"=== {sys.executable} -- {HERE}")
    os.chdir(HERE)
    sys.path.insert(0, HERE)
    sys.argv = [os.path.join(HERE, "ARPES_viewer.py")] + sys.argv[1:]
    code = 0
    try:
        runpy.run_path(sys.argv[0], run_name="__main__")
    except SystemExit as exit_:
        code = exit_.code if isinstance(exit_.code, int) else (0 if exit_.code is None else 1)
    except BaseException:                                  # noqa: BLE001
        traceback.print_exc()
        code = 1
    log.flush()
    if code:
        try:
            with open(LOG, encoding="utf-8", errors="replace") as handle:
                tail = "".join(handle.readlines()[-12:])
        except OSError:
            tail = ""
        _say(f"ARPES viewer stopped with an error (code {code}). The end of its "
             f"log ({LOG}):\n\n{tail}")
    return code


if __name__ == "__main__":
    sys.exit(main())
