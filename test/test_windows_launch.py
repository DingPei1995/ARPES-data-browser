"""Starting without a console window on Windows: ``ARPES viewer.pyw`` runs
the program with its output in the launcher log, and a double-clicked
``ARPES_viewer.py`` hides the console opened for it alone (never one a
command prompt owns)."""
import os
import shutil
import subprocess
import sys

from arpes_viewer.tools.runtime import console

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _run_pyw(tmp_path, program: str):
    folder = tmp_path / "program"
    folder.mkdir()
    shutil.copy(os.path.join(ROOT, "ARPES viewer.pyw"), folder / "ARPES viewer.pyw")
    (folder / "ARPES_viewer.py").write_text(program)
    home = tmp_path / "home"
    env = dict(os.environ, HOME=str(home), USERPROFILE=str(home))
    result = subprocess.run([sys.executable, str(folder / "ARPES viewer.pyw")],
                            env=env, capture_output=True, text=True, timeout=60)
    log = (home / ".arpes_viewer" / "launcher.log").read_text(encoding="utf-8")
    return result, log


def test_the_pyw_runs_the_program_with_its_output_in_the_log(tmp_path):
    result, log = _run_pyw(tmp_path, "import os, sys\n"
                                     "print('started in', os.getcwd())\n"
                                     "sys.exit(0)\n")
    assert result.returncode == 0 and result.stdout == "" and result.stderr == ""
    assert "started in" in log and "program" in log


def test_the_pyw_reports_an_error_in_the_log(tmp_path):
    result, log = _run_pyw(tmp_path, "import no_such_package_here\n")
    assert result.returncode == 1
    assert "ModuleNotFoundError" in log


def test_off_windows_or_from_a_shell_the_console_stays(monkeypatch):
    if sys.platform != "win32":
        assert console.console_processes() is None
        assert console.hide_own_console() is False
    monkeypatch.setattr(console, "console_processes",
                        lambda: ["cmd.exe", "python.exe"])
    assert console.hide_own_console() is False          # the user's prompt
    monkeypatch.setattr(console, "console_processes", lambda: ["py.exe", ""])
    assert console.hide_own_console() is False          # a name unread: in doubt


def test_errors_are_kept_beside_the_program_and_only_errors(tmp_path):
    import io

    from arpes_viewer.tools.runtime import errorlog

    console_out = io.StringIO()
    tee = errorlog.ErrorTee(console_out, str(tmp_path), "1.2.1")
    tee.write("RuntimeWarning: invalid value\n")
    assert tee.path is None                              # a warning: no file
    tee.write("Traceback (most recent call last):\n  File \"x.py\", line 1\n")
    tee.write("ValueError: boom\n")
    assert tee.path and os.path.dirname(tee.path) == str(tmp_path / "error_logs")
    text = open(tee.path, encoding="utf-8").read()
    assert "ARPES viewer 1.2.1" in text and "RuntimeWarning" in text
    assert "ValueError: boom" in text
    assert "ValueError: boom" in console_out.getvalue()   # the console still has it


def test_an_error_in_a_slot_is_logged_and_the_program_carries_on(tmp_path):
    script = (
        "import sys\n"
        "from PyQt5.QtWidgets import QApplication, QPushButton\n"
        "from PyQt5.QtCore import QTimer\n"
        "app = QApplication([])\n"
        "from arpes_viewer.ui.common import crashguard\n"
        f"crashguard.install({str(tmp_path)!r}, 'test')\n"
        "crashguard._say = lambda *a: None\n"
        "b = QPushButton()\n"
        "def bad(): raise ValueError('a bug in a button')\n"
        "b.clicked.connect(bad)\n"
        "QTimer.singleShot(0, b.click)\n"
        "QTimer.singleShot(200, lambda: (print('still alive'), app.quit()))\n"
        "app.exec_()\n")
    env = dict(os.environ, QT_QPA_PLATFORM="offscreen",
               PYTHONPATH=ROOT + os.pathsep + os.environ.get("PYTHONPATH", ""))
    result = subprocess.run([sys.executable, "-c", script], env=env, cwd=ROOT,
                            capture_output=True, text=True, timeout=60)
    assert "still alive" in result.stdout, result.stderr
    logs = list((tmp_path / "error_logs").glob("errors_*.txt"))
    assert len(logs) == 1
    assert "a bug in a button" in logs[0].read_text(encoding="utf-8")


def test_the_panel_foot_copies_the_newest_error_log(tmp_path, monkeypatch):
    from PyQt5.QtWidgets import QApplication

    from arpes_viewer.tools.runtime import errorlog
    from arpes_viewer.ui.launcher.main_window import about_link

    app = QApplication.instance() or QApplication([])
    monkeypatch.setattr(errorlog, "program_folder", lambda: str(tmp_path))
    assert about_link("errors:copy") == ""                # nothing yet
    folder = errorlog.log_folder()
    assert folder == str(tmp_path / "error_logs")
    path = os.path.join(folder, "errors_2026-10-29_120000.txt")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("ValueError: the bug")
    assert about_link("errors:copy") == path
    assert app.clipboard().text() == "ValueError: the bug"
