"""
tools/runtime/errorlog.py
=========================
Errors kept as a text file in the program's folder, by themselves.

What the program writes to the console when something goes wrong -- a
traceback, an error message -- is easily lost: the console is hidden (a
double-clicked ``.py`` on Windows, ``ARPES viewer.pyw``, ``ARPES
viewer.app``), closed, or scrolled away. :func:`install` copies the error
stream into ``error_logs/errors_<date>_<time>.txt`` in the program's folder
(the home folder's ``.arpes_viewer`` when that one is read-only), to send
with a bug report. The file is made only when an error comes: a session
without one leaves nothing behind. What came just before the error on the
same stream (a warning, say) is written with it, for context. The console,
where there is one, still shows everything as before.

:func:`report` writes an exception there (the uncaught-exception hook calls
it). No Qt.
"""
from __future__ import annotations

import datetime
import os
import sys
import traceback

__all__ = ["install", "report", "ErrorTee", "FOLDER_NAME", "program_folder",
           "log_folder", "newest"]

FOLDER_NAME = "error_logs"
#: what in a line says an error has come (a warning alone makes no file)
_MARKS = ("Traceback (most recent call last)", "Error", "Exception", "FATAL",
          "Fatal", "error:")
#: lines kept in memory before the first error, written with it as context
_CONTEXT = 40

_tee = None


def _folders(program_folder: str):
    yield os.path.join(program_folder, FOLDER_NAME)
    yield os.path.join(os.path.expanduser("~"), ".arpes_viewer", FOLDER_NAME)


class ErrorTee:
    """A text stream that passes everything to ``stream`` (if any) and,
    from the first error on, also to a file of its own."""

    def __init__(self, stream, program_folder: str, version: str = ""):
        self.stream = stream
        self.program_folder = program_folder
        self.version = version
        self.path = None
        self._file = None
        self._before = []
        self._partial = ""

    # -- the stream interface ----------------------------------------------
    def write(self, text):
        if self.stream is not None:
            try:
                self.stream.write(text)
            except Exception:                              # noqa: BLE001
                pass
        if not text:
            return 0
        if self._file is None:
            self._partial += text
            lines = self._partial.split("\n")
            self._partial = lines.pop()
            self._before.extend(lines)
            if any(mark in line for line in lines + [self._partial] for mark in _MARKS):
                self._open()
            else:
                del self._before[:-_CONTEXT]
        else:
            self._file.write(text)
            self._file.flush()
        return len(text)

    def flush(self):
        if self.stream is not None:
            try:
                self.stream.flush()
            except Exception:                              # noqa: BLE001
                pass
        if self._file is not None:
            self._file.flush()

    def isatty(self):
        return False

    def fileno(self):
        if self.stream is not None:
            return self.stream.fileno()
        raise OSError("no file descriptor")

    @property
    def encoding(self):
        return getattr(self.stream, "encoding", "utf-8") or "utf-8"

    # -- the file -------------------------------------------------------------
    def _open(self):
        stamp = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")
        for folder in _folders(self.program_folder):
            try:
                os.makedirs(folder, exist_ok=True)
                path = os.path.join(folder, f"errors_{stamp}.txt")
                self._file = open(path, "a", encoding="utf-8")
                self.path = path
                break
            except OSError:
                continue
        if self._file is None:
            return
        self._file.write(
            f"ARPES viewer {self.version} -- errors of the session started "
            f"{stamp.replace('_', ' ')}\nPython {sys.version.split()[0]} on "
            f"{sys.platform}; {sys.executable}\n{'-' * 72}\n")
        self._file.write("\n".join(self._before) + ("\n" if self._before else ""))
        self._file.write(self._partial)
        self._before, self._partial = [], ""
        self._file.flush()


def install(program_folder: str, version: str = "") -> ErrorTee:
    """Copy the error stream into ``program_folder``'s error log from now
    on; returns the tee (its ``path`` once an error has come)."""
    global _tee
    if isinstance(sys.stderr, ErrorTee):
        _tee = sys.stderr
        return _tee
    _tee = ErrorTee(sys.stderr, program_folder, version)
    sys.stderr = _tee
    return _tee


def report(exc_type, exc, tb, where: str = "") -> str:
    """Write an exception to the error stream (so to the log); returns the
    log's path ("" when there is none)."""
    stamp = datetime.datetime.now().strftime("%H:%M:%S")
    text = "".join(traceback.format_exception(exc_type, exc, tb))
    stream = sys.stderr
    if stream is not None:
        stream.write(f"\n[{stamp}] Error{(' in ' + where) if where else ''}:\n{text}")
        stream.flush()
    return getattr(_tee, "path", "") or ""


def program_folder() -> str:
    """The folder ``ARPES_viewer.py`` is in."""
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def log_folder(folder: str = None) -> str:
    """Where the error logs go: ``error_logs`` in the program's folder, or
    the home folder's when that one cannot be written. Made if missing."""
    for candidate in _folders(folder or program_folder()):
        try:
            os.makedirs(candidate, exist_ok=True)
            if os.access(candidate, os.W_OK):
                return candidate
        except OSError:
            continue
    return ""


def newest(folder: str = None) -> str:
    """The most recent error log (of either place), or ""."""
    found = []
    for candidate in _folders(folder or program_folder()):
        try:
            found += [os.path.join(candidate, name) for name in os.listdir(candidate)
                      if name.startswith("errors_") and name.endswith(".txt")]
        except OSError:
            continue
    return max(found, key=os.path.getmtime) if found else ""
