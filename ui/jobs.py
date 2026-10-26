"""
ui/jobs.py
===========
Running one long operation off the GUI thread, with a progress bar that
moves and a Cancel that works.

Why this exists
---------------
Every computation in this program used to run on the GUI thread. Reading a
map, converting to k-space, resampling a cube for the 3-D view: while any of
them ran, the event loop was blocked, so the window stopped repainting and
the desktop marked it "not responding". From the outside that is
indistinguishable from a crash, and the natural response -- kill it -- threw
away everything computed in the session. The freezing and the data loss were
the same bug wearing two hats.

Four places had patched around it with ``QApplication.processEvents()``
inside their own loops. That keeps the window painting but is not
concurrency: it re-enters the event loop from inside the computation, so a
click on another button runs *that* handler on top of the half-finished one.
This module replaces that pattern with a real worker thread.

The rules it enforces
---------------------
* **One job at a time.** HDF5 here is not built thread-safe, and several
  jobs reading the same file from several threads is a crash rather than a
  slowdown. :func:`run_job` refuses to start a second job while one is
  running, and says so.
* **The worker never touches Qt.** A job function gets a plain
  ``report(fraction, message)`` callable and returns a plain object; the
  result is delivered back on the GUI thread by signal. Anything that builds
  a widget must happen in ``on_done``.
* **Cancel is cooperative.** ``report`` raises :class:`JobCancelled` once the
  user has pressed Cancel, which unwinds the job function wherever it
  happens to be. A job that never calls ``report`` cannot be cancelled --
  that is the job's fault, not the runner's, and is worth knowing when
  writing one.

A job that finishes faster than :data:`DIALOG_DELAY_MS` never shows a
dialog, so the common quick case does not flash a window.
"""
from __future__ import annotations

import threading
import traceback

from PyQt5.QtCore import QObject, QThread, Qt, QTimer, pyqtSignal
from PyQt5.QtWidgets import QApplication, QProgressDialog

__all__ = ["JobCancelled", "JobRunner", "run_job", "run_blocking", "busy",
           "HDF5_LOCK"]

#: Serialises every read through h5py. The handles in
#: ``loader.nxs_file.HANDLES`` are shared between the GUI thread (a viewer
#: asking for one frame) and a worker (an operation reading a whole cube),
#: and HDF5 is not safe against that on its own.
HDF5_LOCK = threading.RLock()

#: How long a job may run before its progress dialog appears.
DIALOG_DELAY_MS = 400


class JobCancelled(Exception):
    """Raised inside a job's own thread when Cancel has been pressed."""


class _Worker(QObject):
    finished = pyqtSignal(object)       # the job's return value
    failed = pyqtSignal(object)         # the exception
    cancelled = pyqtSignal()
    progressed = pyqtSignal(float, str)  # 0..1, message

    def __init__(self, work):
        super().__init__()
        self._work = work
        self._cancel = threading.Event()

    def cancel(self):
        self._cancel.set()

    def report(self, fraction: float = None, message: str = ""):
        """What a job calls to say where it is -- and the only place a
        cancellation can take effect."""
        if self._cancel.is_set():
            raise JobCancelled()
        self.progressed.emit(-1.0 if fraction is None else float(fraction),
                             str(message))

    def run(self):
        try:
            result = self._work(self.report)
        except JobCancelled:
            self.cancelled.emit()
        except Exception as exc:                    # noqa: BLE001 -- reported
            traceback.print_exc()
            self.failed.emit(exc)
        else:
            self.finished.emit(result)


class JobRunner(QObject):
    """Owns the one worker thread that may be running.

    Kept as an object rather than loose functions so that a test (or a
    headless script) can make its own and drive it without touching global
    state.
    """

    def __init__(self):
        super().__init__()
        self._thread = None
        self._worker = None
        self._dialog = None
        self._timer = None

    def is_busy(self) -> bool:
        return self._thread is not None and self._thread.isRunning()

    def run(self, parent, title: str, work, on_done=None, on_error=None,
            on_cancel=None, cancellable: bool = True):
        """Start ``work`` on the worker thread.

        ``work`` is called as ``work(report)`` off the GUI thread and must
        not touch any widget. ``on_done(result)``, ``on_error(exc)`` and
        ``on_cancel()`` are called back *on* the GUI thread.

        Returns False (without starting anything) if a job is already
        running, since two at once is the unsafe case this class exists to
        prevent.
        """
        if self.is_busy():
            return False

        self._thread = QThread()
        self._worker = _Worker(work)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)

        dialog = QProgressDialog(title, "Cancel" if cancellable else None,
                                 0, 100, parent)
        dialog.setWindowTitle(title)
        dialog.setWindowModality(Qt.WindowModal)
        dialog.setAutoClose(False)
        dialog.setAutoReset(False)
        dialog.setMinimumDuration(0)
        dialog.reset()                    # do not show it yet
        if cancellable:
            dialog.canceled.connect(self._worker.cancel)
        self._dialog = dialog

        # Only show the bar for something that is actually taking a while.
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.timeout.connect(self._show_dialog)
        self._timer.start(DIALOG_DELAY_MS)

        self._worker.progressed.connect(self._on_progress)
        # Tear down first, then call back: the callback routinely opens a
        # window or starts the next job, and it must not find this one still
        # holding the thread and the progress dialog.
        self._worker.finished.connect(lambda result: self._deliver(on_done, (result,)))
        self._worker.failed.connect(lambda exc: self._deliver(on_error, (exc,)))
        self._worker.cancelled.connect(lambda: self._deliver(on_cancel, ()))
        self._thread.start()
        return True

    # -- plumbing -----------------------------------------------------------
    def _show_dialog(self):
        if self.is_busy() and self._dialog is not None:
            self._dialog.show()

    def _on_progress(self, fraction: float, message: str):
        dialog = self._dialog
        if dialog is None:
            return
        if fraction < 0:
            dialog.setRange(0, 0)         # indeterminate: a busy sweep
        else:
            dialog.setRange(0, 100)
            dialog.setValue(int(max(0.0, min(1.0, fraction)) * 100))
        if message:
            dialog.setLabelText(message)

    def _deliver(self, callback, args: tuple):
        self._teardown()
        if callback is not None:
            callback(*args)

    def _teardown(self):
        if self._timer is not None:
            self._timer.stop()
            self._timer = None
        if self._dialog is not None:
            self._dialog.close()
            self._dialog = None
        if self._thread is not None:
            self._thread.quit()
            self._thread.wait(5000)
            self._thread = None
        self._worker = None

    def run_blocking(self, parent, title: str, work, cancellable: bool = True):
        """Run ``work`` on the worker thread and return its result.

        The caller waits, but the window does not freeze: the event loop is
        pumped while the worker runs, so the application keeps painting and
        the progress bar moves. The dialog is **application**-modal for the
        duration, which is what makes pumping the loop safe here -- with
        every other widget refusing input, the re-entrancy that makes
        ``processEvents`` inside a computation dangerous has nowhere to
        happen. The only thing that can be clicked is Cancel.

        This exists for the one shape :meth:`run` cannot serve: a value that
        a synchronous caller needs *now* and that eight call sites expect to
        be returned rather than delivered to a callback. Reading a CASSIOPEE
        folder is the case -- a hundred text files, fifteen seconds -- and
        rewriting every caller to be asynchronous to accommodate it would be
        a far larger change than the problem justifies.

        Raises :class:`JobCancelled` if the user cancels, and re-raises
        whatever the job raised.
        """
        if self.is_busy():
            raise RuntimeError(
                "another operation is already running; wait for it to finish")

        outcome = {}
        worker = _Worker(work)
        thread = QThread()
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        # Held on the runner, not just locally: the event loop is being
        # pumped below, so a timer or a stray click can run another handler
        # while this job is in flight. is_busy() has to say yes to it, or
        # that handler starts a second job and the one-at-a-time rule that
        # keeps HDF5 safe is broken by the very mechanism that keeps the
        # window alive.
        self._thread = thread
        self._worker = worker

        dialog = QProgressDialog(title, "Cancel" if cancellable else None,
                                 0, 100, parent)
        dialog.setWindowTitle(title)
        dialog.setWindowModality(Qt.ApplicationModal)
        dialog.setAutoClose(False)
        dialog.setAutoReset(False)
        dialog.setMinimumDuration(0)
        dialog.reset()
        if cancellable:
            dialog.canceled.connect(worker.cancel)

        def progressed(fraction, message):
            if fraction < 0:
                dialog.setRange(0, 0)
            else:
                dialog.setRange(0, 100)
                dialog.setValue(int(max(0.0, min(1.0, fraction)) * 100))
            if message:
                dialog.setLabelText(message)

        worker.progressed.connect(progressed)
        worker.finished.connect(lambda value: outcome.update(result=value))
        worker.failed.connect(lambda exc: outcome.update(error=exc))
        worker.cancelled.connect(lambda: outcome.update(error=JobCancelled()))

        timer = QTimer(self)
        timer.setSingleShot(True)
        timer.timeout.connect(
            lambda: dialog.show() if thread.isRunning() else None)
        timer.start(DIALOG_DELAY_MS)

        thread.start()
        try:
            # Waiting on the *outcome*, not on the thread. A QThread runs an
            # event loop, so it stays "running" long after the job function
            # has returned -- until quit() below. Waiting for it to stop
            # would deadlock: quit() only comes after the wait.
            while not outcome:
                QApplication.processEvents()
                if not thread.isRunning():
                    # It ended without emitting anything, which _Worker.run
                    # does not do. Rather than spin forever, drain whatever
                    # is queued and give up.
                    QApplication.processEvents()
                    break
                # A 20 ms sleep, not a wait for completion (which never
                # comes): without it this is a busy loop burning a core for
                # the length of the job.
                thread.wait(20)
            QApplication.processEvents()
        finally:
            timer.stop()
            dialog.close()
            thread.quit()
            thread.wait(5000)
            self._thread = None
            self._worker = None

        if "error" in outcome:
            raise outcome["error"]
        return outcome.get("result")

    def wait(self, milliseconds: int = 30000) -> bool:
        """Block until the running job is done. For tests and for shutdown,
        never for ordinary use -- blocking is the thing this avoids."""
        thread = self._thread
        if thread is None:
            return True
        deadline = milliseconds
        while thread.isRunning() and deadline > 0:
            QApplication.processEvents()
            thread.wait(20)
            deadline -= 20
        QApplication.processEvents()
        return not thread.isRunning()


#: The application's runner. One per process, because the "one job at a
#: time" rule is a property of the HDF5 library, not of any one window.
_RUNNER = None


def runner() -> JobRunner:
    global _RUNNER
    if _RUNNER is None:
        _RUNNER = JobRunner()
    return _RUNNER


def run_job(parent, title: str, work, on_done=None, on_error=None,
            on_cancel=None, cancellable: bool = True) -> bool:
    """Run ``work`` in the background; see :meth:`JobRunner.run`."""
    return runner().run(parent, title, work, on_done=on_done, on_error=on_error,
                        on_cancel=on_cancel, cancellable=cancellable)


def run_blocking(parent, title: str, work, cancellable: bool = True):
    """Run ``work`` off the GUI thread and return its result; see
    :meth:`JobRunner.run_blocking`."""
    return runner().run_blocking(parent, title, work, cancellable=cancellable)


def busy() -> bool:
    return runner().is_busy()
