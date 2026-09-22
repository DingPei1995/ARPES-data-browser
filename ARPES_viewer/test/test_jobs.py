"""Tests for ui/jobs.py: that a long operation really does run off the GUI
thread, that Cancel unwinds it, that a failure comes back as a failure, and
that two jobs cannot run at once (which is the rule that keeps HDF5 safe).

Needs a Qt application; run headless with QT_QPA_PLATFORM=offscreen.
"""
import os
import threading
import time

import pytest

# Before any QApplication exists: on a machine with no display (a build
# server, an ssh session) Qt otherwise aborts the whole interpreter rather
# than raising, which takes the rest of the suite down with it.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt5")
from PyQt5.QtWidgets import QApplication, QWidget      # noqa: E402

from ui.jobs import JobRunner, JobCancelled           # noqa: E402


@pytest.fixture(scope="module")
def app():
    application = QApplication.instance() or QApplication([])
    yield application


@pytest.fixture()
def runner(app):
    made = JobRunner()
    yield made
    made.wait(5000)


@pytest.fixture()
def parent(app):
    return QWidget()


def test_a_job_runs_and_returns_its_result(runner, parent):
    seen = {}
    assert runner.run(parent, "Working", lambda report: 6 * 7,
                      on_done=lambda value: seen.update(value=value))
    assert runner.wait()
    assert seen["value"] == 42


def test_the_job_runs_on_another_thread(runner, parent):
    """The whole point: the GUI thread must be free while it runs."""
    seen = {}

    def work(report):
        return threading.current_thread().ident

    runner.run(parent, "Working", work, on_done=lambda ident: seen.update(ident=ident))
    runner.wait()
    assert seen["ident"] != threading.current_thread().ident


def test_progress_reports_reach_the_gui_thread(runner, parent):
    reports = []

    def work(report):
        for i in range(4):
            report(i / 4.0, f"step {i}")
        return None

    runner.run(parent, "Working", work)
    runner._worker.progressed.connect(lambda f, m: reports.append((f, m)))
    runner.wait()
    assert reports                       # at least the later ones land
    assert all(0.0 <= f <= 1.0 for f, _m in reports)


def test_a_second_job_is_refused_while_one_runs(runner, parent):
    def slow(report):
        for _ in range(40):
            report(0.5)
            time.sleep(0.01)
        return "slow"

    assert runner.run(parent, "Slow", slow) is True
    assert runner.run(parent, "Other", lambda report: 1) is False
    assert runner.wait()


def test_cancel_unwinds_the_job(runner, parent):
    state = {}

    def forever(report):
        for _ in range(10000):
            report(0.1)
            time.sleep(0.002)
        return "finished after all"

    runner.run(parent, "Cancel me", forever,
               on_done=lambda value: state.update(done=value),
               on_cancel=lambda: state.update(cancelled=True))
    QApplication.processEvents()
    time.sleep(0.05)
    runner._worker.cancel()
    assert runner.wait()
    assert state.get("cancelled") is True
    assert "done" not in state


def test_a_job_that_never_reports_cannot_be_cancelled_but_still_finishes(runner, parent):
    """Cancellation is cooperative, and this documents the consequence: a
    job that never calls report() runs to the end. Worth knowing when
    writing one."""
    def uninterruptible(report):
        time.sleep(0.05)
        return "done anyway"

    seen = {}
    runner.run(parent, "Busy", uninterruptible, on_done=lambda v: seen.update(v=v))
    runner._worker.cancel()
    runner.wait()
    assert seen["v"] == "done anyway"


def test_a_failure_is_delivered_not_raised(runner, parent):
    seen = {}

    def boom(report):
        raise ValueError("kaboom")

    runner.run(parent, "Boom", boom, on_error=lambda exc: seen.update(exc=exc))
    assert runner.wait()
    assert isinstance(seen["exc"], ValueError)
    assert "kaboom" in str(seen["exc"])


def test_the_runner_is_idle_again_afterwards(runner, parent):
    runner.run(parent, "Quick", lambda report: None)
    runner.wait()
    assert not runner.is_busy()
    # and can take another job
    assert runner.run(parent, "Again", lambda report: 2) is True
    runner.wait()


# --------------------------------------------------------------------------
# run_blocking: the caller waits, the window does not freeze
# --------------------------------------------------------------------------
def test_run_blocking_returns_the_value(runner, parent):
    assert runner.run_blocking(parent, "Working", lambda report: 6 * 7) == 42


def test_run_blocking_still_runs_off_the_gui_thread(runner, parent):
    """"Blocking" means the caller waits for a value, not that the work
    happens on the event-loop thread -- which is what would freeze the
    window and is the thing this module exists to stop."""
    ident = runner.run_blocking(
        parent, "Working", lambda report: threading.current_thread().ident)
    assert ident != threading.current_thread().ident


def test_run_blocking_pumps_the_event_loop(runner, parent):
    """The window keeps painting while the job runs. Checked by seeing that
    a timer queued before the job fires during it."""
    from PyQt5.QtCore import QTimer

    fired = []
    QTimer.singleShot(10, lambda: fired.append(True))

    def slow(report):
        time.sleep(0.3)
        return "done"

    assert runner.run_blocking(parent, "Slow", slow) == "done"
    assert fired, "the event loop was not running while the job was"


def test_run_blocking_reports_progress(runner, parent):
    def work(report):
        for step in range(4):
            report(step / 4.0, f"step {step}")
        return "done"

    assert runner.run_blocking(parent, "Stepping", work) == "done"


def test_run_blocking_re_raises_the_jobs_failure(runner, parent):
    def boom(report):
        raise ValueError("kaboom")

    with pytest.raises(ValueError, match="kaboom"):
        runner.run_blocking(parent, "Boom", boom)


def test_run_blocking_raises_when_cancelled(runner, parent):
    """Cancel has to reach the caller, or it would get a half-built object
    back and treat it as a successful load."""
    def work(report):
        for step in range(1000):
            report(step / 1000.0, "working")
            time.sleep(0.005)
        return "should not get here"

    from PyQt5.QtCore import QTimer
    QTimer.singleShot(50, lambda: runner._worker.cancel())
    with pytest.raises(JobCancelled):
        runner.run_blocking(parent, "Cancel me", work)


def test_the_runner_is_busy_during_a_blocking_job(runner, parent):
    """The event loop is pumped while it runs, so another handler can fire.
    If is_busy() said no, that handler could start a second job -- breaking
    the one-at-a-time rule by way of the mechanism that keeps the window
    responsive."""
    from PyQt5.QtCore import QTimer

    seen = {}
    QTimer.singleShot(30, lambda: seen.update(
        busy=runner.is_busy(),
        second_job_refused=(runner.run(parent, "Nope", lambda r: 1) is False)))

    def slow(report):
        time.sleep(0.25)
        return "done"

    runner.run_blocking(parent, "Slow", slow)
    assert seen.get("busy") is True
    assert seen.get("second_job_refused") is True


def test_run_blocking_leaves_the_runner_idle(runner, parent):
    runner.run_blocking(parent, "Quick", lambda report: 1)
    assert not runner.is_busy()
    assert runner.run_blocking(parent, "Again", lambda report: 2) == 2
