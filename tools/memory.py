"""
tools/memory.py
===============
How much memory this process is using, and a best-effort way to give some
of it back.

Kept free of Qt so it can be tested on its own; the launcher's "Memory" box
(``ui/main_window.py``, wired in ``ARPES_viewer.py``) is what shows it.

``psutil`` is used when it is installed, since it answers the question the
same way on every platform. Without it the process figure comes from the
operating system directly (``GetProcessMemoryInfo`` on Windows,
``/proc/self/status`` on Linux, ``getrusage`` elsewhere -- the last being the
*peak*, not the current use, and labelled as such).
"""
from __future__ import annotations

import gc
import os
import sys

try:                                   # optional: the numbers work without it
    import psutil
except ImportError:                    # pragma: no cover - depends on the machine
    psutil = None


def format_bytes(n) -> str:
    """``1536 -> '1.5 KB'``; ``None -> 'n/a'``."""
    if n is None:
        return "n/a"
    n = float(n)
    for unit in ("B", "KB", "MB", "GB"):
        if abs(n) < 1024.0 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.2f} {unit}"
        n /= 1024.0
    return f"{n:.2f} GB"


def _windows_rss():
    import ctypes
    from ctypes import wintypes

    class PROCESS_MEMORY_COUNTERS(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD),
                    ("PageFaultCount", wintypes.DWORD),
                    ("PeakWorkingSetSize", ctypes.c_size_t),
                    ("WorkingSetSize", ctypes.c_size_t),
                    ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                    ("PagefileUsage", ctypes.c_size_t),
                    ("PeakPagefileUsage", ctypes.c_size_t)]

    counters = PROCESS_MEMORY_COUNTERS()
    counters.cb = ctypes.sizeof(counters)
    process = ctypes.windll.kernel32.GetCurrentProcess()
    ok = ctypes.windll.psapi.GetProcessMemoryInfo(process, ctypes.byref(counters),
                                                  counters.cb)
    return int(counters.WorkingSetSize) if ok else None


def _linux_rss():
    with open("/proc/self/status", "r", encoding="ascii", errors="replace") as fh:
        for line in fh:
            if line.startswith("VmRSS:"):
                return int(line.split()[1]) * 1024
    return None


def process_memory():
    """``(bytes, is_peak)`` for this process's resident memory, or
    ``(None, False)`` if it cannot be found out. ``is_peak`` is True only on
    the last-resort path, which can report the high-water mark alone."""
    if psutil is not None:
        try:
            return int(psutil.Process(os.getpid()).memory_info().rss), False
        except Exception:                                   # noqa: BLE001
            pass
    try:
        if sys.platform.startswith("win"):
            return _windows_rss(), False
        if os.path.exists("/proc/self/status"):
            return _linux_rss(), False
    except Exception:                                       # noqa: BLE001
        pass
    try:
        import resource
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # kilobytes on Linux, bytes on macOS
        return (int(peak) if sys.platform == "darwin" else int(peak) * 1024), True
    except Exception:                                       # noqa: BLE001
        return None, False


def system_memory():
    """``(total, available)`` bytes for the whole machine, or ``(None, None)``."""
    if psutil is not None:
        try:
            vm = psutil.virtual_memory()
            return int(vm.total), int(vm.available)
        except Exception:                                   # noqa: BLE001
            pass
    try:
        if sys.platform.startswith("win"):
            import ctypes

            class MEMORYSTATUSEX(ctypes.Structure):
                _fields_ = [("dwLength", ctypes.c_ulong),
                            ("dwMemoryLoad", ctypes.c_ulong),
                            ("ullTotalPhys", ctypes.c_ulonglong),
                            ("ullAvailPhys", ctypes.c_ulonglong),
                            ("ullTotalPageFile", ctypes.c_ulonglong),
                            ("ullAvailPageFile", ctypes.c_ulonglong),
                            ("ullTotalVirtual", ctypes.c_ulonglong),
                            ("ullAvailVirtual", ctypes.c_ulonglong),
                            ("sullAvailExtendedVirtual", ctypes.c_ulonglong)]

            status = MEMORYSTATUSEX()
            status.dwLength = ctypes.sizeof(status)
            if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
                return int(status.ullTotalPhys), int(status.ullAvailPhys)
        elif os.path.exists("/proc/meminfo"):
            values = {}
            with open("/proc/meminfo", "r", encoding="ascii", errors="replace") as fh:
                for line in fh:
                    name, _, rest = line.partition(":")
                    values[name] = int(rest.split()[0]) * 1024
            return values.get("MemTotal"), values.get("MemAvailable")
    except Exception:                                       # noqa: BLE001
        pass
    return None, None


def release_memory() -> int:
    """Collect every reference cycle, then ask the C allocator to hand freed
    pages back to the operating system where it can be asked (glibc's
    ``malloc_trim``). Returns the number of objects the collector freed.

    This only frees what nothing refers to any more; data still held by an
    open window stays. The caller drops its own caches first.
    """
    freed = gc.collect()
    if sys.platform.startswith("linux"):
        try:
            import ctypes
            ctypes.CDLL("libc.so.6").malloc_trim(0)
        except Exception:                                   # noqa: BLE001
            pass
    return freed
