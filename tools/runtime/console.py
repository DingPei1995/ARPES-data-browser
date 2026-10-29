"""
tools/runtime/console.py
========================
The black console window Windows opens when ``ARPES_viewer.py`` is
double-clicked (``py.exe`` / ``python.exe``): closing it -- easily done by
mistake, it looks like nothing -- closes the program with it, and whatever
was not saved.

``ARPES viewer.pyw`` beside the program is the clean way: Windows runs a
``.pyw`` with ``pythonw.exe``, which opens no console at all. This is the
fallback for a ``.py`` double-clicked all the same: when the console was
opened for this program alone (no command prompt, PowerShell or other shell
is attached to it), it is hidden. Started from a command prompt, the console
is the user's and stays as it is.

Nothing here does anything off Windows: macOS has ``ARPES viewer.app``,
which starts the program with no Terminal window, and Linux desktops start
a ``.py`` without one. No Qt.
"""
from __future__ import annotations

import os
import sys

__all__ = ["hide_own_console", "console_processes", "SHELLS"]

#: a console with one of these attached was opened by the user, not for us
SHELLS = frozenset({"cmd.exe", "powershell.exe", "pwsh.exe", "bash.exe",
                    "sh.exe", "zsh.exe", "wsl.exe", "windowsterminal.exe",
                    "conemu64.exe", "conemuc64.exe", "far.exe", "tcc.exe"})


def _process_name(kernel32, pid: int) -> str:
    """The executable's file name of process ``pid`` (lower case), or ""."""
    from ctypes import byref, create_unicode_buffer, wintypes

    handle = kernel32.OpenProcess(0x1000, False, pid)   # QUERY_LIMITED_INFORMATION
    if not handle:
        return ""
    try:
        size = wintypes.DWORD(1024)
        buffer = create_unicode_buffer(size.value)
        if kernel32.QueryFullProcessImageNameW(handle, 0, buffer, byref(size)):
            return os.path.basename(buffer.value).lower()
        return ""
    finally:
        kernel32.CloseHandle(handle)


def console_processes():
    """The file names of every process attached to this program's console
    (Windows), or None when there is no console or it cannot be told."""
    if sys.platform != "win32":
        return None
    import ctypes
    from ctypes import wintypes

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.OpenProcess.restype = wintypes.HANDLE
    kernel32.OpenProcess.argtypes = (wintypes.DWORD, wintypes.BOOL, wintypes.DWORD)
    kernel32.GetConsoleWindow.restype = wintypes.HWND
    if not kernel32.GetConsoleWindow():
        return None
    pids = (wintypes.DWORD * 64)()
    count = kernel32.GetConsoleProcessList(pids, 64)
    if not 0 < count <= 64:
        return None
    return [_process_name(kernel32, pid) for pid in pids[:count]]


def hide_own_console() -> bool:
    """Hide the console window when Windows opened it for this program
    alone (a double-clicked ``.py``); True when it was hidden. Any doubt
    -- a name that cannot be read, a shell attached -- leaves it shown."""
    try:
        names = console_processes()
        if not names or any(not name or name in SHELLS for name in names):
            return False
        import ctypes

        kernel32 = ctypes.WinDLL("kernel32")
        kernel32.GetConsoleWindow.restype = ctypes.c_void_p
        window = kernel32.GetConsoleWindow()
        user32 = ctypes.WinDLL("user32")
        user32.ShowWindow.argtypes = (ctypes.c_void_p, ctypes.c_int)
        user32.ShowWindow(window, 0)                         # SW_HIDE
        return True
    except Exception:                                      # noqa: BLE001
        return False
