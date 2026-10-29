"""Start ARPES viewer from a copy of the repository:

    python ARPES_viewer.py

(what ``ARPES viewer.pyw`` on Windows and ``ARPES viewer.app`` on macOS run
too). The program itself is the ``arpes_viewer`` package beside this file;
an installed copy (``pip install arpes-viewer``) starts with the
``arpes-viewer`` command instead.

Copyright (c) 2026 Ding Pei. Licensed under the PolyForm Noncommercial
License 1.0.0 (LICENSE.md).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

if __name__ == "__main__":
    from arpes_viewer.app import main

    sys.exit(main())
