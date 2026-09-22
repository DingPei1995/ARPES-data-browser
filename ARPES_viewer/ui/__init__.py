"""
ui -- the windows, and the plumbing only they need.

``main_window`` lays out the launcher; ``windows`` and ``widgets`` are the
viewers and the pieces they are built from; ``figure``, ``fit``, ``process``
and ``volume`` are the panels; ``loader_dialog`` is the Load-data window;
``jobs`` runs one long operation at a time off the GUI thread.

This is the only package that imports Qt.
"""
