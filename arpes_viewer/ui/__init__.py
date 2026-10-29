"""
ui -- the windows, and the plumbing only they need.

``launcher``  the main window, the Load-data window, the list's menu
``viewers``   the per-dataset viewer windows and the widgets they are
              built from
``panels``    the tool windows a viewer opens: processing, fitting, kz,
              combining, core levels, the figure composer
``common``    running long jobs off the GUI thread, messages, GC guard

This is the only package that imports Qt.
"""
