"""
loader -- reading data in, and keeping it.

``registry`` decides who reads a file and applies the load-time axis
options; ``beamlines`` holds one reader per beamline (copy
``beamlines/soleil.py`` to add an HDF5 one); ``formats`` what several of
them share (NeXus layouts, Scienta SES files); ``native`` reads this
program's own saved format; ``nxs_file`` holds the SOLEIL parser, the saved
format and the lazy arrays that keep a measurement on disk; ``session`` is
where a computed dataset is written the moment it exists.

Nothing in here imports Qt.
"""
