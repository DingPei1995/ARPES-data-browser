"""
loader -- reading data in, and keeping it.

``registry`` decides who reads a file and applies the load-time axis
options; ``soleil`` and ``native`` are the two readers that ship (copy
``soleil`` to add a beamline); ``nxs_file`` holds the SOLEIL parser, this
program's own saved format and the lazy arrays that keep a measurement on
disk; ``session`` is where a computed dataset is written the moment it
exists.

Nothing in here imports Qt.
"""
