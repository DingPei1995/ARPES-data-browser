"""
tools -- the algorithms, with no Qt in any of them.

Everything a viewer does to data rather than with widgets, grouped by what
it is for:

``momentum``    angle and photon energy to momentum (maps, cuts, k_z)
``crystal``     lattices, space groups, Brillouin zones, moire, cleavage
``processing``  2-D and 3-D operations, combining datasets, de-gridding
``spectra``     EDCs/MDCs, spin, peak and Fermi-edge fits, core levels
``output``      export, the figure model, colormaps
``runtime``     memory and system stand-ins

Each is usable from a script and tested directly, which is why they are
kept apart from ``ui``.
"""
