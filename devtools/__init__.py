"""
devtools -- one-off generators, run by hand and never by the program.

``gen_colormaps`` needs scipy and the lab's ``Colormap.mat``;
``gen_spacegroups`` needs ``spglib``. Both write a self-contained table into
``tools/``, which is what the running program imports -- so neither
dependency is one the program has.
"""
