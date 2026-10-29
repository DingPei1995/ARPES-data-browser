"""The 3-D views' Qt-free geometry (tools/processing/volume.py): which box
edge each coordinate axis is drawn on, where it starts, and cropping a
lazy cube before it is read."""
import numpy as np

import arpes_viewer.tools.processing.volume as V

BOUNDS = [(-10.0, 10.0), (-5.0, 5.0), (54.0, 56.0)]


def test_an_axis_runs_along_its_edge_from_the_end_it_is_told():
    low, high = V.edge_ends(BOUNDS, 0, (1, -1))
    np.testing.assert_allclose(low, [-10, 5, 54])
    np.testing.assert_allclose(high, [10, 5, 54])
    (d, start, end, ticks, step), = V.axis_segments(
        BOUNDS, [(1, -1)] * 3, [-1, 1, 1], shown=(True, False, False))
    assert d == 0 and start[0] == 10 and end[0] == -10
    assert ticks.min() >= -10 and ticks.max() <= 10 and step > 0


def test_automatic_edges_are_on_the_outline_facing_the_viewer():
    camera = V.Camera(azimuth=45.0, elevation=25.0)
    edges = V.auto_edges(camera, BOUNDS)
    view = camera.basis()[2]
    for axis, edge in enumerate(edges):
        others = [d for d in range(3) if d != axis]
        facing = [sign * view[d] > 0 for d, sign in zip(others, edge)]
        assert facing[0] != facing[1]                  # one face each way
    # Looking down from above, the two in-plane axes sit on the bottom
    # (z min) of the box, not across its top.
    assert edges[0][1] == -1 and edges[1][1] == -1


class _Lazy:
    """Slices like an array and records what was read."""

    def __init__(self, array):
        self.array, self.shape, self.reads = array, array.shape, []

    def __getitem__(self, key):
        self.reads.append(key)
        return self.array[key]

    def __array__(self, dtype=None):
        raise AssertionError("the whole cube was read")


def test_a_lazy_cube_is_cropped_before_it_is_read():
    axes = [np.linspace(0, 1, 20), np.linspace(0, 1, 30), np.linspace(0, 1, 40)]
    cube = _Lazy(np.random.default_rng(0).random((20, 30, 40)))
    values, out = V.reduce_volume(cube, axes, [(0, 0.5), None, None], [5, 30, 40])
    assert len(cube.reads) == 1 and values.shape == (5, 30, 40)
    assert out[0].max() <= 0.5
