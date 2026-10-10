"""Electrode geometry must not depend on the lead direction beyond a rigid motion."""

import numpy as np
import pytest

from ossdbs.electrodes import ELECTRODES
from ossdbs.electrodes.utilities import get_rotation_axis

# Directions that broke meshing when primitives were built along the lead:
# nearly in the xz-plane (OCC's default cylinder seam then lies within a
# fraction of a degree of the rotated contact seams), axis-aligned, reversed.
DIRECTIONS = [
    (-0.27710985281597844, -0.0028000995593097785, 0.9608341630660125),
    (-0.4565710794209321, 0.2172862364392653, 0.8627453511265445),
    (0.3, 1e-3, 0.95),
    (1.0, 0.0, 0.0),
    (0.0, 0.0, -1.0),
]
POSITION = (0.794, -32.671, -0.618)
THICKNESS = 0.1


def _signature(name, direction):
    electrode = ELECTRODES[name](direction=direction, position=POSITION)
    shapes = [electrode.geometry]
    try:
        shapes.append(electrode.encapsulation_geometry(THICKNESS))
    except NotImplementedError:
        pass
    topology = [(len(s.faces), len(s.edges), len(s.vertices)) for s in shapes]
    volumes = [s.mass for s in shapes]
    return topology, volumes


@pytest.mark.parametrize(
    ("direction", "expected"),
    [
        ((0, 0, 1), (1.0, 0.0, 0.0)),
        ((0, 0, -1), (1.0, 0.0, 0.0)),
        ((1, 0, 0), (0.0, 1.0, 0.0)),
        ((0, 1, 0), (-1.0, 0.0, 0.0)),
    ],
)
def test_get_rotation_axis(direction, expected):
    """Axis is finite and unit length, also for (anti)parallel directions."""
    assert np.allclose(get_rotation_axis(direction), expected)


@pytest.mark.parametrize("name", sorted(ELECTRODES))
def test_topology_independent_of_direction(name):
    """Tilted electrodes have the topology and volume of the upright one.

    Extra faces or edges indicate seams that were split by a boolean
    operation, i.e. sliver faces that Netgen may fail to mesh.
    """
    upright_topology, upright_volumes = _signature(name, (0, 0, 1))
    for direction in DIRECTIONS:
        topology, volumes = _signature(name, direction)
        assert topology == upright_topology, direction
        np.testing.assert_allclose(volumes, upright_volumes, rtol=1e-6)
