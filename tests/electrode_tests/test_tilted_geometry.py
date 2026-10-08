"""Regression tests for electrode geometry along tilted and reversed leads."""

import netgen.occ as occ
import numpy as np
import pytest

from ossdbs.electrodes import ELECTRODES, Medtronic3387, MicroProbesSNEX100
from ossdbs.electrodes.utilities import (
    get_highest_edge,
    get_lowest_edge,
    get_rotation_axis,
)

TILTED_DIRECTIONS = [
    (0.0, 0.5, 0.866),
    (0.6, 0.0, 0.8),
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.302, 0.302, -0.905),
    (0.0, 0.0, -1.0),
]


def _unit(direction):
    direction = np.asarray(direction, dtype=float)
    return tuple(direction / np.linalg.norm(direction))


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
    axis = get_rotation_axis(direction)
    assert np.allclose(axis, expected)


@pytest.mark.parametrize("direction", TILTED_DIRECTIONS)
def test_edge_selection_along_direction(direction):
    """The lowest / highest edge along the lead is a rim circle."""
    direction = _unit(direction)
    start = np.array([1.0, 2.0, 3.0])
    height = 2.0
    cylinder = occ.Cylinder(p=tuple(start), d=direction, r=0.5, h=height)
    for edge, rim_center in (
        (get_lowest_edge(cylinder, direction), start),
        (get_highest_edge(cylinder, direction), start + height * np.array(direction)),
    ):
        center = edge.center
        assert np.allclose((center.x, center.y, center.z), rim_center)


@pytest.mark.parametrize("direction", TILTED_DIRECTIONS)
def test_snex100_tilted_encapsulation(direction):
    """SNEX100 encapsulation fillets the rims, so it builds for any tilt."""
    electrode = MicroProbesSNEX100(direction=_unit(direction))
    encapsulation = electrode.encapsulation_geometry(0.1)
    upright = MicroProbesSNEX100().encapsulation_geometry(0.1)
    assert np.isclose(encapsulation.mass, upright.mass, rtol=1e-6)


def test_medtronic3387_oblique_encapsulation():
    """Tip sphere + lead cylinder fuse keeps the cylinder (PR 142 report)."""
    direction = (-0.4565710794209321, 0.2172862364392653, 0.8627453511265445)
    thickness = 0.1
    electrode = Medtronic3387(direction=direction, position=(0.794, -32.671, -0.618))
    encapsulation = electrode.encapsulation_geometry(thickness)
    parameters = electrode._parameters
    radius = parameters.lead_diameter * 0.5
    outer = radius + thickness
    height = parameters.total_length - parameters.tip_length
    expected = np.pi * (outer**2 - radius**2) * height + 2.0 / 3.0 * np.pi * (
        outer**3 - radius**3
    )
    assert np.isclose(encapsulation.mass, expected, rtol=1e-6)


@pytest.mark.parametrize("name", sorted(ELECTRODES))
def test_reversed_direction(name):
    """Electrodes pointing along -z build with the same volume as along +z."""
    upright = ELECTRODES[name](direction=(0, 0, 1))
    reversed_lead = ELECTRODES[name](direction=(0, 0, -1))
    assert np.isclose(reversed_lead.geometry.mass, upright.geometry.mass, rtol=1e-6)
