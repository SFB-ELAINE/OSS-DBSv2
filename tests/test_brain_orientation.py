"""Sphere and ellipsoid brains keep their poles and seam away from the leads."""

import numpy as np
import pytest
from netgen.occ import Box, Pnt

import ossdbs
from ossdbs.electrodes import ELECTRODES

TIP = np.array([0.794, -32.671, -0.618])
# sheared and rotated voxel-to-world map, as from an MRI header
TRAFO = np.array([[0.9, 0.2, 0.0], [-0.1, 1.1, 0.1], [0.05, 0.0, 0.8]])


def _lead_direction(tilt_degrees):
    tilt = np.radians(tilt_degrees)
    return np.array([0.0, -np.sin(tilt), np.cos(tilt)])


def _centred_brain(shape, direction, size=40.0):
    """Brain region centred 10 mm above the tip, as set up by Lead-DBS."""
    center = TIP + 10 * direction
    region = ossdbs.create_bounding_box(
        {
            "Center": dict(zip(("x[mm]", "y[mm]", "z[mm]"), center, strict=True)),
            "Dimension": dict.fromkeys(("x[mm]", "y[mm]", "z[mm]"), size),
        }
    )
    return ossdbs.BrainGeometry(shape, region)


def _unit_sphere_map(brain):
    """Map world points of a Sphere/Ellipsoid brain onto its unit sphere."""
    half = np.subtract(brain._bbox.end, brain._bbox.start) / 2
    center = np.add(brain._bbox.start, half)
    if brain._shape == "Sphere":
        half = np.full(3, np.min(half))
    to_world = brain._trafo_matrix @ np.diag(half)
    world_center = brain._trafo_matrix @ center + brain._translation
    to_unit = np.linalg.inv(to_world)
    return lambda point: to_unit @ (np.asarray(point) - world_center)


def _check_orientation(brain, direction, tip=TIP):
    to_unit = _unit_sphere_map(brain)
    unit_tip = to_unit(tip)
    axis = to_unit(tip + direction) - unit_tip
    axis /= np.linalg.norm(axis)
    b = unit_tip @ axis
    exit_point = unit_tip + (np.sqrt(b * b - unit_tip @ unit_tip + 1) - b) * axis
    geometry = brain.geometry
    poles = [to_unit(tuple(vertex.p)) for vertex in geometry.vertices]
    np.testing.assert_allclose([pole @ exit_point for pole in poles], 0, atol=1e-6)
    # the seam is listed once per orientation; the centre of mass of the half
    # meridian opposite the exit point lies at -2/pi * exit_point
    seams = [edge for edge in geometry.edges if edge.mass > 1e-6]
    assert seams
    for seam in seams:
        assert to_unit(tuple(seam.center)) @ exit_point < -0.6


@pytest.mark.parametrize("shape", ["Sphere", "Ellipsoid"])
@pytest.mark.parametrize("rotate_initial_geo", [False, True])
def test_poles_perpendicular_to_lead(shape, rotate_initial_geo):
    direction = _lead_direction(1.9)
    brain = _centred_brain(shape, direction)
    brain._trafo_matrix = TRAFO
    brain._affine_trafo = ossdbs.model_geometry.brain_geometry.netgen.occ.gp_GTrsf(
        mat=TRAFO.ravel(), vec=(0, 0, 0)
    )
    brain._rotate_initial_geo = rotate_initial_geo
    volume = brain.geometry.mass
    brain.set_lead_axes([(TIP, direction)])
    _check_orientation(brain, direction)
    np.testing.assert_allclose(brain.geometry.mass, volume, rtol=1e-9)


def test_custom_geometry_is_kept():
    brain = _centred_brain("Sphere", _lead_direction(0))
    brain.set_geometry(Box(Pnt(-10, -10, -10), Pnt(10, 10, 10)))
    custom = brain.geometry
    brain.set_lead_axes([(TIP, _lead_direction(0))])
    assert brain.geometry is custom


@pytest.mark.parametrize("tilt", [1.82, 1.9])
def test_near_vertical_lead(tilt):
    """Lead leaving the brain next to a pole of the global sphere frame.

    Without orientation, 1.82 degrees leaves sub-micrometre edges that Netgen
    cannot mesh and at 1.90 degrees the contact faces are lost.
    """
    direction = _lead_direction(tilt)
    electrode = ELECTRODES["Medtronic3387"](direction=direction, position=tuple(TIP))
    electrode.encapsulation_thickness = 0.1
    brain = _centred_brain("Ellipsoid", direction)
    geometry = ossdbs.ModelGeometry(brain, [electrode])
    micro_edges = [e.mass for e in geometry._shape.edges if 1e-10 < e.mass < 1e-3]
    assert micro_edges == []
    _check_orientation(brain, direction)
