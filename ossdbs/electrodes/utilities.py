# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import logging

import netgen.occ as occ
import numpy as np
from scipy.spatial.transform import Rotation

_logger = logging.getLogger(__name__)


def _edge_height(edge: occ.Edge, direction: tuple) -> float:
    """Position of the edge centre along ``direction``."""
    center = edge.center
    return float(np.dot((center.x, center.y, center.z), direction))


def get_lowest_edge(contact: occ.Face, direction: tuple = (0, 0, 1)) -> occ.Edge:
    """Get lowest edge along ``direction`` (default: z-direction).

    Pass the lead direction for shapes built along a tilted axis; comparing
    global z there picks a cylinder seam line instead of the rim circle.
    """
    return min(contact.edges, key=lambda edge: _edge_height(edge, direction))


def get_highest_edge(contact: occ.Face, direction: tuple = (0, 0, 1)) -> occ.Edge:
    """Get highest edge along ``direction`` (default: z-direction).

    Pass the lead direction for shapes built along a tilted axis; comparing
    global z there picks a cylinder seam line instead of the rim circle.
    """
    return max(contact.edges, key=lambda edge: _edge_height(edge, direction))


def get_rotation_axis(direction: tuple) -> tuple:
    """Axis to rotate the upright z axis onto ``direction``.

    Returns the unit vector of ``z x direction``. When ``direction`` is
    (anti)parallel to z the cross product vanishes; any axis perpendicular to
    z is then valid (the rotation angle is 0 or 180 degrees), so the x axis is
    returned instead of dividing by zero.
    """
    cross = np.cross((0, 0, 1), np.asarray(direction, dtype=float))
    norm = np.linalg.norm(cross)
    if np.isclose(norm, 0.0):
        return (1.0, 0.0, 0.0)
    return tuple(cross / norm)


def rotate_sphere_seam(sphere, center: tuple, direction: tuple):
    """Move a sphere's BREP pole off the electrode axis.

    ``occ.Sphere`` takes no direction, so its two pole vertices always sit at
    global z relative to the centre. When the electrode points along z, those
    poles land on the lead axis, and the resulting degenerate topology can
    make Netgen's surface mesher fail or crash. Rotating about an axis
    perpendicular to both z and ``direction``, by the z-to-direction angle
    plus 45 degrees, leaves the poles at 45 degrees to ``direction`` for any
    direction. A sphere is symmetric about its centre, so the shape itself is
    unchanged -- only the seam moves.

    The offset must not be 90 degrees: that puts the poles on the great circle
    where the tip sphere is tangent to the lead cylinder, and OCC booleans on
    that configuration fail for some directions (e.g. the fused tip + lead
    silently loses the cylinder).

    Call this on the freshly constructed sphere, before combining it with
    anything else.
    """
    direction = np.asarray(direction, dtype=float)
    direction = direction / np.linalg.norm(direction)
    axis = get_rotation_axis(direction)
    angle = np.degrees(np.arccos(np.clip(direction[2], -1.0, 1.0))) + 45.0
    return sphere.Rotate(occ.Axis(p=occ.Pnt(*center), d=occ.Dir(*axis)), angle)


def get_signed_angle(
    v_in: np.ndarray, v_out: np.ndarray, rotation_axis: np.ndarray
) -> None | float:
    """Get signed angle between two vectors.

    Parameters
    ----------
    v_in: np.ndarray
        First vector which needs rotation
    v_out: np.ndarray
        Second vector, which should be matched
    rotation_axis: np.ndarray
        Axis around which the vectors will be rotated
    """
    len_v_in = np.linalg.norm(v_in)
    len_v_out = np.linalg.norm(v_out)
    len_r_axis = np.linalg.norm(rotation_axis)
    # catch zero-length
    if np.isclose(len_v_in, 0.0) or np.isclose(len_v_out, 0.0):
        return None
    if np.isclose(len_r_axis, 0.0):
        raise ValueError("Rotation axis has length zero")

    # determine rotation angle
    rotation_angle = np.degrees(np.arccos(np.dot(v_in / len_v_in, v_out / len_v_out)))
    # apply rotation to input vector
    rotation = Rotation.from_rotvec(
        np.radians(rotation_angle) * rotation_axis / len_r_axis
    )
    corrected_direction = rotation.apply(v_in / len_v_in)
    # if the correction does not match, flip the sign
    if not np.all(np.isclose(corrected_direction, v_out / len_v_out, atol=1e-5)):
        rotation_angle = -rotation_angle
        rotation = Rotation.from_rotvec(
            np.radians(rotation_angle) * rotation_axis / len_r_axis
        )
        corrected_direction = rotation.apply(v_in / len_v_in)

        if not np.all(np.isclose(corrected_direction, v_out / len_v_out, atol=1e-5)):
            raise RuntimeError(
                "Could not determine signed angle between vectors. "
                "Possible reasons: numerical accuracy or wrong geometry "
                "information."
            )
    return rotation_angle


def get_electrode_spin_angle(
    rotation: tuple, angle: float, direction: np.ndarray
) -> float:
    """Determine angle that directed electrode needs to be spinned."""
    # adjust contact angle
    # tilted y-vector marker is in YZ-plane and orthogonal to _direction
    # note that this comes from Lead-DBS
    desired_direction = np.array([0, direction[2], -direction[1]])
    rotate_vector = Rotation.from_rotvec(np.radians(angle) * np.array(rotation))
    current_direction = rotate_vector.apply((0, 1, 0))
    # get angle between current and desired direction
    # current direction is normal
    rotation_angle = get_signed_angle(
        current_direction, desired_direction, np.array(direction)
    )
    if rotation_angle is None:
        _logger.warning(
            "Could not determine rotation angle for "
            "correct spin as per Lead-DBS convention."
            "Returning angle of zero."
        )
        # to return unrotated geo
        rotation_angle = 0.0
    return rotation_angle
