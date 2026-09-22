# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import logging

import netgen.occ as occ
import numpy as np
from scipy.spatial.transform import Rotation

_logger = logging.getLogger(__name__)


def get_lowest_edge(contact: occ.Face) -> occ.Edge:
    """Get lowest edge (i.e., int z-direction)."""
    min_edge_val = float("inf")
    for edge in contact.edges:
        if edge.center.z < min_edge_val:
            min_edge_val = edge.center.z
            min_edge = edge
    return min_edge


def get_highest_edge(contact: occ.Face) -> occ.Edge:
    """Get highest edge (i.e., in z-direction)."""
    max_edge_val = float("-inf")
    for edge in contact.edges:
        if edge.center.z > max_edge_val:
            max_edge_val = edge.center.z
            max_edge = edge
    return max_edge


def rotate_sphere_seam(sphere, center: tuple, direction: tuple):
    """Move a sphere's BREP pole off the electrode axis.

    ``occ.Sphere`` takes no direction, so its two pole vertices always sit at
    global z relative to the centre. When the electrode points along z, those
    poles land on the lead axis, and the resulting degenerate topology can
    make Netgen's surface mesher fail or crash. Rotating about an axis
    perpendicular to both z and ``direction``, by the z-to-direction angle
    plus 90 degrees, leaves the poles perpendicular to ``direction`` for any
    direction. A sphere is symmetric about its centre, so the shape itself is
    unchanged -- only the seam moves.

    Call this on the freshly constructed sphere, before combining it with
    anything else.
    """
    direction = np.asarray(direction, dtype=float)
    direction = direction / np.linalg.norm(direction)
    cross = np.cross((0, 0, 1), direction)
    norm = np.linalg.norm(cross)
    # direction is (anti)parallel to z: any axis perpendicular to z will do
    axis = (1.0, 0.0, 0.0) if np.isclose(norm, 0.0) else tuple(cross / norm)
    angle = np.degrees(np.arccos(np.clip(direction[2], -1.0, 1.0))) + 90.0
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
