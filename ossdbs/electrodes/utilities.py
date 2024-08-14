# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

from typing import Union

import netgen.occ as occ
import numpy as np


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


def get_signed_angle(
    v1: np.ndarray, v2: np.ndarray, vn: np.ndarray
) -> Union[None, float]:
    """Get signed angle between two vectors.

    Parameters
    ----------
    v1: np.ndarray
        First vector
    v2: np.ndarray
        Second vector
    vn: np.ndarray
        Normal vector in common plane

    Notes
    -----
    See more details on StackOverflow:
    https://stackoverflow.com/questions/5188561/signed-angle-between-two-3d-vectors-with-same-origin-within-the-same-plane
    """
    len_v1 = np.linalg.norm(v1)
    len_v2 = np.linalg.norm(v2)
    len_vn = np.linalg.norm(vn)
    # catch zero-length
    if np.isclose(len_v1, 0.0) or np.isclose(len_v2, 0.0):
        return None

    rotation_angle = np.degrees(np.arccos(np.dot(v1 / len_v1, v2 / len_v2)))
    cross = np.cross(v1, v2)
    if np.dot(vn / len_vn, cross) < 0:
        return -rotation_angle
    return rotation_angle
