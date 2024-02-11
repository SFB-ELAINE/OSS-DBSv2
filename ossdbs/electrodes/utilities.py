# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import netgen.occ as occ


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
