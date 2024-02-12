# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np


def get_collapsed_VTA(
    field_on_points: np.ndarray,
    implantation_coordinate: np.ndarray,
    lead_direction: np.ndarray,
    lead_diam: float,
):
    """Postprocess probing points to remove the electrode by inward sideways collapse.

    Notes
    -----
    The point coordinates are pulled into the electrode space to 'counteract'
    tissue misplacement due to the implantation.
    This breaks the grid regularity, so all nifti files
    will be created externally.

    Parameters
    ----------
    field_on_points: numpy.ndarray
        Scalar values on the lattice, Nx7
    implantation_coordinate: numpy array
        Center of the first contact
    lead_direction: numpy array
        Direction of the lead (from head to tail)
    lead_diam: float
        Diameter of the electrode that will be compensated by the inward collapse

    Returns
    -------
    field_on_points_collided: numpy.ndarray
        Array of scalar values on adjusted lattice, Nx7

    """
    # get unit vector for direction
    unit_lead_direction = lead_direction / np.linalg.norm(lead_direction)

    # copy the field values (not sure if it makes sense for vector fields)
    field_on_points_collapsed = np.zeros(field_on_points.shape, float)
    field_on_points_collapsed[:, 3:] = field_on_points[:, 3:]

    # find orthonormal vectors from the points to the lead axis
    relative_coords = (
        np.dot((field_on_points[:, :3] - implantation_coordinate), unit_lead_direction)
        * unit_lead_direction[:, np.newaxis]
    )
    points_on_lead = relative_coords.T + implantation_coordinate[:]

    # now shift points into the direction of the lead axis
    norm_term = 1.0 / np.linalg.norm(points_on_lead - field_on_points[:, :3], axis=1)
    unit_clpse_directions = (points_on_lead - field_on_points[:, :3]) * norm_term[
        :, np.newaxis
    ]
    field_on_points_collapsed[:, :3] = field_on_points[
        :, :3
    ] + unit_clpse_directions * (lead_diam / 2.0)

    return field_on_points_collapsed
