import numpy as np


def get_collapsed_VTA(field_on_points, implantation_coordinate, lead_direction, lead_diam):
    """Postprocess probing points to remove the electrode by inward sideways collapse. The point coordinates are pulled
    into the electrode space to 'counteract' tissue misplacement due to the implantation. This breaks the grid regularity, so
    all nifti files will be created externally.

    Parameters
    ----------
    field_on_points : Nx7 numpy.ndarray of scalar values on the lattice
    implantation_coordinate: 1-D numpy array, center of the first contact
    lead_direction: 1-D numpy array, direction of the lead (from head to tail)
    lead_diam: float, diameter of the electrode that will be compensated by the inward collapse

    Returns
    -------
    field_on_points_collided: Nx7 numpy.ndarray of scalar values on adjusted lattice

    """
    # get unit vector for direction
    unit_lead_direction = lead_direction / np.linalg.norm(lead_direction)

    # copy the field values (not sure if it makes sense for vector fields)
    field_on_points_collapsed = np.zeros(field_on_points.shape, float)
    field_on_points_collapsed[:, 3:] = field_on_points[:, 3:]

    for point_inx in range(field_on_points.shape[0]):
        # find the orthonormal vector from the point to the lead axis
        point_on_lead = implantation_coordinate + np.dot((field_on_points[point_inx, :3] - implantation_coordinate), unit_lead_direction) * unit_lead_direction

        # now shift the point into the direction of the lead axis
        unit_clpse_direction = (point_on_lead - field_on_points[point_inx, :3]) / np.linalg.norm(
            point_on_lead - field_on_points[point_inx, :3])
        field_on_points_collapsed[point_inx, :3] = field_on_points[point_inx, :3] + unit_clpse_direction * (lead_diam / 2.0)

    return field_on_points_collapsed
