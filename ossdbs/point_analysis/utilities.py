import numpy as np

from ossdbs.fem import VolumeConductor


def label_points(points: np.ndarray, volume_conductor_model: VolumeConductor):
    """Label points by region.

    Notes
    -----
    To be discussed:
    Which format shall be returned

    """
    # points outside domain
    volume_conductor_model.mesh.is_included(points)
    # points in csf
    points_in_csf(points, volume_conductor_model)
    # points in encapsulation layer
    points_in_encapsulation_layer(points, volume_conductor_model)


def points_in_encapsulation_layer(
    points: np.ndarray, volume_conductor_model: VolumeConductor
) -> np.ndarray:
    """Return points in encapsulation layer."""
    # TODO finalize implementation
    # use NGSolve regions to identify points

    pass


def points_in_csf(
    points: np.ndarray, volume_conductor_model: VolumeConductor
) -> np.ndarray:
    """Return points in csf."""
    # use MRI image and find a unique label of all points based on conductivity CF
    pass
