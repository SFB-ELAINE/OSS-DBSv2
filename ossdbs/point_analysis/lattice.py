# Copyright 2023, 2024 Konstantin Butenko, Jan Philipp Payonk, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import json
import logging
from typing import Optional

import nibabel
import numpy as np

from .point_model import PointModel

_logger = logging.getLogger(__name__)


class Lattice(PointModel):
    """Matrix of point coordinates.

    Attributes
    ----------
    shape : tuple
        Number of points in each direction (x, y, z).
    center : tuple
        Center position of cuboid matrix.
    distance : float
        Distance between adjacent points.
    direction : tuple
        Orientation of cuboid in 3d space.
    """

    def __init__(
        self,
        shape: tuple,
        center: tuple,
        distance: float,
        direction: tuple,
        collapse_vta: bool = False,
        export_field: bool = True,
    ) -> None:
        if distance < 0:
            raise ValueError("The spacing between points must be positive.")
        if len(shape) != 3:
            raise ValueError("Pass a 3-valued tuple as the lattice shape.")
        self._distance = distance
        self._shape = shape
        self._collapse_VTA = collapse_vta
        self._export_field = export_field
        self._center = center
        norm = np.linalg.norm(direction)
        # TODO why can norm be not be there?
        self._direction = tuple(direction / norm) if norm else (0, 0, 1)
        self._location = np.full(shape[0] * shape[1] * shape[2], "")
        self._coordinates = self._initialize_coordinates()

        # identifiers
        self._name = "Lattice"

        # never compute time-domain signal
        self._time_domain_conversion = False

        # VTA volume
        self._vta_volume = None

    @property
    def VTA_volume(self) -> Optional[float]:
        """Return VTA volume in mm^3."""
        return self._vta_volume

    @VTA_volume.setter
    def VTA_volume(self, value: float) -> None:
        """Set VTA volume in mm^3."""
        self._vta_volume = value

    def _initialize_coordinates(self) -> np.ndarray:
        """Generates coordinates of points.

        Returns
        -------
        np.ndarray
        """
        m, n, o = self._shape
        x_values = (np.arange(m) - ((m - 1) / 2)) * self._distance
        y_values = (np.arange(n) - ((n - 1) / 2)) * self._distance
        z_values = (np.arange(o) - ((o - 1) / 2)) * self._distance

        alpha, beta = self._rotation_angles_xz()
        coordinates = [
            self._rotation((x, y, z), alpha, beta)
            for x in x_values
            for y in y_values
            for z in z_values
        ]

        return np.array(coordinates) + self._center

    def _rotation(self, point, alpha, beta) -> np.ndarray:
        cos_a = np.cos(alpha)
        sin_a = np.sin(alpha)
        r_x = np.array([[1, 0, 0], [0, cos_a, -sin_a], [0, sin_a, cos_a]])

        cos_b = np.cos(beta)
        sin_b = np.sin(beta)
        r_z = np.array([[cos_b, -sin_b, 0], [sin_b, cos_b, 0], [0, 0, 1]])

        return np.dot(r_z, np.dot(r_x, point))

    def _rotation_angles_xz(self) -> tuple[float]:
        x_d, y_d, z_d = self._direction

        if not x_d and not y_d:
            return 0.0, 0.0
        if not y_d:
            return -np.pi / 2, -np.arctan(z_d / x_d)
        if not x_d:
            return 0.0, -np.arctan(z_d / y_d)

        return -np.arctan(y_d / x_d), -np.arctan(z_d / y_d)

    def save_as_nifti(
        self,
        scalar_field: np.ndarray,
        filename: str,
        binarize: bool = False,
        activation_threshold: Optional[float] = None,
    ):
        """Save scalar field in abstract orthogonal space in nifti format.

        Parameters
        ----------
        scalar_field : numpy.ndarray
            Nx1 array of scalar values on the lattice
        filename: str
            Name for the nifti file that should contain full path
        binarize: bool
            Choose to threshold the scalar field and save the binarized result
        activation_threshold: float
            Activation threshold for VTA estimate
        """
        # Assuming data is in the same format as it was generated,
        # you can just reshape it
        nifti_grid = scalar_field.reshape(self._shape)

        nifti_output = np.zeros(nifti_grid.shape, float)
        if binarize:
            if activation_threshold is None:
                raise ValueError("Provide an activation threshold.")
            nifti_output[nifti_grid >= activation_threshold] = 1
            nifti_output[nifti_grid < activation_threshold] = 0
        else:
            nifti_output = nifti_grid  # V/mm

        # create an abstract nifti
        # define affine transform with the correct resolution and offset
        affine = np.eye(4)
        affine[0:3, 3] = [
            self.coordinates[0][0],
            self.coordinates[0][1],
            self.coordinates[0][2],
        ]
        affine[0, 0] = self._distance
        affine[1, 1] = self._distance
        affine[2, 2] = self._distance

        nibabel.save(nibabel.Nifti1Image(nifti_output, affine), filename)

    def export_point_model_information(self, filename: str) -> None:
        """Export all relevant information about the model to JSON."""
        if not filename.endswith(".json"):
            _logger.warning(
                "Filename for export did not end with `json`, "
                "added `json` as fileending."
            )
            filename += ".json"
        vta_info = {
            "distance": self._distance,
            "shape": self._shape,
            "collapse_VTA": self.collapse_VTA,
            "export_field": self._export_field,
            "center": self._center,
            "direction": self._direction,
            "volume": self.VTA_volume,
        }
        # write to file
        with open(filename, "w") as fp:
            json.dump(vta_info, fp)
