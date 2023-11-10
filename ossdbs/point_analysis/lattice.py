from typing import Tuple

import h5py
import nibabel as nib
import numpy as np

from .point_model import PointModel
from .time_results import TimeResult


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
        self, shape: tuple, center: tuple, distance: float, direction: tuple, collapse_vta: bool
    ) -> None:
        self._distance = abs(distance)
        self._shape = shape
        self._center = center
        norm = np.linalg.norm(direction)
        self._direction = tuple(direction / norm) if norm else (0, 0, 1)
        self._location = np.full(shape[0] * shape[1] * shape[2], "")
        self._coordinates = self._initialize_coordinates()
        self._collapse_vta = collapse_vta

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

    def _rotation_angles_xz(self) -> Tuple[float]:
        x_d, y_d, z_d = self._direction

        if not x_d and not y_d:
            return 0.0, 0.0
        if not y_d:
            return -np.pi / 2, -np.arctan(z_d / x_d)
        if not x_d:
            return 0.0, -np.arctan(z_d / y_d)

        return -np.arctan(y_d / x_d), -np.arctan(z_d / y_d)

    def save(self, data: TimeResult, file_name: str) -> None:
        with h5py.File(file_name, "w") as file:
            self._write_file(data, file)

    def save_as_nifti(
        self, scalar_field, filename, binarize=False, activation_threshold=None
    ):
        """Save scalar field (e.g. electric potential or E-field magnitude) in abstract orthogonal space using nifti
         format.

        Parameters
        ----------
        settings: dict of parameters
        scalar_field : Nx1 numpy.ndarray of scalar values on the lattice
        filename: str, name for the nifti file that should contain full path
        binarize: bool, thresholds the scalar field and saves the binarized result

        """
        # Assuming data is in the same format as it was generated,
        # you can just reshape it
        nifti_grid = scalar_field.reshape(self._shape)

        nifti_output = np.zeros(nifti_grid.shape, float)
        if binarize:
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

        img = nib.Nifti1Image(nifti_output, affine)
        nib.save(img, filename)

    def set_location_names(self, names: np.ndarray) -> None:
        self._location = names

    def _write_file(self, data, file):
        file.create_dataset("TimeSteps[s]", data=data.time_steps)
        file.create_dataset("Points[mm]", data=data.points)
        file.create_dataset("Location", data=self._location.astype("S"))
        file.create_dataset("Potential[V]", data=data.potential)
        file.create_dataset(
            "Electric field magnitude[Vm^(-1)]", data=data.electric_field_magnitude
        )
        file.create_dataset(
            "Electric field vector x[Vm^(-1)]", data=data.electric_field_vector[0]
        )
        file.create_dataset(
            "Electric field vector y[Vm^(-1)]", data=data.electric_field_vector[1]
        )
        file.create_dataset(
            "Electric field vector z[Vm^(-1)]", data=data.electric_field_vector[2]
        )

    def filter_csf_encap(self, inside_csf: np.ndarray, inside_encap: np.ndarray):
        # TODO implement
        raise NotImplementedError("Filtering for lattice not implemented yet.")

    def save_hdf5(
        self,
        axon_mask: list,
        lattice: np.ndarray,
        potentials: np.ndarray,
        fields: np.ndarray,
        field_mags: np.ndarray,
        output_path: str,
    ) -> None:
        raise NotImplementedError("Lattice results can not be stored in HDF5 format.")

    def collapse_VTA(self, field_on_points, implantation_coordinate, lead_direction, lead_diam):
        """Postprocess probing points to remove the electrode by inward sideways collapse. The point coordinates are pulled
        into the electrode space to 'counteract' tissue misplacement due to the implantation. This breaks the grid regularity, so
        all nifti files will be created externally.

        Parameters
        ----------
        field_on_points : Nx4 numpy.ndarray of scalar values on the lattice
        implantation_coordinate: 1-D numpy array, center of the first contact
        lead_direction: 1-D numpy array, direction of the lead (from head to tail)
        lead_diam: float, diameter of the electrode that will be compensated by the inward collapse

        Returns
        -------
        field_on_points_collided: Nx4 numpy.ndarray of scalar values on adjusted lattice

        """
        # get unit vector for direction
        unit_lead_direction = lead_direction / np.linalg.norm(lead_direction)

        # copy the field values (not sure if it makes sense for vector fields)
        field_on_points_collapsed = np.zeros(field_on_points.shape, float)
        field_on_points_collapsed[:, 3:] = field_on_points[:, 3:]

        for point_inx in range(field_on_points.shape[0]):
            # find the orthonormal vector from the point to the lead axis
            point_on_lead = implantation_coordinate + np.dot((field_on_points[point_inx, :3] - implantation_coordinate),
                                                            unit_lead_direction) * unit_lead_direction

            # now shift the point into the direction of the lead axis
            unit_clpse_direction = (point_on_lead - field_on_points[point_inx, :3]) / np.linalg.norm(
                point_on_lead - field_on_points[point_inx, :3])
            field_on_points_collapsed[point_inx, :3] = field_on_points[point_inx, :3] + unit_clpse_direction * (lead_diam / 2.0)

        return field_on_points_collapsed
