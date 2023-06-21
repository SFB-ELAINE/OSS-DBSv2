

from typing import Tuple
import numpy as np
import h5py
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

    def __init__(self,
                 shape: tuple,
                 center: tuple,
                 distance: float,
                 direction: tuple
                 ) -> None:
        self.__distance = abs(distance)
        self.__shape = shape
        self.__center = center
        norm = np.linalg.norm(direction)
        self.__direction = tuple(direction / norm) if norm else (0, 0, 1)
        self.__location = np.full(shape[0] * shape[1] * shape[2], '')

    def coordinates(self) -> np.ndarray:
        """Generates coordinates of points.

        Returns
        -------
        np.ndarray
        """
        m, n, o = self.__shape
        x_values = (np.arange(m) - ((m - 1) / 2)) * self.__distance
        y_values = (np.arange(n) - ((n - 1) / 2)) * self.__distance
        z_values = (np.arange(o) - ((o - 1) / 2)) * self.__distance

        alpha, beta = self.__rotation_angles_xz()
        coordinates = [self.__rotation((x, y, z), alpha, beta)
                       for x in x_values for y in y_values for z in z_values]

        return np.array(coordinates) + self.__center

    def __rotation(self, point, alpha, beta) -> np.ndarray:
        cos_a = np.cos(alpha)
        sin_a = np.sin(alpha)
        r_x = np.array([[1, 0, 0], [0, cos_a, -sin_a], [0, sin_a, cos_a]])

        cos_b = np.cos(beta)
        sin_b = np.sin(beta)
        r_z = np.array([[cos_b, -sin_b, 0], [sin_b, cos_b, 0], [0, 0, 1]])

        return np.dot(r_z, np.dot(r_x, point))

    def __rotation_angles_xz(self) -> Tuple[float]:
        x_d, y_d, z_d = self.__direction

        if not x_d and not y_d:
            return 0.0, 0.0
        if not y_d:
            return -np.pi / 2, -np.arctan(z_d / x_d)
        if not x_d:
            return 0.0, -np.arctan(z_d / y_d)

        return -np.arctan(y_d / x_d), -np.arctan(z_d / y_d)

    def save(self, data: TimeResult, file_name: str) -> None:
        with h5py.File(file_name, "w") as file:
            self.__write_file(data, file)

    def set_location_names(self, names: np.ndarray) -> None:
        self.__location = names

    def __write_file(self, data, file):
        file.create_dataset('TimeSteps[s]', data=data.time_steps)
        file.create_dataset('Points[mm]', data=data.points)
        file.create_dataset('Location', data=self.__location.astype('S'))
        file.create_dataset('Potential[V]', data=data.potential)
        file.create_dataset('Current_density[A|m2]', data=data.current_density)
