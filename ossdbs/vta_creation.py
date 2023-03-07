import numpy as np


class VTAPointMatrix:

    def __init__(self, shape: tuple, center: tuple, distance: float) -> None:
        self.__distance = abs(distance)
        self.__shape = shape
        self.__center = center

    def coordinates(self) -> np.ndarray:
        x_c, y_c, z_c = self.__center
        m, n, o = self.__shape
        x_values = x_c + (np.arange(m) - ((m - 1) / 2)) * self.__distance
        y_values = y_c + (np.arange(n) - ((n - 1) / 2)) * self.__distance
        z_values = z_c + (np.arange(o) - ((o - 1) / 2)) * self.__distance

        return np.array([(x, y, z)
                         for x in x_values
                         for y in y_values
                         for z in z_values])
