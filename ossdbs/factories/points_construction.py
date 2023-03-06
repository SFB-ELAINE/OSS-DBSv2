import h5py
import numpy as np


class PointsFactory:

    @staticmethod
    def create(path: str) -> np.ndarray:
        with h5py.File(path, "r") as file:
            points = [np.array(file[key]) for key in file.keys()]
        return np.concatenate(points, axis=0)

    @staticmethod
    def categories(path: str) -> list:
        with h5py.File(path, "r") as file:
            categories = [(key, file[key].shape[0]) for key in file.keys()]
        return categories
