from ossdbs.vta_points import VTAPointMatrix
import numpy as np


class VTAPointsFactory:

    @staticmethod
    def create(vta_parameters: dict) -> np.ndarray:
        shape_par = vta_parameters['Shape']
        shape = shape_par['x'], shape_par['y'], shape_par['z']
        center_par = vta_parameters['Center']
        center = center_par['x[mm]'], center_par['y[mm]'], center_par['z[mm]']
        dir_par = vta_parameters['Direction']
        direction = dir_par['x[mm]'], dir_par['y[mm]'], dir_par['z[mm]']
        distance = vta_parameters['PointDistance[mm]']
        return VTAPointMatrix(shape, center, distance, direction).coordinates()
