

from ossdbs.point_models import PointModel
from ossdbs.point_models import PathwayActivationModelling
from ossdbs.point_models import Lattice


class ActivationModelFactory:

    def create(self, parameters: dict) -> PointModel:
        if parameters['PathwayActivationModelling']['Active']:
            file_name = parameters['PathwayActivationModelling']['FileName']
            return PathwayActivationModelling(file_name)

        shape_par = parameters['VTA']['Shape']
        shape = shape_par['x'], shape_par['y'], shape_par['z']
        center_par = parameters['VTA']['Center']
        center = center_par['x[mm]'], center_par['y[mm]'], center_par['z[mm]']
        dir_par = parameters['VTA']['Direction']
        direction = dir_par['x[mm]'], dir_par['y[mm]'], dir_par['z[mm]']
        distance = parameters['VTA']['PointDistance[mm]']

        return Lattice(shape=shape,
                                              center=center,
                                              distance=distance,
                                              direction=direction)
