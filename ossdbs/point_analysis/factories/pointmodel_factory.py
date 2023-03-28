

from ..point_models import PointModel
from ..point_models import Pathway
from ..point_models import Lattice


class PointModelFactory:

    def create(self, parameters: dict) -> PointModel:
        if parameters['Pathway']['Active']:
            file_name = parameters['Pathway']['FileName']
            return Pathway(file_name)

        shape_par = parameters['Lattice']['Shape']
        shape = shape_par['x'], shape_par['y'], shape_par['z']
        center_par = parameters['Lattice']['Center']
        center = center_par['x[mm]'], center_par['y[mm]'], center_par['z[mm]']
        dir_par = parameters['Lattice']['Direction']
        direction = dir_par['x[mm]'], dir_par['y[mm]'], dir_par['z[mm]']
        distance = parameters['Lattice']['PointDistance[mm]']

        return Lattice(shape=shape,
                       center=center,
                       distance=distance,
                       direction=direction)
