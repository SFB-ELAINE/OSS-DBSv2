
from attr import dataclass
import ngsolve
import os
import pickle
import numpy as np


@dataclass
class FrequencyComponent:
    fourier_coefficient: float
    frequency: float


class Output:

    def __init__(self, mesh, potential) -> None:
        self.__mesh = mesh
        self.__potential = potential

    def save_mesh(self, path: str = '') -> None:
        file_base_name = os.path.basename(path)
        file_dir = os.path.dirname(path)

        if not file_base_name:
            file_base_name = 'mesh.vol'

        if not file_dir:
            file_dir = 'result'

        if not os.path.exists(file_dir):
            os.makedirs(file_dir)

        filename = os.path.join(file_dir, file_base_name)
        self.__mesh.ngmesh.Save(filename)

    def save(self, path: str = '') -> None:

        file_base_name = os.path.basename(path)
        file_dir = os.path.dirname(path)

        if not file_base_name:
            file_base_name = 'result'

        if not file_dir:
            file_dir = 'result'

        if not os.path.exists(file_dir):
            os.makedirs(file_dir)

        filename = os.path.join(file_dir, 'potential.data')
        pickler = pickle.Pickler(open(filename, 'wb'))
        pickler.dump([self.__potential, self.__mesh])

        filename = os.path.join(file_dir, file_base_name)
        ngsolve.VTKOutput(ma=self.__mesh,
                          coefs=[self.__potential.real, self.__potential.imag],
                          names=["potential_real", "potential_imag"],
                          filename=filename,
                          subdivision=0
                          ).Do()
