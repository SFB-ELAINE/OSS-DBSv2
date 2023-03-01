

from ossdbs.conductivity import Conductivity
from ossdbs.electrode_collection import Electrodes
from ossdbs.mesh import Mesh
from ossdbs.solver import Solver
from ossdbs.volume_conductor.floating import VolumeConductorFloating
from ossdbs.volume_conductor.floating_impedance import VolumeConductorFloatingImpedance
from ossdbs.volume_conductor.nonfloating import VolumeConductorNonFloating
from ossdbs.volume_conductor.volume_conductor_model import VolumeConductor


class VolumeConductorFactory:

    def __init__(self, mesh: Mesh,
                 conductivity: Conductivity,
                 electrodes: Electrodes,
                 solver: Solver) -> None:
        self.__conductivity = conductivity
        self.__mesh = mesh
        self.__solver = solver
        self.__electrodes = electrodes

    def create(self, floating_parameters: dict) -> VolumeConductor:

        floating = floating_parameters['Active']
        floating_impedance = floating_parameters['FloatingImpedance']
        volume_conductor = self.__select_type(floating, floating_impedance)
        return volume_conductor(conductivity=self.__conductivity,
                                mesh=self.__mesh,
                                electrodes=self.__electrodes,
                                solver=self.__solver)

    @staticmethod
    def __select_type(floating, floating_impedance) -> VolumeConductor:
        if not floating:
            return VolumeConductorNonFloating

        if not floating_impedance:
            return VolumeConductorFloating

        return VolumeConductorFloatingImpedance
