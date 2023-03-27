
from ossdbs.conductivity import Conductivity
from ossdbs.mesh import Mesh
from ossdbs.solver import Solver
from ossdbs.volume_conductor import VolumeConductor
from ossdbs.volume_conductor import VolumeConductorFloating
from ossdbs.volume_conductor import VolumeConductorFloatingImpedance
from ossdbs.volume_conductor import VolumeConductorNonFloating


class VolumeConductorFactory:

    def __init__(self, mesh: Mesh,
                 conductivity: Conductivity,
                 solver: Solver) -> None:
        self.__conductivity = conductivity
        self.__mesh = mesh
        self.__solver = solver

    def create(self, floating_parameters: dict) -> VolumeConductor:

        floating = floating_parameters['Active']
        floating_impedance = floating_parameters['FloatingImpedance']
        volume_conductor = self.__select_type(floating, floating_impedance)
        return volume_conductor(conductivity=self.__conductivity,
                                mesh=self.__mesh,
                                solver=self.__solver)

    @staticmethod
    def __select_type(floating: bool,
                      floating_impedance: bool
                      ) -> VolumeConductor:
        if not floating:
            return VolumeConductorNonFloating

        if not floating_impedance:
            return VolumeConductorFloating

        return VolumeConductorFloatingImpedance
