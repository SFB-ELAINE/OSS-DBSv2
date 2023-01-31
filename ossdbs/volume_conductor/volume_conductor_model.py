
from abc import ABC, abstractmethod
from dataclasses import dataclass
from ossdbs.electrode_contacts import ContactCollection
from ossdbs.conductivity import Conductivity
from ossdbs.mesh import Mesh
from ossdbs.solver import Solver
import ngsolve


@dataclass
class Potential:
    gridfunction: ngsolve.GridFunction
    frequency: float
    floating_values: dict


class VolumeConductor(ABC):

    @abstractmethod
    def potential(self, frequency: float) -> Potential:
        pass
