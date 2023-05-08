
from abc import ABC, abstractmethod
from dataclasses import dataclass
import ngsolve
from ossdbs.conductivity import Conductivity
from ossdbs.fem.mesh import Mesh
from ossdbs.fem.solver import Solver
from ossdbs.model_geometry import ModelGeometry


@dataclass
class Solution:
    potential: ngsolve.GridFunction
    current_density: ngsolve.GridFunction
    conductivity: ngsolve.CoefficientFunction
    frequency: float
    floating_values: dict


class VolumeConductor(ABC):
    """Template class of a volume conductor

    Parameters
    ----------
    mesh : Mesh
    conductivity : Conductivity
    model_geometry : ModelGeometry
    solver : Solver

    # TODO add more abstractmethod ?
    """

    def __init__(self,
                 mesh: Mesh,
                 conductivity: Conductivity,
                 solver: Solver,
                 model_geometry: ModelGeometry) -> None:
        self.conductivity = conductivity
        self.mesh = mesh
        self.solver = solver
        self.model_geometry = model_geometry

    @abstractmethod
    def compute_solution(self,
                         frequency: float) -> Solution:
        pass
