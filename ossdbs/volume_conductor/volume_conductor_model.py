
from abc import ABC, abstractmethod
from dataclasses import dataclass
import ngsolve


@dataclass
class Solution:
    potential: ngsolve.GridFunction
    current_density: ngsolve.GridFunction
    conductivity: ngsolve.CoefficientFunction
    frequency: float
    floating_values: dict


class VolumeConductor(ABC):

    @abstractmethod
    def compute_solution(self, frequency: float) -> Solution:
        pass
