from dataclasses import dataclass
import numpy as np
from ossdbs.utils.vtk_export import FieldSolution


@dataclass
class TimeResult:
    # TODO check what is really needed
    points: np.ndarray
    time_steps: np.ndarray
    potential: np.ndarray
    current_density: np.ndarray
    field_solution: FieldSolution
