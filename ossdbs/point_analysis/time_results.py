from dataclasses import dataclass
import numpy as np
from ossdbs.point_analysis.field_solution import FieldSolution


@dataclass
class TimeResult:
    points: np.ndarray
    time_steps: np.ndarray
    potential: np.ndarray
    current_density: np.ndarray
    field_solution: FieldSolution
