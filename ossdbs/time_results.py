from dataclasses import dataclass
import numpy as np


@dataclass
class TimeResult:
    points: np.ndarray
    time_steps: np.ndarray
    potential: np.ndarray
    current_density: np.ndarray
