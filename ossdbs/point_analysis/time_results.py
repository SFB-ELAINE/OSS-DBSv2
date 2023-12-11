from dataclasses import dataclass

import numpy as np


@dataclass
class TimeResult:
    """TODO format of electric_field_vector needs to be clarified."""

    points: np.ndarray
    time_steps: np.ndarray
    potential: np.ndarray
    electric_field_magnitude: np.ndarray
    electric_field_vector: np.ndarray
    inside_csf: np.ndarray
    inside_encap: np.ndarray
