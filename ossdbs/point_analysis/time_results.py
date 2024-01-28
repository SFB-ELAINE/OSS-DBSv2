from dataclasses import dataclass

import numpy as np


@dataclass
class TimeResult:
    """Compile time-domain solution."""

    points: np.ndarray
    time_steps: np.ndarray
    potential: np.ndarray
    electric_field_magnitude: np.ndarray
    electric_field_vector_x: np.ndarray
    electric_field_vector_y: np.ndarray
    electric_field_vector_z: np.ndarray
    inside_csf: np.ndarray
    inside_encap: np.ndarray
