# Copyright 2023, 2024 Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

from dataclasses import dataclass

import numpy as np


@dataclass
class TimeResult:
    """Compile time-domain solution.

    Notes
    -----
    The electric field is not always required
    and is not added by default.
    """

    time_steps: np.ndarray
    points: np.ndarray
    inside_csf: np.ndarray
    inside_encap: np.ndarray
    potential: np.ndarray
    electric_field_magnitude: np.ndarray = None
    electric_field_vector_x: np.ndarray = None
    electric_field_vector_y: np.ndarray = None
    electric_field_vector_z: np.ndarray = None
