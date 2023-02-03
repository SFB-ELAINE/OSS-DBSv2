
from dataclasses import dataclass
import numpy as np


@dataclass
class Region:
    start: tuple = (0., 0., 0.)
    end: tuple = (0., 0., 0.)
    shape: tuple = (1, 1, 1)
    step_size: tuple = (1., 1., 1.)

    def coordinates(self) -> np.ndarray:
        x_start, y_start, z_start = self.start
        x_steps, y_steps, z_steps = self.shape
        x_dim, y_dim, z_dim = self.step_size

        x_pos = x_start + np.arange(x_steps) * x_dim
        y_pos = y_start + np.arange(y_steps) * y_dim
        z_pos = z_start + np.arange(z_steps) * z_dim

        return np.array([[x, y, z]
                         for x in x_pos for y in y_pos for z in z_pos])
