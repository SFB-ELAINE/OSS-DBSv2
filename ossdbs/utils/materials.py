# Copyright 2023, 2024 Konstantin Butenko, Jan Philipp Payonk, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later
import numpy as np

MATERIALS = {"Blood": 4, "Gray matter": 1, "White matter": 2, "CSF": 3, "Unknown": 0}


def have_dielectric_properties_changed(
    dielectric_properties: dict, is_complex, old_freq, new_freq, threshold
):
    """Check if dielectric properties have changed compared to last run."""
    max_error = 0.0
    for _material, model in dielectric_properties.items():
        if is_complex:
            old_value = model.conductivity(2.0 * np.pi * old_freq)
            new_value = model.conductivity(2.0 * np.pi * new_freq)
            error_real = np.abs(old_value.real - new_value.real)
            error_real /= np.abs(old_value)
            error_imag = np.abs(old_value.imag - new_value.imag)
            error_imag /= np.abs(old_value)
            error = np.maximum(error_real, error_imag)
        else:
            old_value = model.conductivity(2.0 * np.pi * old_freq)
            new_value = model.conductivity(2.0 * np.pi * new_freq)
            error = np.abs((old_value - new_value) / old_value)
        if error > max_error:
            max_error = error
    if max_error > threshold:
        return True
    return False
