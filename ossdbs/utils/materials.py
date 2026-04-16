# Copyright 2023, 2024 Konstantin Butenko, Jan Philipp Payonk, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

"""Material-related utilities."""

import numpy as np

MATERIALS = {"Blood": 4, "Gray matter": 1, "White matter": 2, "CSF": 3, "Unknown": 0}


def have_dielectric_properties_changed(
    dielectric_properties: dict,
    is_complex: bool,
    old_freq: float,
    new_freq: float,
    threshold: float,
) -> bool:
    """Check if dielectric properties have changed significantly between frequencies.

    This function compares the dielectric properties at two frequencies
    and determines whether the change exceeds a threshold. This is useful
    for optimizing FEM solves by reusing the mesh when properties haven't
    changed significantly.

    Parameters
    ----------
    dielectric_properties : dict
        Dictionary mapping material names to dielectric models.
        Each model should have a `conductivity(omega)` method.
    is_complex : bool
        Whether to compare complex conductivity values.
    old_freq : float
        Previous frequency in Hz.
    new_freq : float
        New frequency in Hz.
    threshold : float
        Relative error threshold (e.g., 0.01 for 1%).

    Returns
    -------
    bool
        True if properties have changed more than threshold, False otherwise.

    Example
    -------
    >>> changed = have_dielectric_properties_changed(
    ...     dielectric_properties={'Gray matter': gray_model},
    ...     is_complex=True,
    ...     old_freq=100.0,
    ...     new_freq=1000.0,
    ...     threshold=0.01,
    ... )
    """
    max_error = 0.0

    for _material, model in dielectric_properties.items():
        old_omega = 2.0 * np.pi * old_freq
        new_omega = 2.0 * np.pi * new_freq

        if is_complex:
            old_value = model.complex_conductivity(old_omega)
            new_value = model.complex_conductivity(new_omega)
            # Compute relative errors for real and imaginary parts
            error_real = np.abs(old_value.real - new_value.real)
            error_real /= np.abs(old_value)
            error_imag = np.abs(old_value.imag - new_value.imag)
            error_imag /= np.abs(old_value)
            error = np.maximum(error_real, error_imag)
        else:
            old_value = model.conductivity(old_omega)
            new_value = model.conductivity(new_omega)
            error = np.abs((old_value - new_value) / old_value)

        if error > max_error:
            max_error = error

    return max_error > threshold
