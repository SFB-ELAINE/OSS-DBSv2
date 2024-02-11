# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import logging

import numpy as np

_logger = logging.getLogger(__name__)


def compute_field_magnitude(fields: np.ndarray):
    """Compute magnitude of field vector."""
    if len(fields.shape) != 2 and fields.shape[1] != 3:
        raise ValueError("Attempt to compute field magnitude from invalid field vector")
    # TODO fix for complex vectors
    field_mags = np.linalg.norm(fields, axis=1).reshape((fields.shape[0], 1))
    _logger.debug(f"Shapes of field_mags: {field_mags.shape}")
    return field_mags


def compute_field_magnitude_from_components(
    Ex: np.ndarray, Ey: np.ndarray, Ez: np.ndarray
):
    """Compute magnitude of field vector."""
    Ex2 = (Ex * np.conj(Ex)).real
    Ey2 = (Ey * np.conj(Ey)).real
    Ez2 = (Ez * np.conj(Ez)).real
    field_mags = np.sqrt(Ex2 + Ey2 + Ez2)
    _logger.debug(f"Shapes of field_mags: {field_mags.shape}")
    return field_mags
