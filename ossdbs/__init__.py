# Copyright 2023, 2024 Konstantin Butenko, Jan Philipp Payonk, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

"""Open-source software for deep brain stimulation."""

import logging

import ngsolve

from ossdbs.api import (
    create_bounding_box,
    generate_brain_model,
    generate_electrodes,
    generate_mesh,
    generate_model_geometry,
    generate_point_models,
    load_images,
    prepare_dielectric_properties,
    prepare_solver,
    prepare_stimulation_signal,
    prepare_volume_conductor_model,
    run_volume_conductor_model,
    set_contact_and_encapsulation_layer_properties,
)
from ossdbs.fem import ConductivityCF, Mesh
from ossdbs.model_geometry import BrainGeometry, ModelGeometry
from ossdbs.utils.nifti1image import (
    DiffusionTensorImage,
    MagneticResonanceImage,
    VTAImage,
)

_logger = logging.getLogger(__name__)
_logger.addHandler(logging.NullHandler())


def log_to_file(output_file: str, level=logging.INFO):
    """Write logging output also to file."""
    # overwrite the previous log
    fh = logging.FileHandler(output_file, mode="w")
    fh.setLevel(level)
    fh.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
    _logger.addHandler(fh)


def set_logger(level=logging.INFO):
    """Set log level."""
    _logger.setLevel(level)
    if level == logging.DEBUG:
        ngsolve.ngsglobals.msg_level = 10
    # to avoid multiple output in Jupyter notebooks
    if len(_logger.handlers) == 1:
        ch = logging.StreamHandler()
        ch.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
        ch.setLevel(level)
        _logger.addHandler(ch)
    else:
        for handler in _logger.handlers:
            if isinstance(handler, logging.StreamHandler):
                handler.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
                handler.setLevel(level)


__all__ = (
    "BrainGeometry",
    "ConductivityCF",
    "DiffusionTensorImage",
    "MagneticResonanceImage",
    "Mesh",
    "ModelGeometry",
    "VTAImage",
    "create_bounding_box",
    "generate_brain_model",
    "generate_electrodes",
    "generate_mesh",
    "generate_model_geometry",
    "generate_point_models",
    "load_images",
    "prepare_dielectric_properties",
    "prepare_solver",
    "prepare_stimulation_signal",
    "prepare_volume_conductor_model",
    "run_volume_conductor_model",
    "set_contact_and_encapsulation_layer_properties",
    "set_logger",
)
