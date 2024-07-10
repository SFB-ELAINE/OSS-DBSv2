# Copyright 2023, 2024 Konstantin Butenko, Jan Philipp Payonk, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

"""Open-source software for deep brain stimulation."""
import logging

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
from ossdbs.main import main_run, set_logger
from ossdbs.model_geometry import BrainGeometry, ModelGeometry
from ossdbs.utils.nifti1image import DiffusionTensorImage, MagneticResonanceImage

_logger = logging.getLogger(__name__)
_logger.addHandler(logging.NullHandler())


__all__ = (
    "set_logger",
    "BrainGeometry",
    "ModelGeometry",
    "ConductivityCF",
    "generate_electrodes",
    "generate_brain_model",
    "generate_model_geometry",
    "generate_mesh",
    "prepare_solver",
    "prepare_volume_conductor_model",
    "prepare_dielectric_properties",
    "prepare_stimulation_signal",
    "load_images",
    "create_bounding_box",
    "MagneticResonanceImage",
    "DiffusionTensorImage",
    "set_contact_and_encapsulation_layer_properties",
    "Mesh",
    "run_volume_conductor_model",
    "generate_point_models",
    "main_run",
)
