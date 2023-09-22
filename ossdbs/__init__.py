"""Open-source software for deep brain stimulation."""
import logging

import ngsolve

from ossdbs.api import (
    create_bounding_box,
    generate_brain_model,
    generate_electrodes,
    generate_mesh,
    generate_model_geometry,
    generate_neuron_grid,
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
from ossdbs.utils.nifti1image import DiffusionTensorImage, MagneticResonanceImage

_logger = logging.getLogger(__name__)
_logger.addHandler(logging.NullHandler())


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
            if type(handler) == logging.StreamHandler:
                handler.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
                handler.setLevel(level)


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
    "generate_neuron_grid",
)
