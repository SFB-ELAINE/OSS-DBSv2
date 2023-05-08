import logging
from .__version__ import __version__
from .nifti1Image import Nifti1Image

from ossdbs.factories import (ElectrodeFactory,
                              BrainRegionFactory)
from ossdbs.model_geometry import (BrainGeometry,
                                   ModelGeometry)

_logger = logging.getLogger(__name__)
_logger.addHandler(logging.NullHandler())


def generate_electrodes(settings):
    """Generate OCC electrode models

    Notes
    -----
    TODO type checking

    """
    electrodes = []
    for electrode_parameters in settings["Electrodes"]:
        name = electrode_parameters["Name"]
        direction = (electrode_parameters["Direction"]["x[mm]"],
                     electrode_parameters["Direction"]["y[mm]"],
                     electrode_parameters["Direction"]["z[mm]"])
        rotation = electrode_parameters["Rotation[Degrees]"]
        position = (electrode_parameters["TipPosition"]["x[mm]"],
                    electrode_parameters["TipPosition"]["y[mm]"],
                    electrode_parameters["TipPosition"]["z[mm]"])

        electrode = ElectrodeFactory.create(name, direction, position, rotation)
        if "EncapsulationLayer" in electrode_parameters:
            electrode.encapsulation_thickness = electrode_parameters["EncapsulationLayer"]["Thickness[mm]"]
        electrodes.append(electrode)
    return electrodes


def generate_brain_model(settings):
    """Generate OCC brain model

    Notes
    -----

    TODO type checking

    """
    brain_region_parameters = settings['BrainRegion']
    brain_shape = brain_region_parameters["Shape"]
    brain_region = BrainRegionFactory.create(brain_region_parameters)
    brain_model = BrainGeometry(brain_region, brain_shape)
    return brain_model


def generate_model_geometry(settings):
    brain = generate_brain_model(settings)
    electrodes = generate_electrodes(settings)
    model_geometry = ModelGeometry(brain, electrodes)
    return model_geometry


def set_logger(level=logging.INFO):
    _logger.setLevel(level)
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
