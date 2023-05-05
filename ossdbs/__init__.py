import logging
from .__version__ import __version__

from ossdbs.factories import ElectrodesFactory
from .brain_geometry import BrainGeometry
from .model_geometry import BrainGeometry as ModelGeometry
from .nifti1Image import Nifti1Image

_logger = logging.getLogger(__name__)
_logger.addHandler(logging.NullHandler())


def generate_electrodes(settings):
    """Generate OCC electrode model

    Notes
    -----
    TODO type checking

    """
    electrodes = ElectrodesFactory.create(settings['Electrodes'],
                                          settings['EncapsulatingLayer'])
    return electrodes


def generate_brain_model(settings):
    """Generate OCC brain model

    Notes
    -----

    TODO type checking

    """
    #TODO define interface
    brain_model = BrainGeometry(settings)
    return brain_model


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
