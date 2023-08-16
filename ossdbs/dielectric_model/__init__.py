from .dielectric_model import DielectricModel
from .colecole4 import MaterialColeColeModel
from .colecole4 import ColeCole4Model
from .colecole3 import ColeCole3Model
from .constant import ConstantModel
from .constant import MaterialConstantModel
from ossdbs.utils.materials import MATERIALS
import numpy as np

import logging
_logger = logging.getLogger(__name__)


def create_cc4_model(parameters: dict):

    cc4_model = ColeCole4Model()

    for key, material in MATERIALS.items():
        if key not in parameters:
            _logger.warning("Parameter {} not a valid parameter for ColeCole4Model".format(key))
            continue

        model = MaterialColeColeModel(eps_inf=parameters[key]['EpsilonInfinite'],
                                      sigma=parameters[key]['Sigma'],
                                      alpha=np.array(parameters[key]['Alpha']),
                                      eps_delta=np.array(parameters[key]['EpsilonDelta']),
                                      tau=np.array(parameters[key]['Tau'])
                                      )
        cc4_model.MODELS.update({material: model})

    return cc4_model


def create_cc3_model(parameters: dict):

    cc3_model = ColeCole3Model()

    for key, material in MATERIALS.items():
        if key not in parameters:
            _logger.warning("Parameter {} not a valid parameter for ColeCole3Model".format(key))
            continue

        model = MaterialColeColeModel(eps_inf=parameters[key]['EpsilonInfinite'],
                                      sigma=parameters[key]['Sigma'],
                                      alpha=np.array(parameters[key]['Alpha']),
                                      eps_delta=np.array(parameters[key]['EpsilonDelta']),
                                      tau=np.array(parameters[key]['Tau'])
                                      )
        cc3_model.MODELS.update({material: model})

    return cc3_model


def create_constant_model(parameters: dict):
    constant_model = ConstantModel()
    for key, material in MATERIALS.items():
        if key not in parameters:
            _logger.warning("Parameter {} not a valid parameter for ConstantModel".format(key))
            continue

        model = MaterialConstantModel(permittivity=parameters[key]['Permittivity'],
                                      conductivity=parameters[key]['Conductivity'],
                                      )
        constant_model.MODELS.update({material: model})
    return constant_model


DIELECTRIC_MODELS = {'ColeCole4': ColeCole4Model(),
                     'Constant': ConstantModel()
                     }


__all__ = ('DielectricModel',
           'ColeCole4Model',
           'MaterialColeColeModel',
           'ConstantModel',
           'MaterialConstantModel',
           )
