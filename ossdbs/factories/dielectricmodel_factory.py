
import json
from ossdbs.dielectric_model import DielectricModel
from ossdbs.dielectric_model import MaterialColeColeModel
from ossdbs.dielectric_model import ColeCole4Model
from ossdbs.dielectric_model import ConstantModel
from ossdbs.dielectric_model import MaterialConstantModel

from ossdbs.materials import Material
import numpy as np


class DielectricModelFactory:

    MATERIALS = {'Blood': Material.BLOOD,
                 'GrayMatter': Material.GRAY_MATTER,
                 'WhiteMatter': Material.WHITE_MATTER,
                 'CerebroSpinalFluid': Material.CSF,
                 'Unknown': Material.UNKNOWN
                 }

    def create(cls, dielectricum_parameters: str) -> DielectricModel:

        if 'Custom' in dielectricum_parameters['Type']:
            path = dielectricum_parameters['PathToCustomParameters']
            with open(path, 'r') as json_file:
                model_parameters = json.load(json_file)
            if 'ColeCole4' in dielectricum_parameters['Type']:
                return cls.__create_cc4_model(model_parameters)
            if 'Constant' in dielectricum_parameters['Type']:
                return cls.__create_constant_model()

        return {'ColeCole4': ColeCole4Model(),
                'Constant': ConstantModel(),
                }[dielectricum_parameters['Type']]

    def __create_cc4_model(cls, parameters):

        cc4_model = ColeCole4Model()

        for key, material in cls.MATERIALS.items():
            if key not in parameters:
                continue

            model = MaterialColeColeModel(
                        eps_inf=parameters[key]['EpsilonInfinite'],
                        sigma=parameters[key]['Sigma'],
                        alpha=np.array(parameters[key]['Alpha']),
                        eps_delta=np.array(parameters[key]['EpsilonDelta']),
                        tau=np.array(parameters[key]['Tau'])
                        )
            cc4_model.MODELS.update({material: model})

        return cc4_model

    def __create_constant_model(cls, parameters):

        constant_model = ConstantModel()
        for key, material in cls.MATERIALS.items():
            if key not in parameters:
                continue

            model = MaterialConstantModel(
                                permitivity=parameters[key]['Permitivity'],
                                conductivity=parameters[key]['Conductivity'],
                                )
            constant_model.MODELS.update({material: model})

        return constant_model
