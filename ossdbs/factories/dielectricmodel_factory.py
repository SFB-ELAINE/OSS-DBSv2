
import json
from ossdbs.dielectric_model import DielectricModel
from ossdbs.dielectric_model import ColeColeModel
from ossdbs.dielectric_model import ColeCole4Model
from ossdbs.dielectric_model import ColeCole4ModelCustom
from ossdbs.dielectric_model import ConstantModel
from ossdbs.dielectric_model import ConstantDielectricModel

from ossdbs.dielectric_model import ConstantModelCustom
from ossdbs.materials import Material
import numpy as np

class DielectricModelFactory:

    def create(cls, dielectricum_parameters: str) -> DielectricModel:

        if 'Custom' in dielectricum_parameters['Type']:
            path = dielectricum_parameters['PathToCustomParameters']
            with open(path, 'r') as json_file:
                model_parameters = json.load(json_file)
            if 'ColeCole4' in dielectricum_parameters['Type']:
                return cls.create_cole_cole_model(model_parameters)
            if 'Constant' in dielectricum_parameters['Type']:
                return cls.create_constant_model()

        return {'ColeCole4': ColeCole4Model(),
                'Constant': ConstantModel(),
                }[dielectricum_parameters['Type']]

    def create_cole_cole_model(cls, parameters):

        wm_paramters = parameters['WhiteMatter']
        wm_model = ColeColeModel(eps_inf=wm_paramters['EpsilonInfinite'],
                                 sigma=wm_paramters['Sigma'],
                                 alpha=np.array(wm_paramters['Alpha']),
                                 eps_delta=np.array(wm_paramters['EpsilonDelta']),
                                 tau=np.array(wm_paramters['Tau'])
                                 )

        gm_paramters = parameters['GrayMatter']
        gm_model = ColeColeModel(eps_inf=gm_paramters['EpsilonInfinite'],
                                 sigma=gm_paramters['Sigma'],
                                 alpha=np.array(gm_paramters['Alpha']),
                                 eps_delta=np.array(gm_paramters['EpsilonDelta']),
                                 tau=np.array(gm_paramters['Tau'])
                                 )

        csf_paramters = parameters['CerebroSpinalFluid']
        csf_model = ColeColeModel(eps_inf=csf_paramters['EpsilonInfinite'],
                                  sigma=csf_paramters['Sigma'],
                                  alpha=np.array(csf_paramters['Alpha']),
                                  eps_delta=np.array(csf_paramters['EpsilonDelta']),
                                  tau=np.array(csf_paramters['Tau'])
                                  )

        blood_paramters = parameters['Blood']
        blood_model = ColeColeModel(eps_inf=blood_paramters['EpsilonInfinite'],
                                    sigma=blood_paramters['Sigma'],
                                    alpha=np.array(blood_paramters['Alpha']),
                                    eps_delta=np.array(blood_paramters['EpsilonDelta']),
                                    tau=np.array(blood_paramters['Tau'])
                                    )


        model = ColeCole4Model()
        model.MODELS = {Material.CSF: csf_model,
                        Material.BLOOD: blood_model,
                        Material.GRAY_MATTER: gm_model,
                        Material.WHITE_MATTER: wm_model}
        return model

    def create_constant_model(cls, parameters):
        wm_paramters = parameters['WhiteMatter']
        wm_model = ConstantDielectricModel(
                                permitivity=wm_paramters['Permitivity'],
                                conductivity=wm_paramters['Conductivity'],
                                )

        gm_paramters = parameters['GrayMatter']
        gm_model = ConstantDielectricModel(
                                permitivity=gm_paramters['Permitivity'],
                                conductivity=gm_paramters['Conductivity'],
                                )

        csf_paramters = parameters['CerebroSpinalFluid']
        csf_model = ConstantDielectricModel(
                                permitivity=csf_paramters['Permitivity'],
                                conductivity=csf_paramters['Conductivity'],
                                )

        blood_paramters = parameters['Blood']
        blood_model = ConstantDielectricModel(
                                permitivity=blood_paramters['Permitivity'],
                                conductivity=blood_paramters['Conductivity'],
                                )

        model = ConstantModel()
        model.MODELS = {Material.CSF: csf_model,
                        Material.BLOOD: blood_model,
                        Material.GRAY_MATTER: gm_model,
                        Material.WHITE_MATTER: wm_model}
        return model
