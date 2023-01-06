from dataclasses import dataclass
from ossdbs.brainsubstance import Material
from ossdbs.dielectric_model import DielectricParamters, DielectricModel
import numpy as np


@dataclass
class WhiteMatterParameters(DielectricParamters):
    alpha: np.ndarray = np.array([0.46, 0.05])
    eps_delta: np.ndarray = np.array([4.29e3, 32.62e6])
    eps_inf: float = 10.27
    sigma: float = 0.0192
    tau: np.ndarray = np.array([3.71e-6, 8.29e-3])


@dataclass
class GrayMatterParameters(DielectricParamters):
    alpha: np.ndarray = np.array([0.45, 0.06])
    eps_delta: np.ndarray = np.array([3.68e3, 42.26e6])
    eps_inf: float = 1.0
    sigma: float = 0.027
    tau: np.ndarray = np.array([1.14e-6, 5.83e-3])


@dataclass
class CerebroSpinalFluidParameters(DielectricParamters):
    alpha: np.ndarray = np.array([0.1, 0.0, 0.0, 0.0])
    eps_delta: np.ndarray = np.array([65.0, 40.0, 0.0, 0.0])
    eps_inf: float = 4.0
    sigma: float = 2.0
    tau: np.ndarray = np.array([7.96e-12, 1.592e-9, 159.155e-6, 5.305e-3])


class ModelCreator():
    @classmethod
    def create(cls, material: Material) -> 'DielectricModel':
        material_parameters = {Material.CSF: CerebroSpinalFluidParameters,
                               Material.WHITE_MATTER: WhiteMatterParameters,
                               Material.GRAY_MATTER: GrayMatterParameters}
        return DielectricModel(material_parameters[material])
