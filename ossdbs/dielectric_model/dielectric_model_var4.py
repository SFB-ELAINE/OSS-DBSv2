from dataclasses import dataclass
from ossdbs.materials import Material
from ossdbs.dielectric_model import DielectricParamters, DielectricModel
import numpy as np


@dataclass
class WhiteMatterParameters(DielectricParamters):
    alpha: np.ndarray = np.array([0.02])
    eps_delta: np.ndarray = np.array([49.03e3])
    eps_inf: float = 7.77
    sigma: float = 0.0611
    tau: np.ndarray = np.array([159.77e-6])


@dataclass
class GrayMatterParameters(DielectricParamters):
    alpha: np.ndarray = np.array([0.0])
    eps_delta: np.ndarray = np.array([95.55e3])
    eps_inf: float = 11.92
    sigma: float = 0.107
    tau: np.ndarray = np.array([170.91e-6])


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
