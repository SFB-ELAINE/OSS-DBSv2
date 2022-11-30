from dataclasses import dataclass
from src.brainsubstance import Material
from src.dielectric_model import DielectricModel, DielectricParamters
import numpy as np


@dataclass
class WhiteMatterParameters(DielectricParamters):
    alpha: np.ndarray = np.array([0.1, 0.1, 0.3, 0.02])
    eps_delta: np.ndarray = np.array([32.0, 100.0, 4.0e4, 3.5e7])
    eps_inf: float = 4.0
    sigma: float = 0.02
    tau: np.ndarray = np.array([7.958e-12, 7.958e-9, 53.052e-6, 7.958e-3])


@dataclass
class GrayMatterParameters(DielectricParamters):
    alpha: np.ndarray = np.array([0.1, 0.15, 0.22, 0.0])
    eps_delta: np.ndarray = np.array([45.0, 400.0, 2.0e5, 4.5e7])
    eps_inf: float = 4.0
    sigma: float = 0.02
    tau: np.ndarray = np.array([7.958e-12, 15.915e-9, 106.103e-6, 5.305e-3])


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
