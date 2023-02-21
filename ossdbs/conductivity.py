
from ossdbs.materials import Material
from ossdbs.dielectric_model import CerebroSpinalFluidModel
from ossdbs.dielectric_model import GrayMatterModel
from ossdbs.dielectric_model import WhiteMatterModel
import numpy as np
import ngsolve


class Conductivity:
    """Represents the conductivity distribution by magnetic resonance imaging.

    Parameters
    ----------
    mri: MagneticResonanceImage
        Image which represents the distributiion of brain substances.
    """

    def __init__(self,
                 material_distribution: np.ndarray,
                 bounding_box: tuple,
                 complex_datatype: bool = False) -> None:
        self.__material_distribution = material_distribution
        self.__bounding_box = bounding_box
        self.__complex = complex_datatype

    def distribution(self, frequency: float) -> ngsolve.VoxelCoefficient:
        """Return the conductivity distribution by the given frequency.

        Parameters
        ----------
        frequency : float

        Returns
        -------
        ngsolve.VoxelCoefficient
            Data structure representing the conductivity distribution and the
            location in space.
        """

        omega = 2 * np.pi * frequency
        default = GrayMatterModel().conductivity(omega)
        data = np.full(self.__material_distribution.shape, default)

        pos_csf = self.__material_distribution == Material.CSF
        pos_wm = self.__material_distribution == Material.WHITE_MATTER
        pos_gm = self.__material_distribution == Material.GRAY_MATTER

        data[pos_csf] = CerebroSpinalFluidModel().conductivity(omega)
        data[pos_wm] = WhiteMatterModel().conductivity(omega)
        data[pos_gm] = GrayMatterModel().conductivity(omega)

        m_to_mm = 1e3
        data = data / m_to_mm
        start, end = self.__bounding_box

        if not self.__complex:
            return ngsolve.VoxelCoefficient(start, end, np.real(data), False)

        return ngsolve.VoxelCoefficient(start, end, data, False)
