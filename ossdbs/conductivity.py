
from ossdbs.materials import Material
from ossdbs.dielectric_model import CerebroSpinalFluidModel
from ossdbs.dielectric_model import GrayMatterModel
from ossdbs.dielectric_model import WhiteMatterModel
from ossdbs.dielectric_model import BloodModel
import numpy as np
import ngsolve

from ossdbs.bounding_box import BoundingBox


class Conductivity:
    """Represents the conductivity distribution by magnetic resonance imaging.

    """

    def __init__(self,
                 material_distribution: np.ndarray,
                 bounding_box: BoundingBox,
                 ) -> None:
        self.__material_distribution = material_distribution
        self.__bounding_box = bounding_box

    def distribution(self,
                     frequency: float,
                     complex_data: bool = False
                     ) -> ngsolve.VoxelCoefficient:
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
        pos_blood = self.__material_distribution == Material.GRAY_MATTER

        data[pos_csf] = CerebroSpinalFluidModel().conductivity(omega)
        data[pos_wm] = WhiteMatterModel().conductivity(omega)
        data[pos_gm] = GrayMatterModel().conductivity(omega)
        data[pos_blood] = BloodModel().conductivity(omega)

        # transform conductivity [S/m] to [S/mm] since the geometry is
        # measured in mm
        values = data * 1e-3
        start, end = self.__bounding_box.start, self.__bounding_box.end

        if not complex_data:
            return ngsolve.VoxelCoefficient(start, end, np.real(values), False)

        return ngsolve.VoxelCoefficient(start, end, values, False)
