
from ossdbs.materials import Material
from ossdbs.brain_imaging.mri import MagneticResonanceImage
from ossdbs.dielectric_model import ColeColeFourModelFactory
import numpy as np
import ngsolve


class Conductivity:
    """Represents the conductivity distribution by magnetic resonance imaging.

    Parameters
    ----------
    mri: MagneticResonanceImage
        Image which represents the distributiion of brain substances.
    """

    __MATERIALS = [Material.CSF, Material.GRAY_MATTER, Material.WHITE_MATTER]

    def __init__(self,
                 mri: MagneticResonanceImage,
                 complex_datatype: bool = False) -> None:
        self.__mri = mri
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
        colecole_model = ColeColeFourModelFactory.create(Material.GRAY_MATTER)
        default = colecole_model.conductivity(omega)
        data = np.full(self.__mri.xyz_shape(), default)

        for material in self.__MATERIALS:
            position = self.__mri.material_distribution(material=material)
            colecole_model = ColeColeFourModelFactory.create(material)
            data[position.data] = colecole_model.conductivity(omega)

        start, end = self.__mri.bounding_box()

        if not self.__complex:
            return ngsolve.VoxelCoefficient(start, end, np.real(data), False)

        return ngsolve.VoxelCoefficient(start, end, data, False)
