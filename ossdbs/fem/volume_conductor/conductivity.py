from ossdbs.utils.nifti1image import (MagneticResonanceImage,
                                      DiffusionTensorImage)
from ossdbs.model_geometry import BoundingBox
from ossdbs.dielectric_model import DielectricModel
from ossdbs.fem.mesh import Mesh
import numpy as np
import ngsolve
import logging

_logger = logging.getLogger(__name__)


class ConductivityCF:
    """Represents the conductivity distribution by magnetic resonance imaging.

    Attributes
    ----------
    material_distribution : np.ndarray
        Matrix represents coding for the distribution of different materials.

    bounding_box : BoundingBox
        Represents a cuboid in real space.

    dielectric_model : DielectricModel
        Model for the dielectric spectrum of a tissues.

    encapsulation_layers : dict
        A dictionary containing the materials of the encapsulation layer
    """

    def __init__(self,
                 mri_image: MagneticResonanceImage,
                 brain_bounding_box: BoundingBox,
                 dielectric_model: DielectricModel,
                 materials: dict,
                 encapsulation_layers=[],  # TODO type hint
                 complex_data: bool = False,
                 dti_image: DiffusionTensorImage = None
                 ) -> None:

        _logger.debug("Crop MRI image")
        brain_bounding_box_voxel = mri_image.get_voxel_bounding_box(brain_bounding_box)
        self._material_distribution, self._mri_voxel_bounding_box = mri_image._crop_image(brain_bounding_box_voxel)
        # account for ordering in NGSolve
        self._material_distribution = np.swapaxes(self._material_distribution, 0, 2)

        if dti_image is not None:
            _logger.debug("Crop DTI image")
            brain_bounding_box_voxel = dti_image.get_voxel_bounding_box(brain_bounding_box)
            self._dti_data, self._dti_voxel_bounding_box = dti_image._crop_image(brain_bounding_box_voxel)
            self._dti_voxel_cf = self.create_dti_voxel_cf(self._dti_data, self._dti_voxel_bounding_box, dti_image)

        self._dielectric_model = dielectric_model
        self._encapsulation_layers = encapsulation_layers
        self._is_complex = complex_data
        self._data = np.zeros(self._material_distribution.shape, dtype=self._get_datatype())
        self._materials = materials
        self._masks = [None] * len(self._materials)
        self._trafo_cf = mri_image.trafo_cf
        # Creates a boolean mask for the indices that material is present in
        for material in self._materials:
            material_idx = self._materials[material]
            self._masks[material_idx] = self._material_distribution == material_idx

    @property
    def is_complex(self) -> bool:
        return self._is_complex

    def _get_datatype(self):
        """Return numpy datatype
        """
        if self.is_complex:
            return np.complex128
        return np.float64

    def __call__(self,
                 mesh: Mesh,
                 frequency: float,
                 ) -> ngsolve.CF:
        omega = 2 * np.pi * frequency
        material_dict = {"Brain": self._distribution(omega)}
        for encapsulation_layer in self._encapsulation_layers:
            if self.is_complex:
                material_dict[encapsulation_layer.name] = encapsulation_layer.dielectric_model.complex_conductivity(encapsulation_layer.material, omega)
            else:
                material_dict[encapsulation_layer.name] = encapsulation_layer.dielectric_model.conductivity(encapsulation_layer.material, omega)
        return mesh.material_coefficients(material_dict)

    def _distribution(self,
                      omega: float,
                      ) -> ngsolve.VoxelCoefficient:
        """Return the conductivity distribution at the given frequency.

        Parameters
        ----------
        omega : float

        Returns
        -------
        ngsolve.VoxelCoefficient
            Data structure representing the conductivity distribution in space.
        """
        for material in self._materials:
            material_idx = self._materials[material]
            # Sets conductivity values based on what indices in self._data are material_idx
            if self.is_complex:
                self._data[self._masks[material_idx]] = self._dielectric_model.complex_conductivity(material, omega)
            else:
                self._data[self._masks[material_idx]] = self._dielectric_model.conductivity(material, omega)
        start = self._mri_voxel_bounding_box.start
        end = self._mri_voxel_bounding_box.end
        return ngsolve.VoxelCoefficient(start, end, self._data, False, trafocf=self._trafo_cf)

    def material_distribution(self, mesh: Mesh) -> ngsolve.VoxelCoefficient:
        """Return MRI image projected onto mesh

        Notes
        -----

        The encapsulation layer material is set to a fixed value.
        """
        start, end = self._mri_voxel_bounding_box.start, self._mri_voxel_bounding_box.end
        materials = ngsolve.VoxelCoefficient(start, end, self._material_distribution, linear=False, trafocf=self._trafo_cf)
        material_dict = {"Brain": materials}
        for encapsulation_layer in self._encapsulation_layers:
            material_dict[encapsulation_layer.name] = self._materials[encapsulation_layer.material]
        return mesh.material_coefficients(material_dict)

    def create_dti_voxel_cf(self, dti_data, dti_voxel_bounding_box, dti_image):
        start, end = dti_voxel_bounding_box.start, dti_voxel_bounding_box.end
        dti_flat_matrix = []
        for component, index in dti_image.components.items():
            dti_component_data = dti_data[:, :, :, index]
            # account for ordering in NGSolve
            dti_component_data = np.swapaxes(dti_component_data, 0, 2)
            dti_flat_matrix.append(ngsolve.VoxelCoefficient(start, end, dti_component_data, linear=False, trafocf=dti_image.trafo_cf))
        return ngsolve.CoefficientFunction(tuple(dti_flat_matrix), dims=(3, 3))
