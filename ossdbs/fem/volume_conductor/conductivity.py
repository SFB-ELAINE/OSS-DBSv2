from ossdbs.utils.nifti1image import (MagneticResonanceImage,
                                      DiffusionTensorImage)
from ossdbs.model_geometry import BoundingBox
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

    dielectric_properties : dict
        Dictionary with dielectric properties of each material

    encapsulation_layers : dict
        A dictionary containing the materials of the encapsulation layer
    """

    def __init__(self,
                 mri_image: MagneticResonanceImage,
                 brain_bounding_box: BoundingBox,
                 dielectric_properties: dict,
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
        self._dti_voxel_cf = None

        if dti_image is not None:
            _logger.debug("Crop DTI image")
            brain_bounding_box_voxel = dti_image.get_voxel_bounding_box(brain_bounding_box)
            self._dti_data, self._dti_voxel_bounding_box = dti_image._crop_image(brain_bounding_box_voxel)
            self._dti_voxel_cf = self.create_dti_voxel_cf(self._dti_data, self._dti_voxel_bounding_box, dti_image)

        self._dielectric_properties = dielectric_properties
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
                encapsulation_layer_properties = encapsulation_layer.dielectric_properties.complex_conductivity(omega)
            else:
                encapsulation_layer_properties = encapsulation_layer.dielectric_properties.conductivity(omega)
            if self._dti_voxel_cf is not None:
                # reshape CoefficientFunction for isotropic encapsulation layer
                encapsulation_layer_properties_cf =\
                        ngsolve.CoefficientFunction((encapsulation_layer_properties, 0, 0,
                                                     0, encapsulation_layer_properties, 0,
                                                     0, 0, encapsulation_layer_properties),
                                                    dims=(3, 3))
            else:
                encapsulation_layer_properties_cf = encapsulation_layer_properties
            material_dict[encapsulation_layer.name] = encapsulation_layer_properties_cf
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
                self._data[self._masks[material_idx]] = self._dielectric_properties[material].complex_conductivity(omega)
            else:
                self._data[self._masks[material_idx]] = self._dielectric_properties[material].conductivity(omega)
        start = self._mri_voxel_bounding_box.start
        end = self._mri_voxel_bounding_box.end
        if self._dti_voxel_cf is None:
            return ngsolve.VoxelCoefficient(start, end, self._data, False, trafocf=self._trafo_cf)
        return self._dti_voxel_cf * ngsolve.VoxelCoefficient(start, end, self._data, False, trafocf=self._trafo_cf)

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
