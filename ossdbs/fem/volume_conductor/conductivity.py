from ossdbs.utils.nifti1image import (MagneticResonanceImage,
                                      DiffusionTensorImage)
from ossdbs.model_geometry import BoundingBox
from ossdbs.dielectric_model import DielectricModel
from ossdbs.fem.mesh import Mesh
import numpy as np
import ngsolve


class ConductivityCF:
    """Represents the conductivity distribution by magnetic resonance imaging.

    Attributes
    ----------
    material_distribution : np.ndarray
        Matrix represents coding for the distribution of different materials.

    bounding_box : BoundingBox
        Represents a cuboid aligned with the cartesian axis.

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
        self._material_distribution, self._bounding_box = self._crop_mri_image(mri_image, brain_bounding_box)
        if dti_image is not None:
            self._diffusion, self._bounding_box = self._crop_dti_image(dti_image, brain_bounding_box)
            # TODO implement slicing
        self._dielectric_model = dielectric_model
        self._encapsulation_layers = encapsulation_layers
        self._is_complex = complex_data

        self._data = np.zeros(self._material_distribution.shape, dtype=self._get_datatype())
        self._materials = materials
        self._masks = [None] * len(self._materials)
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
        # TODO integrate diffusion if applicable
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
            if self.is_complex:
                self._data[self._masks[material_idx]] = self._dielectric_model.complex_conductivity(material, omega)
            else:
                self._data[self._masks[material_idx]] = self._dielectric_model.conductivity(material, omega)

        # transform conductivity [S/m] to [S/mm] since the geometry is
        # measured in mm
        values = self._data * 1e-3
        start, end = self._bounding_box.start, self._bounding_box.end

        return ngsolve.VoxelCoefficient(start, end, values, False)

    def _crop_mri_image(self, nifti, brain_bounding_box) -> np.ndarray:
        """Crop the Nifti image of the conductivity to match the brain geometry.

        Returns
        -------
        TODO
        """

        voxel_size = nifti.voxel_size
        offset = nifti.offset
        bbox = brain_bounding_box.intersection(nifti.bounding_box)
        start, end = bbox.start, bbox.end
        start_index = np.floor(np.subtract(start, offset) / voxel_size)
        x_s, y_s, z_s = start_index.astype(int)
        end_index = np.ceil(np.subtract(end, offset) / voxel_size)
        x_e, y_e, z_e = end_index.astype(int)

        data = nifti.data_map[x_s:x_e, y_s:y_e, z_s:z_e]
        new_start = tuple(start_index * voxel_size + offset)
        new_end = tuple(end_index * voxel_size + offset)

        bounding_box = BoundingBox(start=new_start, end=new_end)

        return data, bounding_box

    def _crop_dti_image(self, nifti, brain_bounding_box) -> np.ndarray:
        """Crop the Nifti image of the diffusion to match the brain geometry.

        Returns
        -------
        TODO
        """

        voxel_size = nifti.voxel_size
        offset = nifti.offset
        bbox = brain_bounding_box.intersection(nifti.bounding_box)
        start, end = bbox.start, bbox.end
        start_index = np.floor(np.subtract(start, offset) / voxel_size)
        x_s, y_s, z_s = start_index.astype(int)
        end_index = np.ceil(np.subtract(end, offset) / voxel_size)
        x_e, y_e, z_e = end_index.astype(int)

        #TODO
        # data = nifti.data_map[x_s:x_e, y_s:y_e, z_s:z_e]
        datatmp = nifti.diffusion()
        data = datatmp[x_s:x_e, y_s:y_e, z_s:z_e]

        new_start = tuple(start_index * voxel_size + offset)
        new_end = tuple(end_index * voxel_size + offset)

        bounding_box = BoundingBox(start=new_start, end=new_end)

        return data, bounding_box
