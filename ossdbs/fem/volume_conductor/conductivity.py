# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import logging
from typing import Optional

import ngsolve
import numpy as np

from ossdbs.fem.mesh import Mesh
from ossdbs.model_geometry import BoundingBox
from ossdbs.utils.nifti1image import DiffusionTensorImage, MagneticResonanceImage

_logger = logging.getLogger(__name__)


# TODO add a linear interpolation
# linear_interpolation: bool = False


class ConductivityCF:
    """Conductivity wrapper."""

    def __init__(
        self,
        mri_image: MagneticResonanceImage,
        brain_bounding_box: BoundingBox,
        dielectric_properties: dict,
        materials: dict,
        encapsulation_layers: Optional[dict] = None,
        complex_data: bool = False,
        dti_image: DiffusionTensorImage = None,
    ) -> None:
        """Convert MRI conductivity distribution to NGSolve.

        Parameters
        ----------
        mri_image: MagneticResonanceImage
            the MRI image holding the material distribution
        dti_image: DiffusionTensorImage
            the DTI image holding the anisotropy information
        brain_bounding_box: BoundingBox
            bounding box of the brain in real space.
        materials: dict
            mapping of materials to integers in MRI
        dielectric_properties: dict
            dictionary with dielectric properties of each material
        encapsulation_layers: dict
            a dictionary containing the materials of the encapsulation layer
        complex_data: bool
            if complex arithmetic is required
        """
        if encapsulation_layers is None:
            encapsulation_layers = []
        self._dielectric_properties = dielectric_properties
        self._encapsulation_layers = encapsulation_layers
        self._is_complex = complex_data

        _logger.debug("Crop MRI image")
        brain_bounding_box_voxel = mri_image.get_voxel_bounding_box(brain_bounding_box)
        (
            self._material_distribution,
            self._mri_voxel_bounding_box,
        ) = mri_image._crop_image(brain_bounding_box_voxel)
        # account for ordering in NGSolve
        self._material_distribution = np.swapaxes(self._material_distribution, 0, 2)
        self._dti_voxel_cf = None

        if dti_image is not None:
            _logger.debug("Crop DTI image")
            brain_bounding_box_voxel = dti_image.get_voxel_bounding_box(
                brain_bounding_box
            )
            dti_data, self._dti_voxel_bounding_box = dti_image._crop_image(
                brain_bounding_box_voxel
            )
            # cast to complex datatype
            if self._is_complex:
                dti_data = dti_data.astype(self._get_datatype())
            self._dti_voxel_cf = self.create_dti_voxel_cf(
                dti_data, self._dti_voxel_bounding_box, dti_image
            )

        self._data = np.zeros(
            self._material_distribution.shape, dtype=self._get_datatype()
        )

        self._materials = materials
        self._masks = [None] * len(self.materials)
        self._trafo_cf = mri_image.trafo_cf
        # Creates a boolean mask for the indices that material is present in
        for material in self.materials:
            material_idx = self.materials[material]
            self._masks[material_idx] = np.isclose(
                self._material_distribution, material_idx
            )

    @property
    def is_complex(self) -> bool:
        """Return if complex arithmetic is used."""
        return self._is_complex

    @property
    def is_tensor(self) -> bool:
        """Return if conductivity is a tensor."""
        return self._dti_voxel_cf is not None

    @property
    def materials(self) -> dict:
        """Return dictionary of material names."""
        return self._materials

    def _get_datatype(self):
        """Return numpy datatype."""
        if self.is_complex:
            return np.complex128
        return np.float64

    def __call__(
        self,
        mesh: Mesh,
        frequency: float,
    ) -> ngsolve.CF:
        """Evaluate conductivity on mesh for a given frequency.

        Parameters
        ----------
        mesh: ossdbs.Mesh
            Mesh on which the computation is done
        frequency: float
            Stimulation frequency

        Notes
        -----
        Evaluates the conductivity and defines it per brain region.

        """
        omega = 2 * np.pi * frequency
        material_dict = {"Brain": self._distribution(omega)}
        for encapsulation_layer in self._encapsulation_layers:
            if self.is_complex:
                encapsulation_layer_properties = (
                    encapsulation_layer.dielectric_properties.complex_conductivity(
                        omega
                    )
                )
            else:
                encapsulation_layer_properties = (
                    encapsulation_layer.dielectric_properties.conductivity(omega)
                )
            if self._dti_voxel_cf is not None:
                # reshape CoefficientFunction for isotropic encapsulation layer
                encapsulation_layer_properties_cf = ngsolve.CoefficientFunction(
                    (
                        encapsulation_layer_properties,
                        0,
                        0,
                        0,
                        encapsulation_layer_properties,
                        0,
                        0,
                        0,
                        encapsulation_layer_properties,
                    ),
                    dims=(3, 3),
                )
            else:
                encapsulation_layer_properties_cf = encapsulation_layer_properties
            material_dict[encapsulation_layer.name] = encapsulation_layer_properties_cf
        return mesh.material_coefficients(material_dict)

    def _distribution(
        self,
        omega: float,
    ) -> ngsolve.VoxelCoefficient:
        """Return the conductivity distribution at the given frequency.

        Parameters
        ----------
        omega : float
            Angular frequency

        Returns
        -------
        ngsolve.VoxelCoefficient
            Data structure representing the conductivity distribution in space.
        """
        for material in self.materials:
            material_idx = self.materials[material]
            # Sets conductivity values based on
            # what indices in self._data are material_idx
            if self.is_complex:
                self._data[self._masks[material_idx]] = self._dielectric_properties[
                    material
                ].complex_conductivity(omega)
            else:
                self._data[self._masks[material_idx]] = self._dielectric_properties[
                    material
                ].conductivity(omega)
        start = self._mri_voxel_bounding_box.start
        end = self._mri_voxel_bounding_box.end
        if self._dti_voxel_cf is None:
            return ngsolve.VoxelCoefficient(
                start, end, self._data, False, trafocf=self._trafo_cf
            )
        return self._dti_voxel_cf * ngsolve.VoxelCoefficient(
            start, end, self._data, False, trafocf=self._trafo_cf
        )

    def material_distribution(self, mesh: Mesh) -> ngsolve.VoxelCoefficient:
        """Return MRI image projected onto mesh.

        Notes
        -----
        The encapsulation layer material is set to a fixed value.
        """
        start, end = (
            self._mri_voxel_bounding_box.start,
            self._mri_voxel_bounding_box.end,
        )
        material_voxelcf = ngsolve.VoxelCoefficient(
            start,
            end,
            self._material_distribution,
            linear=False,
            trafocf=self._trafo_cf,
        )
        material_dict = {"Brain": material_voxelcf}
        for encapsulation_layer in self._encapsulation_layers:
            material_dict[encapsulation_layer.name] = self.materials[
                encapsulation_layer.material
            ]
        return mesh.material_coefficients(material_dict)

    def create_dti_voxel_cf(self, dti_data, dti_voxel_bounding_box, dti_image):
        """Map DTI image on NGSolve VoxelCoefficient function.

        Notes
        -----
        The DTI images assume symmetry of the tensor.
        This is not (yet) the case in NGSolve.
        """
        start, end = dti_voxel_bounding_box.start, dti_voxel_bounding_box.end
        dti_flat_matrix = []
        for _component, index in dti_image.components.items():
            dti_component_data = dti_data[:, :, :, index]
            # account for ordering in NGSolve
            dti_component_data = np.swapaxes(dti_component_data, 0, 2)
            dti_flat_matrix.append(
                ngsolve.VoxelCoefficient(
                    start,
                    end,
                    dti_component_data,
                    linear=False,
                    trafocf=dti_image.trafo_cf,
                )
            )
        return ngsolve.CoefficientFunction(tuple(dti_flat_matrix), dims=(3, 3))

    @property
    def dti_voxel_distribution(self) -> ngsolve.VoxelCoefficient:
        """Return DTI data as VoxelCoefficient before scaling by conductivity."""
        return self._dti_voxel_cf

    @property
    def dielectric_properties(self) -> dict:
        """Return dielectric properties."""
        return self._dielectric_properties
