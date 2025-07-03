# Copyright 2023, 2024 Konstantin Butenko, Jan Philipp Payonk, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

from typing import Optional

import nibabel
import numpy as np

from .point_model import PointModel


class VoxelLattice(PointModel):
    """Matrix of grid points in the centers of MRI voxels.

    Attributes
    ----------
    imp_coord : np.ndarray
        Implantation coordinates in mm space, this will be the center of the grid
    affine : np.ndarray
        Affine transformation of the MRI image.
    shape : np.ndarray
        Number of points in each direction (x, y, z).

    Notes
    -----
    Each dimension of shape must be odd to ensure the grid is centered.
    """

    def __init__(
        self,
        imp_coord: np.ndarray,
        affine: np.ndarray,
        shape: np.ndarray,
        header: nibabel.Nifti1Header,
        export_field: bool = True,
    ) -> None:
        self._imp_coord = imp_coord
        self._affine = affine
        self._shape = shape
        self._header = header
        self._collapse_VTA = False
        self._export_field = export_field

        # TODO is that correct?
        self._location = None

        # Check on dimension condition on shape input
        if np.any(self._shape % 2 == 0):
            raise ValueError("Each dimension of the shape must be an odd number")

        self._coordinates = self._initialize_coordinates()

        # identifiers
        self._name = "VoxelLattice"

        # never compute time-domain signal by default
        self._time_domain_conversion = False

    def _initialize_coordinates(self) -> np.ndarray:
        """Generate coordinates for voxel lattice centered at implantation coordinate.

        Returns
        -------
        np.ndarray
            Array of voxel center coordinates in MRI space.
        """
        # CALCULATION OF GRID CENTER

        # Calculate the voxel index corresponding to the implantation coordinate
        inv_aff = np.linalg.inv(self.affine)
        imp_vox_idx = inv_aff @ np.append(self.imp_coord, 1)

        # DOUBLE CHECK AFFINE MAPS COORDINATE TO VOXEL CENTER
        # a coordinate such as (i,j,k) == np.round((i,j,k)) is the center of a voxel
        # https://spinalcordtoolbox.com/overview/concepts/spaces-and-coordinates.html
        # Therefore the implantation coordinate will lie within
        # the voxel indexed by imp_vox_idx
        imp_vox_idx = np.round(imp_vox_idx)

        # This is the coordinate of the center of the voxel
        # that contains the implantation coordinate
        imp_vox_center_coord = self.affine @ imp_vox_idx.T

        # GRID CONSTRUCTION AROUND GRID CENTER ###
        base_xg, base_yg, base_zg = self._gen_grid()
        base_xg = base_xg.reshape((np.prod(self.shape), 1))
        base_yg = base_yg.reshape((np.prod(self.shape), 1))
        base_zg = base_zg.reshape((np.prod(self.shape), 1))

        # These points are the centers of the voxels to
        # transform from voxel space into coordinate space
        points = np.concatenate([base_xg, base_yg, base_zg], axis=1)
        # Homogenize points for affine transformation
        points = np.concatenate([points, np.ones_like(base_xg)], axis=1).T

        _coordinates = (
            (self.affine @ points)[0:3, :].T
            - self.affine[:3, 3]
            + imp_vox_center_coord[0:3]
        )

        # Apply affine to homogenized points, center around center, and unhomogenize
        return _coordinates

    def save_as_nifti(
        self,
        scalar_field: np.ndarray,
        filename: str,
        binarize: bool = False,
        activation_threshold: Optional[float] = None,
    ):
        """Save scalar field in abstract orthogonal space in nifti format.

        Parameters
        ----------
        scalar_field : numpy.ndarray
            Nx1 array of scalar values on the lattice
        filename: str
            Name for the nifti file that should contain full path
        binarize: bool
            Choose to threshold the scalar field and save the binarized result
        activation_threshold: float
            Activation threshold for VTA estimate

        Notes
        -----
        TODO add support for upsampled lattice
        """
        if binarize and activation_threshold is None:
            raise ValueError("Need to provide activation_threshold to binarize")
        # get affine of the segmented MRI image to use as a template
        affine_grid = self.affine.copy()
        affine_grid[0:3, 3] = [
            self.coordinates[0][0],
            self.coordinates[0][1],
            self.coordinates[0][2],
        ]

        # Assuming data is in the same format as it was generated,
        # you can just reshape it
        nifti_grid = scalar_field.reshape(self._shape)

        nifti_output = np.zeros(nifti_grid.shape, float)
        if binarize:
            nifti_output[nifti_grid >= activation_threshold] = 1
            nifti_output[nifti_grid < activation_threshold] = 0
        else:
            nifti_output = nifti_grid  # V/mm

        nibabel.save(
            nibabel.Nifti1Image(nifti_output, affine_grid, self._header), filename
        )

    def _gen_grid(self):
        """Return list of ndarrays (coordinate matrices from coordinate vectors)."""
        begin = -((self.shape - 1) / 2).astype(int)
        end = -begin
        base_x = np.linspace(begin[0], end[0], self.shape[0])
        base_y = np.linspace(begin[1], end[1], self.shape[1])
        base_z = np.linspace(begin[2], end[2], self.shape[2])

        # iteration ordering is z,y,x
        return np.meshgrid(base_x, base_y, base_z, indexing="ij")

    @property
    def shape(self):
        """Shape of MRI data."""
        return self._shape

    @property
    def imp_coord(self):
        """Implantation coordinates."""
        return self._imp_coord

    @property
    def affine(self):
        """Affine transformation."""
        return self._affine
