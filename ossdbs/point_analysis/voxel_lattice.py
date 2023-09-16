import numpy as np
import nibabel
import h5py
from .point_model import PointModel
from .time_results import TimeResult
import nibabel as nib


class VoxelLattice(PointModel):
    """Matrix of grid points in the centers of MRI voxels

    Attributes
    ----------

    imp_coord : np.ndarray
        (3,) np.ndarray specifying the implantation coordinate in mm space (this will be the center of the grid)

    affine : np.ndarray
        (4,4) np.ndarray specifying the affine of the MRI image

    shape : np.ndarray
        (3,) np.ndarray specifying the number of points in each direction (x, y, z)
        Each dimension must be odd. Otherwise, the point corresponding to the voxel containing the implantation
        coordinate will not be at the center of the grid.

    """

    def __init__(self,
                 imp_coord: np.ndarray,
                 affine: np.ndarray,
                 shape: np.ndarray,
                 header: nibabel.Nifti1Header
                 ) -> None:
        self._imp_coord = imp_coord
        self._affine = affine
        self._shape = shape
        self._header = header

        # Check on dimension condition on shape input
        if np.sum(shape[shape / 2 == 0]) > 0:
            raise Exception("Each dimension of the shape argument must be odd number")

        self._coordinates = self._initialize_coordinates()

    def _initialize_coordinates(self) -> np.ndarray:
        """Generates grids of points in the MRI voxels.

        Returns
        -------
        np.ndarray
        """

        # CALCULATION OF GRID CENTER

        # Calculate the voxel index corresponding to the implantation coordinate
        inv_aff = np.linalg.inv(self.affine)
        imp_vox_idx = inv_aff @ np.append(self.imp_coord, 1)

        # DOUBLE CHECK AFFINE MAPS COORDINATE TO VOXEL CENTER
        # "[...] a coordinate such as (i,j,k) == np.round((i,j,k)) expresses the center of a voxel"
        # https://spinalcordtoolbox.com/overview/concepts/spaces-and-coordinates.html
        # Therefore the implantation coordinate will lie within the voxel indexed by imp_vox_idx
        imp_vox_idx = np.round(imp_vox_idx)

        # This is the coordinate of the center of the voxel that contains the implantation coordinate
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

        _coordinates = (self.affine @ points)[0:3, :].T - self.affine[:3, 3] + imp_vox_center_coord[0:3]

        # Apply affine to homogenized points, center around center, and unhomogenize
        return _coordinates

    def save_as_nifti(self, scalar_field, filename, binarize=False, activation_threshold=None):
        # TODO add support for upsampled lattice

        if binarize and activation_threshold is None:
            raise ValueError("Need to provide activation_threshold to binarize")
        # get affine of the segmented MRI image to use as a template
        affine_grid = self.affine.copy()
        affine_grid[0:3, 3] = [self.coordinates[0][0], self.coordinates[0][1], self.coordinates[0][2]]

        # Assuming data is in the same format as it was generated,
        # you can just reshape it
        nifti_grid = scalar_field.reshape(self._shape)

        nifti_output = np.zeros(nifti_grid.shape, float)
        if binarize:
            nifti_output[nifti_grid >= activation_threshold] = 1
            nifti_output[nifti_grid < activation_threshold] = 0
        else:
            nifti_output = nifti_grid  # V/mm

        nib.save(nib.Nifti1Image(nifti_output, affine_grid, self._header), filename)

    def _gen_grid(self):
        """ Return list of ndarrays (coordinate matrices from coordinate vectors).
        """
        begin = -((self.shape - 1) / 2).astype(int)
        end = -begin
        base_x = np.linspace(begin[0], end[0], self.shape[0])
        base_y = np.linspace(begin[1], end[1], self.shape[1])
        base_z = np.linspace(begin[2], end[2], self.shape[2])

        # iteration ordering is z,y,x
        return np.meshgrid(base_x, base_y, base_z, indexing='ij')

    # Time result still needs to be implemented
    def save(self, data: TimeResult, file_name: str) -> None:
        with h5py.File(file_name, "w") as file:
            self._write_file(data, file)

    def set_location_names(self, names: np.ndarray) -> None:
        self._location = names

    def _write_file(self, data, file):
        file.create_dataset('TimeSteps[s]', data=data.time_steps)
        file.create_dataset('Points[mm]', data=data.points)
        file.create_dataset('Location', data=self._location.astype('S'))
        file.create_dataset('Potential[V]', data=data.potential)
        file.create_dataset('Current_density[A|m2]', data=data.current_density)

    @property
    def shape(self):
        return self._shape

    @property
    def imp_coord(self):
        return self._imp_coord

    @property
    def affine(self):
        return self._affine