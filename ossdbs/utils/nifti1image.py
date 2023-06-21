from ossdbs.model_geometry import BoundingBox
import nibabel
import numpy as np


class Nifti1Image:
    """Intrerface for Nifti1 image.

    Attributes
    ----------
    file_path : str
        File path of a Nifti1 image.

    Notes
    -----

    TODO: add link to Nibabel documentation

    """

    _N_DIMENSION = 3

    def __init__(self, file_path: str) -> None:
        self._image = self._load_image(file_path)

    @property
    def data_map(self) -> np.memmap:
        """Return the data of the nifti1 image.

        Returns
        -------
        np.memmap
        """
        return self._image.get_fdata()

    @property
    def bounding_box(self) -> BoundingBox:
        """Return the bounding box of the voxel data.

        Returns
        -------
        BoundingBox
        """
        start = self.offset
        shape = np.array(self.xyz_shape, dtype=np.float64)
        ends = start + shape * self.voxel_size
        return BoundingBox(tuple(start), tuple(ends))

    @property
    def header(self) -> nibabel.nifti1.Nifti1Header:
        """Return the header of the nifti1 image.

        Returns
        -------
        nibabel.nifti1.Nifti1Header
        """
        return self._image.header

    @property
    def offset(self) -> tuple:
        """Returns the lowest cartesian coordinates of the voxel data.

        Returns
        -------
        tuple
        """
        offset = np.array([self._image.header['qoffset_x'],
                           self._image.header['qoffset_y'],
                           self._image.header['qoffset_z']
                           ], dtype=np.float64)
        return offset * self._scaling

    @property
    def voxel_size(self) -> tuple:
        """Returns the sizes of a voxel in x-, y- and z-direction.

        Returns
        -------
        tuple
        """
        x, y, z = self._image.header.get_zooms()[:self._N_DIMENSION]
        return tuple(np.array((x, y, z), dtype=np.float64) * self._scaling)

    @property
    def xyz_shape(self) -> tuple:
        """Returns the number of voxels in x-, y- and z-direction.

        Returns
        -------
        tuple
        """
        return self._image.header.get_data_shape()[:self._N_DIMENSION]

    @staticmethod
    def _load_image(file_path: str) -> nibabel.nifti1.Nifti1Image:
        try:
            return nibabel.load(file_path)
        except FileNotFoundError:
            raise IOError('File Not Found.')

    @property
    def _scaling(self) -> float:
        xyz_unit = self._image.header.get_xyzt_units()[0]
        return {'unknown': 1.0,
                'meter': 1.0e3,
                'mm': 1.0,
                'micron': 1.0e-3}[xyz_unit]


class MagneticResonanceImage(Nifti1Image):

    def __init__(self, file_path: str) -> None:
        super().__init__(file_path)
        if not self.data_map.ndim == 3:
            raise IOError('MRI Data shape is not three dimensional.')


class DiffusionTensorImage(Nifti1Image):

    def __init__(self, file_path: str) -> None:
        super().__init__(file_path)
        if not self._image.get_fdata().ndim == 4:
            raise IOError('Data Shape not four dimensional.')
        if not self._image.get_fdata().shape[-1] == 6:
            raise IOError('Data Shape is not (x,y,z,6).')

    def diffusion(self) -> np.ndarray:
        x_max, y_max, z_max = self.xyz_shape
        diffusion = self.data_map.reshape((x_max * y_max * z_max, 6))
        return np.array([(xx, xy, xz, xy, yy, yz, xz, yz, zz)
                         for xx, xy, xz, yy, yz, zz in diffusion]
                        ).reshape(x_max, y_max, z_max, 3, 3).astype(float)
