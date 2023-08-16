from ossdbs.model_geometry import BoundingBox
import nibabel
import numpy as np
import numpy.linalg as npl
import netgen.occ
import ngsolve


class Nifti1Image:
    """Interface for Nifti1 image.

    Attributes
    ----------
    file_path : str
        File path of a Nifti1 image.

    Notes
    -----

    TODO: add link to Nibabel documentation

    """

    def __init__(self, file_path: str) -> None:
        self._image = self._load_image(file_path)
        self._trafo_cf = None

    @property
    def affine(self):
        return self._scaling * self._image.affine

    @property
    def data(self) -> np.ndarray:
        """Return the data of the nifti1 image.

        Returns
        -------
        np.memmap
        """
        return self._image.get_fdata()

    @property
    def bounding_box(self) -> BoundingBox:
        """Return the bounding box of the voxel data in voxel space.

        Returns
        -------
        BoundingBox
        """
        start = (0, 0, 0)
        ends = self.xyz_shape
        return BoundingBox(start, ends)

    @property
    def header(self) -> nibabel.nifti1.Nifti1Header:
        """Return the header of the nifti1 image.

        Returns
        -------
        nibabel.nifti1.Nifti1Header
        """
        return self._image.header

    @property
    def trafo_matrix(self) -> np.ndarray:
        return np.array([self.affine[0, :3],
                         self.affine[1, :3],
                         self.affine[2, :3]])

    @property
    def translation(self) -> np.ndarray:
        """Returns the lowest cartesian coordinates of the voxel data.

        Returns
        -------
        tuple
        """
        return self.affine[:3, 3]

    @property
    def xyz_shape(self) -> tuple:
        """Returns the number of voxels in x-, y- and z-direction.

        Returns
        -------
        tuple
        """
        return self._image.header.get_data_shape()

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

    def _crop_image(self, voxel_bounding_box: BoundingBox) -> np.ndarray:
        """Crop the image to match, for example, the brain geometry.

        Parameters:
            voxel_bounding_box:
                    Defined in voxel space!
        Returns
        -------
        TODO
        """

        start_index = np.floor(np.array(voxel_bounding_box.start))
        x_s, y_s, z_s = start_index.astype(int)
        end_index = np.floor(np.array(voxel_bounding_box.end))
        x_e, y_e, z_e = end_index.astype(int)
        bbox = BoundingBox((x_s, y_s, z_s), (x_e, y_e, z_e))

        return self.data[x_s:x_e, y_s:y_e, z_s:z_e], bbox

    def get_voxel_bounding_box(self, bounding_box: BoundingBox):
        inv_trafo = npl.inv(self.trafo_matrix)
        inv_affine_trafo = netgen.occ.gp_GTrsf(mat=inv_trafo.ravel(), vec=-inv_trafo.dot(self.translation))
        box = netgen.occ.Box(bounding_box.start, bounding_box.end)
        box_voxel = inv_affine_trafo(box)
        voxel_bounding_box = box_voxel.bounding_box
        start = voxel_bounding_box[0]
        end = voxel_bounding_box[1]
        # convert to float
        bbox_start = [start.x, start.y, start.z]
        bbox_end = [end.x, end.y, end.z]
        for i in range(3):
            if bbox_start[i] < 0:
                bbox_start[i] = 0
            if bbox_end[i] > self.xyz_shape[i]:
                bbox_end[i] = self.xyz_shape[i]
        return BoundingBox(np.floor(np.array(bbox_start)).astype(int),
                           np.floor(np.array(bbox_end)).astype(int))

    @property
    def trafo_cf(self):
        if self._trafo_cf is None:
            # ravel because ngsolve needs a vector in a tuple
            inv_trafo_matrix = tuple(npl.inv(self.trafo_matrix).ravel())
            mm_space_coordinates = ngsolve.CoefficientFunction((ngsolve.x, ngsolve.y, ngsolve.z), dims=(3,))
            mm_space_to_voxel_space = ngsolve.CoefficientFunction(inv_trafo_matrix, dims=(3, 3))
            # Getting image offset for translation
            translation = ngsolve.CoefficientFunction(tuple(self.translation), dims=(3,))
            self._trafo_cf = mm_space_to_voxel_space * (mm_space_coordinates - translation)
            self._trafo_cf.Compile()
        return self._trafo_cf


class MagneticResonanceImage(Nifti1Image):

    def __init__(self, file_path: str) -> None:
        super().__init__(file_path)
        if not self.data.ndim == 3:
            raise IOError('MRI Data shape is not three dimensional.')


class DiffusionTensorImage(Nifti1Image):

    def __init__(self, file_path: str) -> None:
        super().__init__(file_path)
        if not self._image.get_fdata().ndim == 4:
            raise IOError('Data Shape not four dimensional.')
        if not self._image.get_fdata().shape[-1] == 6:
            raise IOError('Data Shape is not (x,y,z,6).')

    @property
    def components(self):
        """Component and respective entry in flattened array
        """
        return {"xx": 0,
                "xy": 1,
                "xz": 2,
                "yx": 1,
                "yy": 3,
                "yz": 4,
                "zx": 2,
                "zy": 4,
                "zz": 5}

    def get_component(self, component: str) -> np.ndarray:
        """Get component of DTI tensor

        Notes
        -----

        Determine which index to extract based on data shape (x,y,z,6)
        Note the pairs (xy,yx), (yz,zy), and (xz,zx) map to the same indices
        """
        if component not in self.components:
            raise ValueError("Component must be in {}".format(self.components.keys()))
        extract_index = self.components[component]
        return self.data[:, :, :, extract_index]
