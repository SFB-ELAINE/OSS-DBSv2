

from ossdbs.bounding_box import BoundingBox
from ossdbs.Nifti1Image import Nifti1Image
from ossdbs.encapsulatin_layer import EncapsulatingLayers
from ossdbs.materials import Material
from ossdbs.conductivity import Conductivity
import numpy as np


class ConductivityFactory:
    """Transform the input json.

    Parameters
    ----------
    json_path : str
    """

    def __init__(self,
                 nifti: Nifti1Image,
                 bounding_box: BoundingBox,
                 encap_layer: EncapsulatingLayers,
                 encap_material: Material = Material.GRAY_MATTER) -> None:
        self.__nifti = nifti
        self.__bbox = bounding_box
        self.__encap_layer = encap_layer
        self.__encap_material = encap_material

    def create(self) -> Conductivity:
        """Return the conductivity.

        Returns
        -------
        Conductivity
            Conductivity distribution in a given space.
        """

        self.__check_mri_data_shape()
        voxel_size = self.__nifti.voxel_size()
        offset = self.__nifti.offset()
        start, end = self.__bbox.start, self.__bbox.end
        start_index = np.floor(np.subtract(start, offset) / voxel_size)
        x_s, y_s, z_s = start_index.astype(int)
        end_index = np.ceil(np.subtract(end, offset) / voxel_size)
        x_e, y_e, z_e = end_index.astype(int)

        data = self.__nifti.data_map()[x_s:x_e, y_s:y_e, z_s:z_e]
        new_start = tuple(start_index * voxel_size + offset)
        new_end = tuple(end_index * voxel_size + offset)
        self.__set_encap_data(data, new_start)

        return Conductivity(data, BoundingBox(start=new_start, end=new_end))

    def __set_encap_data(self, data: np.ndarray, offset: tuple) -> None:

        voxel_size = self.__nifti.voxel_size()
        bounding_boxes = [bbox.intersection(self.__bbox)
                          for bbox in self.__encap_layer.bounding_boxes()]
        points = np.concatenate([bbox.points(offset, voxel_size)
                                 for bbox in bounding_boxes])

        included = self.__encap_layer.is_included(points=points)
        encap_points = points[included]
        point_indices = (encap_points - offset) / self.__nifti.voxel_size()

        for index in point_indices:
            x, y, z = index.astype(int)
            data[x, y, z] = self.__encap_material

    def __check_mri_data_shape(self):
        if not self.__nifti.data_map().ndim == 3:
            raise IOError('MRI Data shape is not three dimensional.')
