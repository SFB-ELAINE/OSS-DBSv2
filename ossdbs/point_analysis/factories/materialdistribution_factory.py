

from ossdbs.nifti1ImageX import Nifti1Image
from ossdbs.bounding_box import BoundingBox
from ossdbs.point_analysis.material_disribution import MaterialDistribution
import numpy as np


class MaterialDistributionFactory:

    def __init__(self,
                 nifti: Nifti1Image,
                 bounding_box: BoundingBox
                 ) -> None:

        self.__nifti = nifti
        self.__bbox = bounding_box

    def create(self) -> MaterialDistribution:

        self.__check_mri_data_shape()
        voxel_size = self.__nifti.voxel_size()
        offset = self.__nifti.offset()
        bbox = self.__bbox.intersection(self.__nifti.bounding_box())
        start, end = bbox.start, bbox.end
        start_index = np.floor(np.subtract(start, offset) / voxel_size)
        x_s, y_s, z_s = start_index.astype(int)
        end_index = np.ceil(np.subtract(end, offset) / voxel_size)
        x_e, y_e, z_e = end_index.astype(int)
        data = self.__nifti.data_map()[x_s:x_e, y_s:y_e, z_s:z_e]
        return MaterialDistribution(data, offset, voxel_size)

    def __check_mri_data_shape(self):
        if not self.__nifti.data_map().ndim == 3:
            raise IOError('MRI Data shape is not three dimensional.')
