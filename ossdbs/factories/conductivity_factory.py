

from ossdbs.dielectric_model import DielectricModel
from ossdbs.bounding_box import BoundingBox
from ossdbs.nifti1Image import Nifti1Image
from ossdbs.electrodes import Electrodes
from ossdbs.materials import Material
from ossdbs.conductivity import Conductivity
import numpy as np


class ConductivityFactory:
    """Transform the input json.

    Parameters
    ----------
    json_path : str
    """

    __MATERIALS = {'Blood': Material.BLOOD,
                   'CSF': Material.CSF,
                   'GrayMatter': Material.GRAY_MATTER,
                   'WhiteMAtter': Material.WHITE_MATTER,
                   'Unknown': Material.UNKNOWN}

    def __init__(self,
                 nifti: Nifti1Image,
                 bounding_box: BoundingBox,
                 dielectrc_model: DielectricModel,
                 electrodes: Electrodes,
                 encap_material: str = 'GrayMatter',
                 ) -> None:
        self.__nifti = nifti
        self.__bbox = bounding_box
        self.__electrodes = electrodes
        self.__encap_material = self.__MATERIALS[encap_material]
        self.__dielectric_model = dielectrc_model

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
        bbox = self.__bbox.intersection(self.__nifti.bounding_box())
        start, end = bbox.start, bbox.end
        start_index = np.floor(np.subtract(start, offset) / voxel_size)
        x_s, y_s, z_s = start_index.astype(int)
        end_index = np.ceil(np.subtract(end, offset) / voxel_size)
        x_e, y_e, z_e = end_index.astype(int)

        data = self.__nifti.data_map()[x_s:x_e, y_s:y_e, z_s:z_e]
        new_start = tuple(start_index * voxel_size + offset)
        new_end = tuple(end_index * voxel_size + offset)
        # TODO
        # self.__set_encapapsulation(data, new_start)

        bounding_box = BoundingBox(start=new_start, end=new_end)

        return Conductivity(data, bounding_box, self.__dielectric_model)

    def __set_encapapsulation(self, data: np.ndarray, offset: tuple) -> None:
        """
        TODO
        """

        encapsulation = self.__electrodes.encapsulation()
        voxel_size = self.__nifti.voxel_size()
        bounding_boxes = [bbox.intersection(self.__bbox)
                          for bbox in encapsulation.bounding_boxes()]
        points = np.concatenate([bbox.points(offset, voxel_size)
                                 for bbox in bounding_boxes])

        included = encapsulation.is_included(points=points)
        encap_points = points[included]
        point_indices = (encap_points - offset) / self.__nifti.voxel_size()

        for index in point_indices:
            x, y, z = index.astype(int)
            data[x, y, z] = self.__encap_material

    def __check_mri_data_shape(self):
        if not self.__nifti.data_map().ndim == 3:
            raise IOError('MRI Data shape is not three dimensional.')
