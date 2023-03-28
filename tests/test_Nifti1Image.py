
from ossdbs.nifti1Image import Nifti1Image
import numpy as np
import pytest
import nibabel

from ossdbs.bounding_box import BoundingBox


class TestNifti1Image:
    @pytest.fixture
    def nifti1_image_3d(self, tmpdir):
        path = tmpdir.mkdir("Test_Nifti1Image").join("test_Nifti1Image_3d.nii")
        data = np.array([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]])
        affine = np.array([[0.5, 0, 0, 2],
                           [0, 0.5, 0, 2],
                           [0, 0, 0.5, 2],
                           [0, 0, 0, 1]])
        nii_image = nibabel.Nifti1Image(dataobj=data, affine=affine)
        nibabel.save(nii_image, path)
        return Nifti1Image(file_path=path)

    @pytest.fixture
    def nifti1_image_4d(self, tmpdir):
        path = tmpdir.mkdir("Test_Nifti1Image").join("test_Nifti1Image_4d.nii")
        data = np.array([[[[1.0], [2.0]], [[3.0], [4.0]]],
                        [[[5.0], [6.0]], [[7.0], [8.0]]]])
        affine = np.array([[0.5, 0, 0, 2],
                           [0, 0.5, 0, 2],
                           [0, 0, 0.5, 2],
                           [0, 0, 0, 1]])
        nii_image = nibabel.Nifti1Image(dataobj=data, affine=affine)
        nibabel.save(nii_image, path)
        return Nifti1Image(file_path=path)

    def test_inavlid_file_path(self):
        with pytest.raises(IOError):
            Nifti1Image(file_path='inavlid_file_path')

    def test_data_file(self, nifti1_image_3d):
        desired = np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]])
        np.testing.assert_equal(nifti1_image_3d.data_map(), desired)

    def test_boundingbox_units_default(self, nifti1_image_3d):
        desired = BoundingBox((2, 2, 2), (3, 3, 3))
        np.testing.assert_equal(nifti1_image_3d.bounding_box(), desired)

    def test_boundingbox_4d_shape(self, nifti1_image_4d):
        desired = BoundingBox((2, 2, 2), (3, 3, 3))
        np.testing.assert_equal(nifti1_image_4d.bounding_box(), desired)

    def test_boundingbox_units_mm(self, nifti1_image_3d):
        nifti1_image_3d.header()['xyzt_units'] = 2
        desired = BoundingBox((2, 2, 2), (3, 3, 3))
        np.testing.assert_equal(nifti1_image_3d.bounding_box(), desired)

    def test_boundingbox_units_meter(self, nifti1_image_3d):
        nifti1_image_3d.header()['xyzt_units'] = 1
        desired = BoundingBox((2e3, 2e3, 2e3), (3e3, 3e3, 3e3))
        np.testing.assert_equal(nifti1_image_3d.bounding_box(), desired)

    def test_boundingbox_units_micron(self, nifti1_image_3d):
        nifti1_image_3d.header()['xyzt_units'] = 3
        desired = BoundingBox((2e-3, 2e-3, 2e-3), (3e-3, 3e-3, 3e-3))
        np.testing.assert_equal(nifti1_image_3d.bounding_box(), desired)

    def test_voxel_size(self, nifti1_image_3d):
        desired = (0.5, 0.5, 0.5)
        np.testing.assert_equal(nifti1_image_3d.voxel_size(), desired)
