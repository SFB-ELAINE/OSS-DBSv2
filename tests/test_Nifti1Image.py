
from ossdbs.utils.nifti1image import (Nifti1Image,
                                      MagneticResonanceImage,
                                      DiffusionTensorImage)
import numpy as np
import pytest
import nibabel

from ossdbs.model_geometry.bounding_box import BoundingBox


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

    def test_invalid_file_path(self):
        with pytest.raises(IOError):
            Nifti1Image(file_path='inavlid_file_path')

    def test_data_file(self, nifti1_image_3d):
        desired = np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]])
        np.testing.assert_equal(nifti1_image_3d.data_map, desired)

    def test_boundingbox_units_default(self, nifti1_image_3d):
        desired = BoundingBox((2, 2, 2), (3, 3, 3))
        np.testing.assert_equal(nifti1_image_3d.bounding_box, desired)

    def test_boundingbox_4d_shape(self, nifti1_image_4d):
        desired = BoundingBox((2, 2, 2), (3, 3, 3))
        np.testing.assert_equal(nifti1_image_4d.bounding_box, desired)

    def test_boundingbox_units_mm(self, nifti1_image_3d):
        nifti1_image_3d.header['xyzt_units'] = 2
        desired = BoundingBox((2, 2, 2), (3, 3, 3))
        np.testing.assert_equal(nifti1_image_3d.bounding_box, desired)

    def test_boundingbox_units_meter(self, nifti1_image_3d):
        nifti1_image_3d.header['xyzt_units'] = 1
        desired = BoundingBox((2e3, 2e3, 2e3), (3e3, 3e3, 3e3))
        np.testing.assert_equal(nifti1_image_3d.bounding_box, desired)

    def test_boundingbox_units_micron(self, nifti1_image_3d):
        nifti1_image_3d.header['xyzt_units'] = 3
        desired = BoundingBox((2e-3, 2e-3, 2e-3), (3e-3, 3e-3, 3e-3))
        np.testing.assert_equal(nifti1_image_3d.bounding_box, desired)

    def test_voxel_size(self, nifti1_image_3d):
        desired = (0.5, 0.5, 0.5)
        np.testing.assert_equal(nifti1_image_3d.voxel_size, desired)


class TestMagneticResonanceImageFalseShape:

    @pytest.fixture
    def mri(self, tmpdir):
        path = tmpdir.mkdir("Test_MRIImage").join("test_Nifti1Image_4d.nii")
        data = np.array([[[[1.0], [2.0]], [[3.0], [4.0]]],
                        [[[5.0], [6.0]], [[7.0], [8.0]]]])
        affine = np.array([[0.5, 0, 0, 2],
                           [0, 0.5, 0, 2],
                           [0, 0, 0.5, 2],
                           [0, 0, 0, 1]])
        nii_image = nibabel.Nifti1Image(dataobj=data, affine=affine)
        nibabel.save(nii_image, path)

        with pytest.raises(IOError):
            MagneticResonanceImage(file_path=path)


class TestDiffusionTensorImage:

    @pytest.fixture
    def dti(self, tmpdir):
        path = tmpdir.mkdir("Test_DTI").join("test_dti.nii")
        data = np.arange(1, 25, dtype=np.float64).reshape((1, 2, 2, 6))
        affine = np.array([[0.5, 0, 0, 2],
                           [0, 0.5, 0, 2],
                           [0, 0, 0.5, 2],
                           [0, 0, 0, 1]])
        nii_image = nibabel.Nifti1Image(dataobj=data, affine=affine)
        nibabel.save(nii_image, path)
        return DiffusionTensorImage(file_path=path)

    def test_diffusion(self, dti):
        order = np.array((1, 2, 3, 2, 4, 5, 3, 5, 6))
        desired = np.array([[[order, order + 6],
                             [order + 12, order + 18]]
                            ], dtype=float).reshape((1, 2, 2, 3, 3))
        np.testing.assert_equal(dti.diffusion(), desired)


class TestDiffusionTensorImageFalseShape:

    def test_data_false_dimension(self, tmpdir):
        path = tmpdir.mkdir("Test_DTI").join("test_dti.nii")
        data = np.arange(0, 27, dtype=np.float64).reshape((3, 3, 3))
        affine = np.array([[0.5, 0, 0, 2],
                           [0, 0.5, 0, 2],
                           [0, 0, 0.5, 2],
                           [0, 0, 0, 1]])
        nii_image = nibabel.Nifti1Image(dataobj=data, affine=affine)
        nibabel.save(nii_image, path)

        with pytest.raises(IOError):
            DiffusionTensorImage(file_path=path)

    def test_data_false_shape(self, tmpdir):
        path = tmpdir.mkdir("Test_DTI").join("test_dti.nii")
        data = np.arange(0, 27, dtype=np.float64).reshape((3, 3, 3, 1))
        affine = np.array([[0.5, 0, 0, 2],
                           [0, 0.5, 0, 2],
                           [0, 0, 0.5, 2],
                           [0, 0, 0, 1]])
        nii_image = nibabel.Nifti1Image(dataobj=data, affine=affine)
        nibabel.save(nii_image, path)

        with pytest.raises(IOError):
            DiffusionTensorImage(file_path=path)
