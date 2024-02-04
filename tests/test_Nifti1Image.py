import nibabel
import numpy as np
import pytest

from ossdbs.utils.nifti1image import (
    DiffusionTensorImage,
    MagneticResonanceImage,
    Nifti1Image,
)


class TestNifti1Image:
    @pytest.fixture
    def nifti1_image_3d(self, tmpdir):
        path = tmpdir.mkdir("Test_Nifti1Image").join("test_Nifti1Image_3d.nii")
        data = np.array([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]])
        affine = np.array(
            [[0.5, 0, 0, 2], [0, 0.5, 0, 2], [0, 0, 0.5, 2], [0, 0, 0, 1]]
        )
        nii_image = nibabel.Nifti1Image(dataobj=data, affine=affine)
        nibabel.save(nii_image, path)
        return Nifti1Image(file_path=path)

    @pytest.fixture
    def nifti1_image_4d(self, tmpdir):
        path = tmpdir.mkdir("Test_Nifti1Image").join("test_Nifti1Image_4d.nii")
        data = np.array(
            [[[[1.0], [2.0]], [[3.0], [4.0]]], [[[5.0], [6.0]], [[7.0], [8.0]]]]
        )
        affine = np.array(
            [[0.5, 0, 0, 2], [0, 0.5, 0, 2], [0, 0, 0.5, 2], [0, 0, 0, 1]]
        )
        nii_image = nibabel.Nifti1Image(dataobj=data, affine=affine)
        nibabel.save(nii_image, path)
        return Nifti1Image(file_path=path)

    def test_invalid_file_path(self):
        with pytest.raises(FileNotFoundError):
            Nifti1Image(file_path="invalid_file_path")

    def test_data_file(self, nifti1_image_3d):
        desired = np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]])
        np.testing.assert_equal(nifti1_image_3d.data, desired)


class TestMagneticResonanceImageFalseShape:
    @pytest.fixture
    def mri(self, tmpdir):
        path = tmpdir.mkdir("Test_MRIImage").join("test_Nifti1Image_4d.nii")
        data = np.array(
            [[[[1.0], [2.0]], [[3.0], [4.0]]], [[[5.0], [6.0]], [[7.0], [8.0]]]]
        )
        affine = np.array(
            [[0.5, 0, 0, 2], [0, 0.5, 0, 2], [0, 0, 0.5, 2], [0, 0, 0, 1]]
        )
        nii_image = nibabel.Nifti1Image(dataobj=data, affine=affine)
        nibabel.save(nii_image, path)

        with pytest.raises(IOError):
            MagneticResonanceImage(file_path=path)


class TestDiffusionTensorImage:
    @pytest.fixture
    def dti(self, tmpdir):
        path = tmpdir.mkdir("Test_DTI").join("test_dti.nii")
        data = np.arange(1, 25, dtype=np.float64).reshape((1, 2, 2, 6))
        affine = np.array(
            [[0.5, 0, 0, 2], [0, 0.5, 0, 2], [0, 0, 0.5, 2], [0, 0, 0, 1]]
        )
        nii_image = nibabel.Nifti1Image(dataobj=data, affine=affine)
        nibabel.save(nii_image, path)
        return DiffusionTensorImage(file_path=path)


class TestDiffusionTensorImageFalseShape:
    def test_data_false_dimension(self, tmpdir):
        path = tmpdir.mkdir("Test_DTI").join("test_dti.nii")
        data = np.arange(0, 27, dtype=np.float64).reshape((3, 3, 3))
        affine = np.array(
            [[0.5, 0, 0, 2], [0, 0.5, 0, 2], [0, 0, 0.5, 2], [0, 0, 0, 1]]
        )
        nii_image = nibabel.Nifti1Image(dataobj=data, affine=affine)
        nibabel.save(nii_image, path)

        with pytest.raises(IOError):
            DiffusionTensorImage(file_path=path)

    def test_data_false_shape(self, tmpdir):
        path = tmpdir.mkdir("Test_DTI").join("test_dti.nii")
        data = np.arange(0, 27, dtype=np.float64).reshape((3, 3, 3, 1))
        affine = np.array(
            [[0.5, 0, 0, 2], [0, 0.5, 0, 2], [0, 0, 0.5, 2], [0, 0, 0, 1]]
        )
        nii_image = nibabel.Nifti1Image(dataobj=data, affine=affine)
        nibabel.save(nii_image, path)

        with pytest.raises(IOError):
            DiffusionTensorImage(file_path=path)
