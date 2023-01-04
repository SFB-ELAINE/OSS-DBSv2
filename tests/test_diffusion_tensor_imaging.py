from ossdbs.brain_imaging.dti import DiffusionTensorImage
import nibabel
import pytest
import numpy as np


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
