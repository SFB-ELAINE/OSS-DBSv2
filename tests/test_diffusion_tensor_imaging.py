from src.diffusion_tensor_imaging import DiffusionTensorImage
import nibabel
import pytest
import numpy as np


class TestDiffusionTensorImage:

    @pytest.fixture
    def dti_path(self, tmpdir):
        path = tmpdir.mkdir("Test_DTI").join("test_dti.nii")
        data = np.array([[[1, 2, 3, 4, 5, 6], [7, 8, 9, 10, 11, 12]]])
        affine = np.array([[0.5, 0, 0, 2], 
                            [0, 0.5, 0, 2], 
                            [0, 0, 0.5, 2], 
                            [0, 0, 0, 1]])
        header = nibabel.Nifti1Header()
        nii_image = nibabel.Nifti1Image(dataobj=data, 
                                        affine=affine, 
                                        header=header)
        nibabel.save(nii_image, path)
        return path
