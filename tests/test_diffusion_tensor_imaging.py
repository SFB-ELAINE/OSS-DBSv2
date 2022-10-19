from src.diffusion_tensor_imaging import DiffusionTensorImage
import nibabel
import pytest
import numpy as np


class TestDiffusionTensorImage:

    @pytest.fixture
    def dti_path(self, tmpdir):
        path = tmpdir.mkdir("Test_DTI").join("test_dti.nii")
        data = np.array([[[(1, 2, 3, 4, 5, 6), 
                           (7, 8, 9, 10, 11, 12)],
                          [(13, 14, 15, 16, 17, 18), 
                           (19, 20, 21, 22, 23, 24)]],
                         [[(1, 2, 3, 4, 5, 6), 
                           (7, 8, 9, 10, 11, 12)],
                          [(13, 14, 15, 16, 17, 18), 
                           (19, 20, 21, 22, 23, 24)]],])
        affine = np.array([[0.5, 0, 0, 2], 
                            [0, 0.5, 0, 2], 
                            [0, 0, 0.5, 2], 
                            [0, 0, 0, 1]])
        nii_image = nibabel.Nifti1Image(dataobj=data, affine=affine)
        nibabel.save(nii_image, path)
        return path

    def test_diffusion_at(self, dti_path):
        dti = DiffusionTensorImage(file_path=dti_path)
        positions = [(0, 0, 0.5), 
                     (0, 0, 0),
                     (0.5, 0.4, 0), 
                     (0.2, 0.6, 0.4)]
        expected = [[(7, 8, 9), (8, 10, 11), (9, 11, 12)],
                    [(1, 2, 3), (2, 4, 5), (3, 5, 6)],
                    [(1, 2, 3), (2, 4, 5), (3, 5, 6)],
                    [(13, 14, 15), (14, 16, 17), (15, 17, 18)]]
        assert dti.diffusion_at(positions) == expected