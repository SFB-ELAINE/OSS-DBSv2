
import nibabel
import pytest
import numpy as np


class TestMagneticResonanceImage:

    @pytest.fixture
    def mri_path(self, tmpdir):
        path = tmpdir.mkdir("Test_MRI").join("test_mri.nii")
        data = np.array([[[1, 0], [3, 2]], [[2, 1], [1, 3]]])
        affine = np.array([[0.5, 0, 0, 2],
                           [0, 0.5, 0, 2],
                           [0, 0, 0.5, 2],
                           [0, 0, 0, 1]])
        header = nibabel.Nifti1Header()
        header['xyzt_units'] = 2
        nii_image = nibabel.Nifti1Image(dataobj=data,
                                        affine=affine,
                                        header=header)
        nibabel.save(nii_image, path)
        return path
