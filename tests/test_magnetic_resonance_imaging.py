
import nibabel
import pytest
import numpy as np

from src.brain_imaging import MagneticResonanceImage
from src.brainsubstance import Material
from src.voxels import Voxels


class TestMagneticResonanceImage:

    @pytest.fixture
    def mri(self, tmpdir):
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
        return MagneticResonanceImage(file_path=path)

    def is_equal(self, actual, desired):
        np.testing.assert_equal(actual.data, desired.data)
        assert actual.start == desired.start and actual.end == desired.end

    def test_material_distribution_CSF(self, mri):
        is_material = np.array([[[True, False], [False, False]],
                                [[False, True], [True, False]]])
        desired = Voxels(start=(2, 2, 2), end=(3, 3, 3), data=is_material)
        actual = mri.material_distribution(Material.CSF)
        self.is_equal(actual, desired)

    def test_material_distribution_WhiteMatter(self, mri):
        is_material = np.array([[[False, False], [False, True]],
                                [[True, False], [False, False]]])
        desired = Voxels(start=(2, 2, 2), end=(3, 3, 3), data=is_material)
        actual = mri.material_distribution(Material.WHITE_MATTER)
        self.is_equal(actual, desired)

    def test_material_distribution_GrayMatter(self, mri):
        is_material = np.array([[[False, False], [True, False]],
                                [[False, False], [False, True]]])
        desired = Voxels(start=(2, 2, 2), end=(3, 3, 3), data=is_material)
        actual = mri.material_distribution(Material.GRAY_MATTER)
        self.is_equal(actual, desired)

    def test_material_distribution_Unknown(self, mri):
        is_material = np.array([[[False, True], [False, False]],
                                [[False, False], [False, False]]])
        desired = Voxels(start=(2, 2, 2), end=(3, 3, 3), data=is_material)
        actual = mri.material_distribution(Material.UNKNOWN)
        self.is_equal(actual, desired)
