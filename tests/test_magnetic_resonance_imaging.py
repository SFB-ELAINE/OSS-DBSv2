
import nibabel
import pytest
import numpy as np

from ossdbs.brain_imaging import MagneticResonanceImage
from ossdbs.materials import Material
from ossdbs.voxels import Voxels


class TestMagneticResonanceImage:

    @pytest.fixture
    def mri(self, tmpdir):
        path = tmpdir.mkdir("Test_MRI").join("test_mri.nii")
        data = np.array([[[3, 0], [1, 2]], [[2, 3], [3, 1]]])
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
        desired = Voxels(start=(.002, .002, .002),
                         end=(.003, .003, .003),
                         data=is_material)
        actual = mri.material_distribution(Material.CSF)
        self.is_equal(actual, desired)

    def test_material_distribution_WhiteMatter(self, mri):
        is_material = np.array([[[False, False], [False, True]],
                                [[True, False], [False, False]]])
        desired = Voxels(start=(.002, .002, .002),
                         end=(.003, .003, .003),
                         data=is_material)
        actual = mri.material_distribution(Material.WHITE_MATTER)
        self.is_equal(actual, desired)

    def test_material_distribution_GrayMatter(self, mri):
        is_material = np.array([[[False, False], [True, False]],
                                [[False, False], [False, True]]])
        desired = Voxels(start=(.002, .002, .002),
                         end=(.003, .003, .003),
                         data=is_material)
        actual = mri.material_distribution(Material.GRAY_MATTER)
        self.is_equal(actual, desired)

    def test_material_distribution_Unknown(self, mri):
        is_material = np.array([[[False, True], [False, False]],
                                [[False, False], [False, False]]])
        desired = Voxels(start=(.002, .002, .002),
                         end=(.003, .003, .003),
                         data=is_material)
        actual = mri.material_distribution(Material.UNKNOWN)
        self.is_equal(actual, desired)


class TestMagneticResonanceImageFalseDataShape:

    def test_data_false_dimension(self, tmpdir):
        path = tmpdir.mkdir("Test_DTI").join("test_dti.nii")
        data = np.arange(0, 27, dtype=np.float64).reshape((3, 3, 3, 1))
        affine = np.array([[0.5, 0, 0, 2],
                           [0, 0.5, 0, 2],
                           [0, 0, 0.5, 2],
                           [0, 0, 0, 1]])
        nii_image = nibabel.Nifti1Image(dataobj=data, affine=affine)
        nibabel.save(nii_image, path)

        with pytest.raises(IOError):
            MagneticResonanceImage(file_path=path)

    def test_data_false_shape(self, tmpdir):
        path = tmpdir.mkdir("Test_DTI").join("test_dti.nii")
        data = np.arange(0, 9, dtype=np.float64).reshape((3, 3))
        affine = np.array([[0.5, 0, 0, 2],
                           [0, 0.5, 0, 2],
                           [0, 0, 0.5, 2],
                           [0, 0, 0, 1]])
        nii_image = nibabel.Nifti1Image(dataobj=data, affine=affine)
        nibabel.save(nii_image, path)

        with pytest.raises(IOError):
            MagneticResonanceImage(file_path=path)
