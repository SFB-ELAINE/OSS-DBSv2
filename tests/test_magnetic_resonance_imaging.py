from src.magnetic_resonance_imaging import MagneticResonanceImage, MaterialMap
from src.magnetic_resonance_imaging import BrainSubstance       
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

    def test_material_at_initial_mapping(self, mri_path):
        mri = MagneticResonanceImage(file_path=mri_path)
        positions = np.array([(0, 0, 0.5, ), 
                              (0, 0, 0),
                              (0.5, 0.4, 0), 
                              (0.2, 0.6, 0.4)])
        expected = [BrainSubstance.UNKNOWN,
                    BrainSubstance.CEREBROSPINAL_FLUID,
                    BrainSubstance.WHITE_MATTER,
                    BrainSubstance.GRAY_MATTER]
        assert mri.material_at(positions) == expected

    def test_material_at_custom_mapping(self, mri_path):
        material_map = MaterialMap(csf_index=3, wm_index=2, gm_index=1)
        mri = MagneticResonanceImage(file_path=mri_path,
                                     material_map=material_map)
        positions = [(0, 0, 0.5),
                     (0, 0, 0),
                     (0.5, 0.4, 0), 
                     (0.2, 0.6, 0.4)]
        expected = [BrainSubstance.UNKNOWN,
                    BrainSubstance.GRAY_MATTER,
                    BrainSubstance.WHITE_MATTER,
                    BrainSubstance.CEREBROSPINAL_FLUID]
        assert mri.material_at(positions) == expected