
from src.brain_model import BrainModel
from src.voxels import Voxels
from src.brainsubstance import Material
import numpy as np
import pytest


class MockMRI():

    def __init__(self) -> None:
        self.__data = np.array([3, 2, 1, 3, 2, 1, 0, 0]).reshape((2, 2, 2))

    def bounding_box(self) -> np.ndarray:
        return (0, 0, 0), (1, 1, 1)

    def material_distribution(self, material: Material) -> Voxels:
        data = self.__data == material
        return Voxels(start=(0, 0, 0), end=(1, 1, 1), data=data)

    def xyz_shape(self) -> tuple:
        return self.__data.shape


@pytest.mark.skip
class TestDielectricDistribution:

    def test_complex_conductivity(self):
        model = BrainModel(mri=MockMRI())
        actual = model.conductivity(frequency=0)
        data = np.array([[[2.0 + 0j, 0.02 + 0j], [0.02 + 0j, 2.0 + 0j]],
                         [[0.02 + 0j, 0.02 + 0j], [0.02 + 0j, 0.02 + 0j]]])
        desired = Voxels(start=(0, 0, 0), end=(1, 1, 1), data=data)
        np.testing.assert_equal(actual.data, desired.data)
        assert actual.start == desired.start and actual.end == desired.end
