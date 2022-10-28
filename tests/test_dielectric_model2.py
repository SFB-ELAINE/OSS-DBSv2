from src.dielectric_model import DielectricModel2
from src.brainsubstance import BrainSubstance
import numpy as np
import pytest


class TestWhiteMatterModel2:

    @pytest.fixture
    def model(self):
        return DielectricModel2.create_model(BrainSubstance.WHITE_MATTER)

    def test_relative_permitivity_0Hz(self, model):
        result = model.relative_permitivity(frequency=0)
        np.testing.assert_allclose(result, 36.999965, atol=0.1)

    def test_relative_permitivity_400MHz(self, model):
        result = model.relative_permitivity(frequency=400e6*2*np.pi)
        np.testing.assert_allclose(result, 36.314021, atol=0.1)

    def test_relative_permitivity_4GHz(self, model):
        result = model.relative_permitivity(frequency=400e9*np.pi)
        np.testing.assert_allclose(result, 6.659155, atol=0.1)

    def test_conductivity_0Hz(self, model):
        result = model.conductivity(frequency=0)
        np.testing.assert_allclose(result, 0.47, atol=0.001)

    def test_conductivity_400MHz(self, model):
        result = model.conductivity(frequency=400e6*2*np.pi)
        np.testing.assert_allclose(result, 0.50382, atol=0.001)

    def test_conductivity_4GHz(self, model):
        result = model.conductivity(frequency=400e9**2*np.pi)
        np.testing.assert_allclose(result, 35839.099752, atol=0.001)


class TestGrayMatterModel2:

    @pytest.fixture
    def model(self):
        return DielectricModel2.create_model(BrainSubstance.GRAY_MATTER)

    def test_relative_permitivity_0Hz(self, model):
        result = model.relative_permitivity(frequency=0)
        np.testing.assert_allclose(result, 55.5, atol=0.1)

    def test_relative_permitivity_400MHz(self, model):
        result = model.relative_permitivity(frequency=400e6*2*np.pi)
        np.testing.assert_allclose(result, 55.152151, atol=0.1)

    def test_relative_permitivity_400GHz(self, model):
        result = model.relative_permitivity(frequency=400e9*np.pi)
        np.testing.assert_allclose(result, 6.092419, atol=0.1)

    def test_conductivity_0Hz(self, model):
        result = model.conductivity(frequency=0)
        np.testing.assert_allclose(result, 1.03, atol=0.001)

    def test_conductivity_400MHz(self, model):
        result = model.conductivity(frequency=400e6*2*np.pi)
        np.testing.assert_allclose(result, 1.064773, atol=0.001)

    def test_conductivity_400GHz(self, model):
        result = model.conductivity(frequency=400e9**2*np.pi)
        np.testing.assert_allclose(result, 1872.869183, atol=0.001)


class TestCerebrospinalFluidModel2:

    @pytest.fixture
    def model(self):
        return DielectricModel2.create_model(
                                            BrainSubstance.CEREBROSPINAL_FLUID)

    def test_relative_permitivity_0Hz(self, model):
        result = model.relative_permitivity(frequency=0)
        np.testing.assert_allclose(result, 80, atol=1)

    def test_relative_permitivity_1Hz(self, model):
        result = model.relative_permitivity(frequency=1*2*np.pi)
        np.testing.assert_allclose(result, 80, atol=1)

    def test_relative_permitivity_1kHz(self, model):
        result = model.relative_permitivity(frequency=1000*2*np.pi)
        np.testing.assert_allclose(result, 80, atol=1)

    def test_conductivity_0Hz(self, model):
        result = model.conductivity(frequency=0)
        np.testing.assert_allclose(result, 1.79, atol=0.001)

    def test_conductivity_1Hz(self, model):
        result = model.conductivity(frequency=1*2*np.pi)
        np.testing.assert_allclose(result, 1.79, atol=0.001)

    def test_conductivity_1kHz(self, model):
        result = model.conductivity(frequency=1000*2*np.pi)
        np.testing.assert_allclose(result, 1.79, atol=0.001)
