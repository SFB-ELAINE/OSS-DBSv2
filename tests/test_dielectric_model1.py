from src.dielectric_model import DielectricModel1
from src.brainsubstance import Material
import numpy as np
import pytest


class TestWhiteMatterModel1:

    @pytest.fixture
    def model(self):
        return DielectricModel1.create_model(Material.WHITE_MATTER)

    def test_permitivity_0Hz(self, model):
        result = model.permitivity(frequency=0)
        np.testing.assert_allclose(result, 0.00031, atol=1e-5)

    def test_permitivity_1Hz(self, model):
        result = model.permitivity(frequency=1*2*np.pi)
        np.testing.assert_allclose(result, 0.000309, atol=1e-5)

    def test_permitivity_1kHz(self, model):
        result = model.permitivity(frequency=1000*2*np.pi)
        np.testing.assert_allclose(result, 6.172316e-07, atol=1e-8)

    def test_conductivity_0Hz(self, model):
        result = model.conductivity(frequency=0)
        np.testing.assert_allclose(result, 0.02, atol=0.001)

    def test_conductivity_1Hz(self, model):
        result = model.conductivity(frequency=1*2*np.pi)
        np.testing.assert_allclose(result, 0.020103, atol=0.001)

    def test_conductivity_1kHz(self, model):
        result = model.conductivity(frequency=1000*2*np.pi)
        np.testing.assert_allclose(result, 0.062575, atol=0.001)


class TestGrayMatterModel1:

    @pytest.fixture
    def model(self):
        return DielectricModel1.create_model(Material.GRAY_MATTER)

    def test_permitivity_0Hz(self, model):
        result = model.permitivity(frequency=0)
        np.testing.assert_allclose(result, 0.0004, atol=1e-5)

    def test_permitivity_1Hz(self, model):
        result = model.permitivity(frequency=1*2*np.pi)
        np.testing.assert_allclose(result, 0.0004, atol=1e-5)

    def test_permitivity_1kHz(self, model):
        result = model.permitivity(frequency=1000*2*np.pi)
        np.testing.assert_allclose(result, 1.452645e-06, atol=1e-7)

    def test_conductivity_0Hz(self, model):
        result = model.conductivity(frequency=0)
        np.testing.assert_allclose(result, 0.02, atol=0.001)

    def test_conductivity_1Hz(self, model):
        result = model.conductivity(frequency=1*2*np.pi)
        np.testing.assert_allclose(result, 0.020103, atol=0.001)

    def test_conductivity_1kHz(self, model):
        result = model.conductivity(frequency=1000*2*np.pi)
        np.testing.assert_allclose(result, 0.098807, atol=0.001)


class TestCerebrospinalFluidModel1:

    @pytest.fixture
    def model(self):
        return DielectricModel1.create_model(
                                        Material.CSF)

    def test_permitivity_0Hz(self, model):
        result = model.permitivity(frequency=0)
        np.testing.assert_allclose(result, 9.651065e-10, atol=1e-11)

    def test_permitivity_1Hz(self, model):
        result = model.permitivity(frequency=1*2*np.pi)
        np.testing.assert_allclose(result, 9.651065e-10, atol=1e-11)

    def test_permitivity_1kHz(self, model):
        result = model.permitivity(frequency=1000*2*np.pi)
        np.testing.assert_allclose(result, 9.651065e-10, atol=1e-11)

    def test_conductivity_0Hz(self, model):
        result = model.conductivity(frequency=0)
        np.testing.assert_allclose(result, 2.00, atol=0.001)

    def test_conductivity_1Hz(self, model):
        result = model.conductivity(frequency=1*2*np.pi)
        np.testing.assert_allclose(result, 2.00, atol=0.001)

    def test_conductivity_1kHz(self, model):
        result = model.conductivity(frequency=1000*2*np.pi)
        np.testing.assert_allclose(result, 2.00, atol=0.001)
