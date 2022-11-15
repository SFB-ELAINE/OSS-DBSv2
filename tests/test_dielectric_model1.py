from ossdbs.dielectric_model import DielectricModel1
from ossdbs.brainsubstance import Material
import numpy as np
import pytest


class TestWhiteMatterModel1:

    @pytest.fixture
    def model(self):
        return DielectricModel1.create_model(Material.WHITE_MATTER)

    def test_relative_permitivity_0Hz(self, model):
        result = model.relative_permitivity(frequency=0)
        np.testing.assert_allclose(result, 35040136, atol=1)

    def test_relative_permitivity_1Hz(self, model):
        result = model.relative_permitivity(frequency=1*2*np.pi)
        np.testing.assert_allclose(result, 34883924, atol=1)

    def test_relative_permitivity_1kHz(self, model):
        result = model.relative_permitivity(frequency=1000*2*np.pi)
        np.testing.assert_allclose(result, 69710, atol=1)

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

    def test_relative_permitivity_0Hz(self, model):
        result = model.relative_permitivity(frequency=0)
        np.testing.assert_allclose(result, 45200449, atol=1)

    def test_relative_permitivity_1Hz(self, model):
        result = model.relative_permitivity(frequency=1*2*np.pi)
        np.testing.assert_allclose(result, 45150280, atol=1)

    def test_relative_permitivity_1kHz(self, model):
        result = model.relative_permitivity(frequency=1000*2*np.pi)
        np.testing.assert_allclose(result, 164064, atol=1)

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

    def test_relative_permitivity_0Hz(self, model):
        result = model.relative_permitivity(frequency=0)
        np.testing.assert_allclose(result, 109, atol=1)

    def test_relative_permitivity_1Hz(self, model):
        result = model.relative_permitivity(frequency=1*2*np.pi)
        np.testing.assert_allclose(result, 109, atol=1)

    def test_relative_permitivity_1kHz(self, model):
        result = model.relative_permitivity(frequency=1000*2*np.pi)
        np.testing.assert_allclose(result, 109, atol=1)

    def test_conductivity_0Hz(self, model):
        result = model.conductivity(frequency=0)
        np.testing.assert_allclose(result, 2.00, atol=0.001)

    def test_conductivity_1Hz(self, model):
        result = model.conductivity(frequency=1*2*np.pi)
        np.testing.assert_allclose(result, 2.00, atol=0.001)

    def test_conductivity_1kHz(self, model):
        result = model.conductivity(frequency=1000*2*np.pi)
        np.testing.assert_allclose(result, 2.00, atol=0.001)
