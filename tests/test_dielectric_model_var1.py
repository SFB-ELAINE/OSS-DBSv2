from src.dielectric_model.dielectric_model_var1 import DielectricModelVariant1
from src.brainsubstance import BrainSubstance
import numpy as np
import pytest


class TestWhiteMatterVariant1:

    @pytest.fixture
    def model(self):
        return DielectricModelVariant1.create_model(
                                                BrainSubstance.WHITE_MATTER)

    def test_permitivity_0Hz(self, model):
        result = model.permitivity(frequency=0)
        np.testing.assert_allclose(result, 35040136, atol=1)

    def test_permitivity_1Hz(self, model):
        result = model.permitivity(frequency=1*2*np.pi)
        np.testing.assert_allclose(result, 34884019, atol=1)

    def test_permitivity_1kHz(self, model):
        result = model.permitivity(frequency=1000*2*np.pi)
        np.testing.assert_allclose(result, 69810, atol=1)

    def test_conductivity_0Hz(self, model):
        result = model.conductivity(frequency=0)
        np.testing.assert_allclose(result, 0.02, atol=0.001)

    def test_conductivity_1Hz(self, model):
        result = model.conductivity(frequency=1*2*np.pi)
        np.testing.assert_allclose(result, 0.020103, atol=0.001)

    def test_conductivity_1kHz(self, model):
        result = model.conductivity(frequency=1000*2*np.pi)
        np.testing.assert_allclose(result, 0.062575, atol=0.001)


class TestGrayMatterVariant1:

    @pytest.fixture
    def model(self):
        return DielectricModelVariant1.create_model(
                                                BrainSubstance.GRAY_MATTER)

    def test_permitivity_0Hz(self, model):
        result = model.permitivity(frequency=0)
        np.testing.assert_allclose(result, 45200449, atol=1)

    def test_permitivity_1Hz(self, model):
        result = model.permitivity(frequency=1*2*np.pi)
        np.testing.assert_allclose(result, 45150280, atol=1)

    def test_permitivity_1kHz(self, model):
        result = model.permitivity(frequency=1000*2*np.pi)
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
