from ossdbs.dielectric_model import ColeColeFourModelFactory
from ossdbs.brainsubstance import Material
import numpy as np
import pytest


class TestModel1Permitivity:

    TESTDATA = [
        (Material.WHITE_MATTER, 0, 0.00031, 1e-5),
        (Material.WHITE_MATTER, 1, 0.000309 - 0.003199j, 1e-5),
        (Material.WHITE_MATTER, 1000, 6.18117e-07 - 9.959125e-06j, 1e-8),
        (Material.GRAY_MATTER, 0, 0.0004, 1e-5),
        (Material.GRAY_MATTER, 1, 0.0004 - 0.003196j, 1e-5),
        (Material.GRAY_MATTER, 1000, 1.452645e-06 - 1.572557e-05j, 1e-7),
        (Material.CSF, 0, 9.651065e-10, 1e-11),
        (Material.CSF, 1, 9.651065e-10 - 0.31831j, 1e-5),
        (Material.CSF, 1000, 9.651064e-10 - 0.000318j, 1e-5)
    ]

    @pytest.mark.parametrize('material, frequency, permitivity, tolerance',
                             TESTDATA)
    def test_permitivity(self, material, frequency, permitivity, tolerance):
        model = ColeColeFourModelFactory.create(material)
        actual = model.permitivity(omega=2*np.pi*frequency)
        np.testing.assert_allclose(actual, permitivity, atol=tolerance)


class TestModel1Conductivity:

    TESTDATA = [
        (Material.WHITE_MATTER, 0, 0.02, 1e-5),
        (Material.WHITE_MATTER, 1, 0.020103 - 0.001941j, 1e-5),
        (Material.WHITE_MATTER, 1000, 0.062575 - 0.003884j, 1e-5),
        (Material.GRAY_MATTER, 0, 0.02, 1e-5),
        (Material.GRAY_MATTER, 1, 0.020083 - 0.002512j, 1e-5),
        (Material.GRAY_MATTER, 1000, 0.098807-0.009127j, 1e-5),
        (Material.CSF, 0, 2.0, 1e-5),
        (Material.CSF, 1, 2.0 - 6.063943e-09j, 1e-10),
        (Material.CSF, 1000, 2.0 - 6.063943e-06j, 1e-7)
    ]

    @pytest.mark.parametrize('material, frequency, conductivity, tolerance',
                             TESTDATA)
    def test_permitivity(self, material, frequency, conductivity, tolerance):
        model = ColeColeFourModelFactory.create(material)
        actual = model.conductivity(omega=2*np.pi*frequency)
        np.testing.assert_allclose(actual, conductivity, atol=tolerance)
