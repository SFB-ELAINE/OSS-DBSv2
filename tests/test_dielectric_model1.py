from ossdbs.dielectric_model import ColeCole4Model
import numpy as np
import pytest


class TestWhiteMatterModel:

    TESTDATA_PERMITTIVITY = [(0, 0.00031, 1e-5),
                             (1, 0.000309 - 0.003199j, 1e-5),
                             (1000, 6.18117e-07 - 9.959125e-06j, 1e-8)]

    TESTDATA_CONDUCTIVITY = [(0, 0.02, 1e-5),
                             (1, 0.020103 + 0.001941j, 1e-5),
                             (1000, 0.062575 + 0.003884j, 1e-5)]

    @pytest.mark.parametrize('frequency, permittivity, tolerance',
                             TESTDATA_PERMITTIVITY)
    def test_permittivity(self, frequency, permittivity, tolerance):
        omega = 2 * np.pi * frequency
        actual = ColeCole4Model().complex_permittivity("White matter", omega)
        np.testing.assert_allclose(actual, permittivity, atol=tolerance)

    @pytest.mark.parametrize('frequency, conductivity, tolerance',
                             TESTDATA_CONDUCTIVITY)
    def test_conductivity(self, frequency, conductivity, tolerance):
        omega = 2 * np.pi * frequency
        actual = ColeCole4Model().complex_conductivity("White matter", omega)
        np.testing.assert_allclose(actual, conductivity, atol=tolerance)


class TestGrayMatterModel:

    TESTDATA_PERMITTIVITY = [(0, 0.0004, 1e-5),
                             (1, 0.0004 - 0.003196j, 1e-5),
                             (1000, 1.452645e-06 - 1.572557e-05j, 1e-7)]

    TESTDATA_CONDUCTIVITY = [(0, 0.02, 1e-5),
                             (1, 0.020083 + 0.002512j, 1e-5),
                             (1000, 0.098807 + 0.009127j, 1e-5)]

    @pytest.mark.parametrize('frequency, permittivity, tolerance',
                             TESTDATA_PERMITTIVITY)
    def test_permittivity(self, frequency, permittivity, tolerance):
        omega = 2 * np.pi * frequency
        actual = ColeCole4Model().complex_permittivity("Gray matter", omega)
        np.testing.assert_allclose(actual, permittivity, atol=tolerance)

    @pytest.mark.parametrize('frequency, conductivity, tolerance',
                             TESTDATA_CONDUCTIVITY)
    def test_conductivity(self, frequency, conductivity, tolerance):
        omega = 2 * np.pi * frequency
        actual = ColeCole4Model().complex_conductivity("Gray matter", omega)
        np.testing.assert_allclose(actual, conductivity, atol=tolerance)


class TestCerebroSpinalFluidModel:

    TESTDATA_PERMITTIVITY = [(0, 9.651065e-10, 1e-11),
                             (1, 9.651065e-10 - 0.31831j, 1e-5),
                             (1000, 9.651064e-10 - 0.000318j, 1e-5)]

    TESTDATA_CONDUCTIVITY = [(0, 2.0, 1e-5),
                             (1, 2.0 + 6.063943e-09j, 1e-10),
                             (1000, 2.0 + 6.063943e-06j, 1e-7)]

    @pytest.mark.parametrize('frequency, permittivity, tolerance',
                             TESTDATA_PERMITTIVITY)
    def test_permittivity(self, frequency, permittivity, tolerance):
        omega = 2 * np.pi * frequency
        actual = ColeCole4Model().complex_permittivity("CSF", omega)
        np.testing.assert_allclose(actual, permittivity, atol=tolerance)

    @pytest.mark.parametrize('frequency, conductivity, tolerance',
                             TESTDATA_CONDUCTIVITY)
    def test_conductivity(self, frequency, conductivity, tolerance):
        omega = 2 * np.pi * frequency
        actual = ColeCole4Model().complex_conductivity("CSF", omega)
        np.testing.assert_allclose(actual, conductivity, atol=tolerance)
