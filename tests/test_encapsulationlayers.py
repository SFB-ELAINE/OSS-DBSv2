import numpy as np
import pytest

from ossdbs.model_geometry.encapsulation_layers import EncapsulationLayer


class TestEncapsulationLayer:
    @pytest.fixture
    def encapsulationLayer(self):
        return EncapsulationLayer(
            name="EncapsulationLayer_0", dielectric_model="ColeCole4", material="Blood"
        )

    def test_dielectric_peoperties(self, encapsulationLayer):
        desired = {
            4.0,
            0.7,
            tuple(np.array([0.1, 0.1, 0.0, 0.0])),
            tuple(np.array([56.0, 5200.0, 0.0, 0.0])),
            tuple(np.array([8.38e-12, 132.63e-9, 0, 0])),
        }

        actual = set()
        properties = encapsulationLayer.dielectric_properties._parameters
        actual.add(properties.eps_inf)
        actual.add(properties.sigma)
        actual.add(tuple(properties.alpha))
        actual.add(tuple(properties.eps_delta))
        actual.add(tuple(properties.tau))

        assert actual == desired
