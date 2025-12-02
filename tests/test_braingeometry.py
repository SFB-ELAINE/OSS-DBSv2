from typing import ClassVar

import numpy as np
import pytest
from netgen.occ import Box, Pnt

import ossdbs


class TestBrainGeometry:
    """Class for testing brain geometry."""

    region_parameters: ClassVar[dict[str, dict[str, float]]] = {
        "Center": {"x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 0.0},
        "Dimension": {"x[mm]": 10.0, "y[mm]": 20.0, "z[mm]": 30.0},
    }

    TESTDATA: ClassVar[list[tuple[dict[str, dict[str, float]], str]]] = [
        # region_parameters, shape
        (region_parameters, "Box"),
        (region_parameters, "Sphere"),
        (region_parameters, "Ellipsoid"),
    ]

    @pytest.mark.parametrize("region_parameters, shape", TESTDATA)
    def test_geometry(self, region_parameters, shape) -> None:
        brain_region = ossdbs.create_bounding_box(region_parameters)
        geometry = ossdbs.BrainGeometry(shape, brain_region)

        actual = geometry.geometry.mass
        tolerance = 1e-3
        if geometry._shape == "Sphere":
            x, y, z = np.subtract(geometry._bbox.end, geometry._bbox.start) / 2
            radius = np.min([x, y, z])
            desired = 4 / 3 * np.pi * radius**3

            np.testing.assert_allclose(actual, desired, atol=tolerance)
        elif geometry._shape == "Ellipsoid":
            x, y, z = np.subtract(geometry._bbox.end, geometry._bbox.start) / 2
            desired = 4 / 3 * np.pi * x * y * z

            np.testing.assert_allclose(actual, desired, atol=tolerance)
        elif geometry._shape == "Box":
            x = region_parameters["Dimension"]["x[mm]"]
            y = region_parameters["Dimension"]["y[mm]"]
            z = region_parameters["Dimension"]["z[mm]"]
            desired = x * y * z

            np.testing.assert_allclose(actual, desired, atol=tolerance)

    @pytest.fixture
    def brain_geometry(self):
        """Sample brain geometry."""
        region_parameters = {
            "Center": {"x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 0.0},
            "Dimension": {"x[mm]": 20.0, "y[mm]": 20.0, "z[mm]": 20.0},
        }

        brain_region = ossdbs.create_bounding_box(region_parameters)
        return ossdbs.BrainGeometry("Sphere", brain_region)

    def test_get_surface_names(self, brain_geometry):
        """Test get_surface_names()."""
        assert set(brain_geometry.get_surface_names()) == {"BrainSurface"}

    def test_set_geometry(self, brain_geometry):
        """Test set_geometry()."""
        external_geometry = Box(Pnt(-10, -10, -10), (10, 10, 10))
        brain_geometry.set_geometry(external_geometry)

        assert {face.name for face in brain_geometry._geometry.faces} == {
            "BrainSurface"
        }
