import netgen
import ngsolve
import numpy as np
import pytest

from ossdbs.electrodes import MicroProbesSNEX100


class TestMicroProbesSNEX_100:
    @pytest.fixture
    def MicroProbesSNEX100_electrode(self):
        return MicroProbesSNEX100()

    def test_rename_boundaries(self, MicroProbesSNEX100_electrode):
        """Test whether set_contact_names() works."""
        electrode = MicroProbesSNEX100_electrode
        electrode.set_contact_names(
            {
                "Body": "RenamedBody",
                "Contact_1": "RenamedContact_1",
                "NonExistingPart": "NonExistingPart",
            }
        )
        geometry = electrode.geometry
        netgen_geometry = netgen.occ.OCCGeometry(geometry)
        with ngsolve.TaskManager():
            mesh = ngsolve.Mesh(netgen_geometry.GenerateMesh())
        desired = {"RenamedBody", "RenamedContact_1", "Contact_2"}
        assert desired == set(mesh.GetBoundaries())

    def test_contacts(self, MicroProbesSNEX100_electrode):
        """Test the number and names of contacts."""
        electrode = MicroProbesSNEX100_electrode
        geometry = electrode.geometry
        netgen_geometry = netgen.occ.OCCGeometry(geometry)
        with ngsolve.TaskManager():
            mesh = ngsolve.Mesh(netgen_geometry.GenerateMesh())
        desired = {
            "Body",
            "Contact_1",
            "Contact_2",
        }
        assert desired == set(mesh.GetBoundaries())

    def test_electrode_volume(self, MicroProbesSNEX100_electrode):
        """Test volume of the entire electrode."""
        electrode = MicroProbesSNEX100_electrode

        radius_1 = electrode._parameters.core_tubing_diameter * 0.5
        height_1 = electrode._parameters.core_tubing_length

        distance_2 = (
            electrode._parameters.core_electrode_length
            + electrode._parameters.core_tubing_length
            + electrode._parameters.outer_electrode_length
        )

        radius_2 = electrode._parameters.outer_tubing_diameter * 0.5
        height_2 = electrode._parameters.total_length - distance_2

        body1 = np.pi * radius_1**2 * height_1
        body2 = np.pi * radius_2**2 * height_2

        outer_electrode_radius = electrode._parameters.outer_electrode_diameter * 0.5
        outer_electrode_length = electrode._parameters.outer_electrode_length

        core_electrode_radius = electrode._parameters.core_electrode_diameter * 0.5
        core_electrode_length = electrode._parameters.core_electrode_length

        contact_1 = np.pi * outer_electrode_radius**2 * outer_electrode_length
        contact_2 = (4 / 3 * np.pi * core_electrode_radius**3 * 0.5) + (
            np.pi
            * core_electrode_radius**2
            * (core_electrode_length - core_electrode_radius)
        )

        desired = body1 + body2 + contact_1 + contact_2
        actual = electrode.geometry.mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)

    def test_contacts_volume(self, MicroProbesSNEX100_electrode):
        """Test volume of all the contacts."""
        electrode = MicroProbesSNEX100_electrode

        outer_electrode_radius = electrode._parameters.outer_electrode_diameter * 0.5
        outer_electrode_length = electrode._parameters.outer_electrode_length

        core_electrode_radius = electrode._parameters.core_electrode_diameter * 0.5
        core_electrode_length = electrode._parameters.core_electrode_length

        contact_1 = np.pi * outer_electrode_radius**2 * outer_electrode_length
        contact_2 = (4 / 3 * np.pi * core_electrode_radius**3 * 0.5) + (
            np.pi
            * core_electrode_radius**2
            * (core_electrode_length - core_electrode_radius)
        )

        desired = contact_1 + contact_2
        actual = electrode._contacts().mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)
