import netgen
import ngsolve
import numpy as np
import pytest

from ossdbs.electrodes import AbbottStJudeDirected6173


class TestAbbottStJudeDirected6173:
    @pytest.fixture
    def AbbottStJudeDirected6173_electrode(self):
        return AbbottStJudeDirected6173()

    def test_rename_boundaries(self, AbbottStJudeDirected6173_electrode):
        """Test whether set_contact_names() works."""
        electrode = AbbottStJudeDirected6173_electrode
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
        desired = {
            "RenamedBody",
            "RenamedContact_1",
            "Contact_2",
            "Contact_3",
            "Contact_4",
            "Contact_5",
            "Contact_6",
            "Contact_7",
            "Contact_8",
        }
        assert desired == set(mesh.GetBoundaries())

    def test_contacts(self, AbbottStJudeDirected6173_electrode):
        """Test the number and names of contacts."""
        electrode = AbbottStJudeDirected6173_electrode
        geometry = electrode.geometry
        netgen_geometry = netgen.occ.OCCGeometry(geometry)
        with ngsolve.TaskManager():
            mesh = ngsolve.Mesh(netgen_geometry.GenerateMesh())
        desired = {
            "Body",
            "Contact_1",
            "Contact_2",
            "Contact_3",
            "Contact_4",
            "Contact_5",
            "Contact_6",
            "Contact_7",
            "Contact_8",
        }
        assert desired == set(mesh.GetBoundaries())

    def test_electrode_volume(self, AbbottStJudeDirected6173_electrode):
        """Test volume of the entire electrode."""
        electrode = AbbottStJudeDirected6173_electrode

        total_length = electrode._parameters.total_length
        tip_length = electrode._parameters.tip_length
        radius = electrode._parameters.lead_diameter * 0.5
        height = total_length - tip_length

        desired = (np.pi * radius**2 * height) + (4 / 3 * np.pi * radius**3 * 0.5)
        actual = electrode.geometry.mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)

    def test_contacts_volume(self, AbbottStJudeDirected6173_electrode):
        """Test volume of all the contacts."""
        electrode = AbbottStJudeDirected6173_electrode

        contact_length = electrode._parameters.contact_length
        radius = electrode._parameters.lead_diameter * 0.5

        desired = (np.pi * radius**2 * contact_length * 2) + (
            np.pi * radius**2 * contact_length * 90 / 360
        ) * 6
        actual = electrode._contacts().mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)
