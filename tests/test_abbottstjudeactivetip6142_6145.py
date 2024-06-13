import netgen
import ngsolve
import numpy as np
import pytest

from ossdbs.electrodes import AbbottStJudeActiveTip6142_6145


class TestAbbottStJudeActiveTip6142_6145:
    @pytest.fixture
    def AbbottStJudeActiveTip6142_6145_electrode(self):
        return AbbottStJudeActiveTip6142_6145()

    def test_rename_boundaries(self, AbbottStJudeActiveTip6142_6145_electrode):
        """Test whether set_contact_names() works."""
        electrode = AbbottStJudeActiveTip6142_6145_electrode
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
        }
        assert desired == set(mesh.GetBoundaries())

    def test_contacts(self, AbbottStJudeActiveTip6142_6145_electrode):
        """Test the number and names of contacts."""
        electrode = AbbottStJudeActiveTip6142_6145_electrode
        geometry = electrode.geometry
        netgen_geometry = netgen.occ.OCCGeometry(geometry)
        with ngsolve.TaskManager():
            mesh = ngsolve.Mesh(netgen_geometry.GenerateMesh())
        desired = {"Body", "Contact_1", "Contact_2", "Contact_3", "Contact_4"}
        assert desired == set(mesh.GetBoundaries())

    def test_electrode_volume(self, AbbottStJudeActiveTip6142_6145_electrode):
        """Test volume of the entire electrode."""
        electrode = AbbottStJudeActiveTip6142_6145_electrode

        total_length = electrode._parameters.total_length
        tip_length = electrode._parameters.tip_length
        radius = electrode._parameters.lead_diameter * 0.5
        height = total_length - tip_length

        desired = (np.pi * radius**2 * height) + (4 / 3 * np.pi * radius**3 * 0.5)
        actual = electrode.geometry.mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)

    def test_contacts_volume(self, AbbottStJudeActiveTip6142_6145_electrode):
        """Test volume of all the contacts."""
        electrode = AbbottStJudeActiveTip6142_6145_electrode

        contact_length = electrode._parameters.contact_length
        radius = electrode._parameters.lead_diameter * 0.5
        n_contacts = electrode._n_contacts

        C1_height = electrode._parameters.tip_length - radius
        C1_volume = (4 / 3 * radius**3 * np.pi * 0.5) + (
            C1_height * radius**2 * np.pi
        )

        desired = (contact_length * radius**2 * np.pi) * (n_contacts - 1) + C1_volume
        actual = electrode._contacts().mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)
