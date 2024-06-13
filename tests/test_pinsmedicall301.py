import netgen
import ngsolve
import numpy as np
import pytest

from ossdbs.electrodes import PINSMedicalL301


class TestPINSMedicalL301:
    @pytest.fixture
    def PINSMedicalL301_electrode(self):
        return PINSMedicalL301()

    def test_rename_boundaries(self, PINSMedicalL301_electrode):
        """Test whether set_contact_names() works."""
        electrode = PINSMedicalL301_electrode
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

    def test_contacts(self, PINSMedicalL301_electrode):
        """Test the number and names of contacts."""
        electrode = PINSMedicalL301_electrode
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
        }
        assert desired == set(mesh.GetBoundaries())

    def test_electrode_volume(self, PINSMedicalL301_electrode):
        """Test volume of the entire electrode."""
        electrode = PINSMedicalL301_electrode

        total_length = electrode._parameters.total_length
        tip_length = electrode._parameters.tip_length
        radius = electrode._parameters.lead_diameter * 0.5
        height = total_length - tip_length

        desired = (np.pi * radius**2 * height) + (4 / 3 * np.pi * radius**3 * 0.5)
        actual = electrode.geometry.mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)

    def test_contacts_volume(self, PINSMedicalL301_electrode):
        """Test volume of all the contacts."""
        electrode = PINSMedicalL301_electrode
        contact_length = electrode._parameters.contact_length
        radius = electrode._parameters.lead_diameter * 0.5
        n_contacts = electrode._n_contacts

        desired = (contact_length * radius**2 * np.pi) * n_contacts
        actual = electrode._contacts().mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)
