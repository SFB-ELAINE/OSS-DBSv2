import numpy as np
import pytest

from ossdbs.electrodes import BostonScientificVercise


class TestBostonScientificVercise:
    @pytest.fixture
    def BostonScientificVercise_electrode(self):
        return BostonScientificVercise()

    def test_rename_boundaries(self, BostonScientificVercise_electrode):
        """Test whether set_contact_names() works."""
        electrode = BostonScientificVercise_electrode
        electrode.set_contact_names(
            {
                "Body": "RenamedBody",
                "Contact_1": "RenamedContact_1",
                "NonExistingPart": "NonExistingPart",
            }
        )
        geometry = electrode.geometry
        faces = [face.name for face in geometry.faces]
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
        assert desired == set(faces)

    def test_contacts(self, BostonScientificVercise_electrode):
        """Test the number and names of contacts."""
        electrode = BostonScientificVercise_electrode
        geometry = electrode.geometry
        faces = [face.name for face in geometry.faces]
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
        assert desired == set(faces)

    def test_electrode_volume(self, BostonScientificVercise_electrode):
        """Test volume of the entire electrode."""
        electrode = BostonScientificVercise_electrode

        total_length = electrode._parameters.total_length
        tip_length = electrode._parameters.tip_length
        radius = electrode._parameters.lead_diameter * 0.5
        height = total_length - tip_length

        desired = (np.pi * radius**2 * height) + (4 / 3 * np.pi * radius**3 * 0.5)
        actual = electrode.geometry.mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)

    def test_contacts_volume(self, BostonScientificVercise_electrode):
        """Test volume of all the contacts."""
        electrode = BostonScientificVercise_electrode

        contact_length = electrode._parameters.contact_length
        radius = electrode._parameters.lead_diameter * 0.5
        n_contacts = electrode._n_contacts

        desired = (contact_length * radius**2 * np.pi) * n_contacts
        actual = electrode._contacts().mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)
