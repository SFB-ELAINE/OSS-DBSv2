import numpy as np
import pytest

from ossdbs.electrodes import ELECTRODE_MODELS, NeuroNexusA1x16_5mm_50_177

from .test_electrodes import TestElectrode


class TestNeuroNexusA1x16_5mm_50_177(TestElectrode):
    @pytest.fixture
    def electrode(self):
        return NeuroNexusA1x16_5mm_50_177()

    @pytest.fixture
    def electrode_name(self):
        return "NeuroNexusA1x16_5mm_50_177"

    def test_rename_boundaries(self, electrode, electrode_name):
        """Test whether set_contact_names() works."""
        self.check_rename_boundaries(electrode, electrode_name)

    def test_contacts(self, electrode, electrode_name):
        """Test the number and names of contacts."""
        self.check_contacts(electrode, electrode_name)

    def test_electrode_volume(self, electrode, electrode_name):
        """Test volume of the entire electrode."""
        self.check_electrode_volume(electrode, electrode_name)

    # Do not test this here because contacts have no volume
    '''
    def test_contacts_volume(self, electrode, electrode_name):
        """Test volume of all the contacts."""
        self.check_contacts_volume(electrode, electrode_name)
    '''

    def test_custom_exists(self, electrode_name):
        customname = electrode_name + "Custom"
        assert customname in ELECTRODE_MODELS.keys()

    def test_contact_area(self, electrode):
        """Check the contact area."""
        geometry = electrode.geometry
        contact_diameter = electrode._parameters.contact_diameter

        desired = np.pi * (0.5 * contact_diameter) ** 2
        actual = []
        for face in geometry.faces:
            if face.name is not None and "Contact" in face.name:
                actual.append(face.mass)
        assert np.all(np.isclose(np.array(actual), desired))
