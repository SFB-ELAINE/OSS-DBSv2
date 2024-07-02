import pytest

from ossdbs.electrodes import NeuroPaceDL344_10

from .test_electrodes import TestElectrode


class TestNeuroPaceDL344_10(TestElectrode):
    @pytest.fixture
    def electrode(self):
        return NeuroPaceDL344_10()

    @pytest.fixture
    def electrode_name(self):
        return "NeuroPaceDL344_10"

    def test_rename_boundaries(self, electrode, electrode_name):
        """Test whether set_contact_names() works."""
        self.check_rename_boundaries(electrode, electrode_name)

    def test_contacts(self, electrode, electrode_name):
        """Test the number and names of contacts."""
        self.check_contacts(electrode, electrode_name)

    def test_electrode_volume(self, electrode, electrode_name):
        """Test volume of the entire electrode."""
        self.check_electrode_volume(electrode, electrode_name)

    def test_contacts_volume(self, electrode, electrode_name):
        """Test volume of all the contacts."""
        self.check_contacts_volume(electrode, electrode_name)
