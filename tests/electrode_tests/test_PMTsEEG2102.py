import pytest

from ossdbs.electrodes import (
    ELECTRODE_MODELS,
    PMTsEEG2102_08,
    PMTsEEG2102_10,
    PMTsEEG2102_12,
    PMTsEEG2102_14,
    PMTsEEG2102_16,
)

from .test_electrodes import TestElectrode


class TestPMTsEEG2102_08(TestElectrode):
    @pytest.fixture
    def electrode(self):
        return PMTsEEG2102_08()

    @pytest.fixture
    def electrode_name(self):
        return "PMTsEEG2102_08"

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

    def test_custom_exists(self, electrode_name):
        customname = electrode_name + "Custom"
        assert customname in ELECTRODE_MODELS.keys()


class TestPMTsEEG2102_10(TestPMTsEEG2102_08):
    @pytest.fixture
    def electrode(self):
        return PMTsEEG2102_10()

    @pytest.fixture
    def electrode_name(self):
        return "PMTsEEG2102_10"


class TestPMTsEEG2102_12(TestPMTsEEG2102_08):
    @pytest.fixture
    def electrode(self):
        return PMTsEEG2102_12()

    @pytest.fixture
    def electrode_name(self):
        return "PMTsEEG2102_12"


class TestPMTsEEG2102_14(TestPMTsEEG2102_08):
    @pytest.fixture
    def electrode(self):
        return PMTsEEG2102_14()

    @pytest.fixture
    def electrode_name(self):
        return "PMTsEEG2102_14"


class TestPMTsEEG2102_16(TestPMTsEEG2102_08):
    @pytest.fixture
    def electrode(self):
        return PMTsEEG2102_16()

    @pytest.fixture
    def electrode_name(self):
        return "PMTsEEG2102_16"
