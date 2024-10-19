import pytest

from ossdbs.electrodes import ELECTRODE_MODELS, BostonScientificCartesiaX

from .test_electrodes import TestElectrode


class TestBostonScientificCartesiaX(TestElectrode):
    @pytest.fixture
    def electrode(self):
        return BostonScientificCartesiaX()

    @pytest.fixture
    def electrode_name(self):
        return "BostonScientificCartesiaX"

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

    # TODO check electrode orientation and test it
    '''
    def test_spin_direction_0(self, electrode_name):
        """Test first spin direction."""
        idx = 0
        ref_geo, model_geo = get_reference_and_model_geo(electrode_name, idx)
        reference_centers = {}
        for face in ref_geo.faces:
            if "E1C" in face.name:
                if face.name not in ["E1C1", "E1C8"]:
                    print(type(face.center))
                    reference_centers[face.name] = np.array([face.center.x, face.center.y, face.center.z])
        for face in model_geo.geometry.shape.faces:
            if "E1C" in face.name:
                if face.name not in ["E1C1", "E1C8"]:
                    print(reference_centers[face.name].dtype)
                    assert np.all(np.isclose(reference_centers[face.name],
                                             np.array([face.center.x, face.center.y, face.center.z]), rtol=1e-4))

    def test_spin_direction_1(self, electrode_name):
        """Test first spin direction."""
        idx = 1
        ref_geo, model_geo = get_reference_and_model_geo(electrode_name, idx)
        reference_centers = {}
        for face in ref_geo.faces:
            if "E1C" in face.name:
                if face.name not in ["E1C1", "E1C8"]:
                    print(type(face.center))
                    reference_centers[face.name] = np.array([face.center.x, face.center.y, face.center.z])
        for face in model_geo.geometry.shape.faces:
            if "E1C" in face.name:
                if face.name not in ["E1C1", "E1C8"]:
                    print(reference_centers[face.name].dtype)
                    assert np.all(np.isclose(reference_centers[face.name],
                                             np.array([face.center.x, face.center.y, face.center.z]), rtol=1e-4))

    def test_spin_direction_2(self, electrode_name):
        """Test first spin direction."""
        idx = 2
        ref_geo, model_geo = get_reference_and_model_geo(electrode_name, idx)
        reference_centers = {}
        for face in ref_geo.faces:
            if "E1C" in face.name:
                if face.name not in ["E1C1", "E1C8"]:
                    print(type(face.center))
                    reference_centers[face.name] = np.array([face.center.x, face.center.y, face.center.z])
        for face in model_geo.geometry.shape.faces:
            if "E1C" in face.name:
                if face.name not in ["E1C1", "E1C8"]:
                    print(reference_centers[face.name].dtype)
                    assert np.all(np.isclose(reference_centers[face.name],
                                             np.array([face.center.x, face.center.y, face.center.z]), rtol=1e-4))

    def test_spin_direction_3(self, electrode_name):
        """Test first spin direction."""
        idx = 3
        ref_geo, model_geo = get_reference_and_model_geo(electrode_name, idx)
        reference_centers = {}
        for face in ref_geo.faces:
            if "E1C" in face.name:
                if face.name not in ["E1C1", "E1C8"]:
                    print(type(face.center))
                    reference_centers[face.name] = np.array([face.center.x, face.center.y, face.center.z])
        for face in model_geo.geometry.shape.faces:
            if "E1C" in face.name:
                if face.name not in ["E1C1", "E1C8"]:
                    print(reference_centers[face.name].dtype)
                    assert np.all(np.isclose(reference_centers[face.name],
                                             np.array([face.center.x, face.center.y, face.center.z]), rtol=1e-4))
    '''
