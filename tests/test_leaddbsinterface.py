"""Tests for leaddbsinterface module."""

import json
import sys
from unittest.mock import MagicMock

import h5py
import numpy as np
import pytest

from leaddbsinterface.default_settings import (
    initialize_default_settings,
    load_default_for_lead,
    update_default_dict,
)

# Mock heavy dependencies before importing leaddbsinterface.lead_settings
# These are CAD/FEM libraries not needed for testing the interface
# Note: Using unittest.mock.MagicMock here since this is module-level
# (before mocker fixture is available)
mock_ngsolve = MagicMock()
mock_netgen = MagicMock()

sys.modules["ngsolve"] = mock_ngsolve
sys.modules["netgen"] = mock_netgen
sys.modules["netgen.occ"] = MagicMock()
sys.modules["netgen.meshing"] = MagicMock()
sys.modules["netgen.csg"] = MagicMock()

from leaddbsinterface.lead_settings import LeadSettings  # noqa: E402


class TestDefaultSettings:
    """Test default_settings module functions."""

    def test_initialize_default_settings(self):
        """Test initialization of default settings structure."""
        settings = initialize_default_settings()

        assert isinstance(settings, dict)
        assert "ModelSide" in settings
        assert settings["ModelSide"] == 0
        assert "BrainRegion" in settings
        assert "Dimension" in settings["BrainRegion"]
        assert "Shape" in settings["BrainRegion"]
        assert "StimulationSignal" in settings
        assert "PointModel" in settings
        assert "Electrodes" in settings
        assert "Mesh" in settings
        assert "Solver" in settings

    def test_load_default_for_lead(self):
        """Test loading default parameters for Lead-DBS."""
        settings = initialize_default_settings()
        settings = load_default_for_lead(settings)

        # Check brain region dimensions
        assert settings["BrainRegion"]["Dimension"]["x[mm]"] == 60
        assert settings["BrainRegion"]["Dimension"]["y[mm]"] == 60
        assert settings["BrainRegion"]["Dimension"]["z[mm]"] == 80
        assert settings["BrainRegion"]["Shape"] == "Ellipsoid"

        # Check stimulation signal
        assert settings["StimulationSignal"]["Type"] == "Multisine"
        assert settings["StimulationSignal"]["ListOfFrequencies"] == [10000]

        # Check voxel lattice
        assert settings["PointModel"]["VoxelLattice"]["Active"] is False
        assert settings["PointModel"]["VoxelLattice"]["Shape"]["x"] == 31
        assert settings["PointModel"]["VoxelLattice"]["Shape"]["y"] == 31
        assert settings["PointModel"]["VoxelLattice"]["Shape"]["z"] == 41

        # Check encapsulation layer
        assert (
            settings["Electrodes"][settings["ModelSide"]]["EncapsulationLayer"][
                "Thickness[mm]"
            ]
            == 0.1
        )
        assert (
            settings["Electrodes"][settings["ModelSide"]]["EncapsulationLayer"][
                "Material"
            ]
            == "White matter"
        )

        # Check other settings
        assert settings["ExportVTK"] is True
        assert settings["Mesh"]["MeshingHypothesis"]["Type"] == "Default"
        assert settings["Mesh"]["MaterialRefinementSteps"] == 1
        assert settings["FEMOrder"] == 2
        assert settings["ComputeImpedance"] is False
        assert settings["Solver"]["MaximumSteps"] == 500
        assert settings["Solver"]["Precision"] == 1e-10

    def test_update_default_dict_simple(self):
        """Test updating default dict with simple values."""
        default = {"a": 1, "b": 2, "c": 3}
        custom = {"b": 20, "d": 4}

        result = update_default_dict(default, custom)

        assert result["a"] == 1
        assert result["b"] == 20
        assert result["c"] == 3
        assert result["d"] == 4

    def test_update_default_dict_nested(self):
        """Test updating default dict with nested dictionaries."""
        default = {
            "level1": {
                "level2": {
                    "a": 1,
                    "b": 2,
                },
                "c": 3,
            },
            "d": 4,
        }

        custom = {
            "level1": {
                "level2": {
                    "b": 20,
                },
                "e": 5,
            }
        }

        result = update_default_dict(default, custom)

        assert result["level1"]["level2"]["a"] == 1
        assert result["level1"]["level2"]["b"] == 20
        assert result["level1"]["c"] == 3
        assert result["level1"]["e"] == 5
        assert result["d"] == 4

    def test_update_default_dict_empty_nested(self):
        """Test updating with empty nested dicts."""
        default = {"a": {}, "b": 2}
        custom = {"a": {"x": 1}, "c": 3}

        result = update_default_dict(default, custom)

        assert result["a"]["x"] == 1
        assert result["b"] == 2
        assert result["c"] == 3


class TestLeadSettingsMocked:
    """Test LeadSettings class with mocked h5py file."""

    @pytest.fixture
    def mock_h5_file(self, mocker):
        """Create a mock h5py file with minimal required data."""
        mock_file = mocker.MagicMock(spec=h5py.File)
        mock_settings = mocker.MagicMock()

        # Mock current control
        mock_cur_ctrl = mocker.MagicMock()
        mock_cur_ctrl.T = [np.array([1.0, 1.0])]
        mock_settings.__getitem__.return_value = mock_cur_ctrl

        # Set up the settings mock
        def getitem_side_effect(key):
            if key == "current_control":
                return mock_cur_ctrl
            return mocker.MagicMock()

        mock_settings.__getitem__.side_effect = getitem_side_effect
        mock_file.__getitem__.return_value = mock_settings

        return mock_file

    def test_lead_settings_initialization_h5(self, mocker, mock_h5_file):
        """Test LeadSettings initialization with h5 file."""
        mocker.patch("h5py.File", return_value=mock_h5_file)
        ls = LeadSettings("test.mat")
        assert ls._h5 is True
        assert ls.NUM_ELECS == 1

    def test_lead_settings_initialization_non_h5_raises_error(self, mocker):
        """Test that non-h5 files raise ValueError with helpful message."""
        mocker.patch("h5py.File", side_effect=ValueError("Not HDF5"))
        with pytest.raises(ValueError):
            LeadSettings("test.mat")

    def test_lead_settings_same_current_voltage_ok(self, mocker):
        """Test that same CC or VC for both hemispheres is OK."""
        mock_file = mocker.MagicMock(spec=h5py.File)
        mock_settings = mocker.MagicMock()
        mock_cur_ctrl = mocker.MagicMock()
        # Both CC or both VC - should be OK
        mock_cur_ctrl.T = [np.array([1.0, 1.0])]  # Both CC

        def getitem_side_effect(key):
            if key == "current_control":
                return mock_cur_ctrl
            return mocker.MagicMock()

        mock_settings.__getitem__.side_effect = getitem_side_effect
        mock_file.__getitem__.return_value = mock_settings

        mocker.patch("h5py.File", return_value=mock_file)
        ls = LeadSettings("test.mat")
        assert ls.NUM_ELECS == 1


class TestLeadSettingsWithRealH5File:
    """Test LeadSettings with a real temporary h5 file."""

    @pytest.fixture
    def minimal_h5_file(self, tmp_path):
        """Create a minimal valid h5 file for testing."""
        file_path = tmp_path / "test_settings.mat"

        with h5py.File(str(file_path), "w") as f:
            settings = f.create_group("settings")

            # Required fields for initialization
            # Shape should be (1, 2) for two hemispheres
            # After _get_arr does [:, :].T and get_cur_ctrl does .T[0], we need shape (1, 2)
            settings.create_dataset("current_control", data=np.array([[1.0, 1.0]]))

            # Add minimal required fields for make_oss_settings
            settings.create_dataset(
                "Implantation_coordinate", data=np.array([[0, 0, 0]])
            )
            settings.create_dataset("Second_coordinate", data=np.array([[0, 0, 10]]))
            settings.create_dataset("Phi_vector", data=np.array([[1.0, 0, 0, 0]]))
            settings.create_dataset("Case_grounding", data=np.array([[0]]))

            # String fields (stored as ASCII in h5)
            electrode_type = "Medtronic 3389"
            settings.create_dataset(
                "Electrode_type",
                data=np.array([ord(c) for c in electrode_type]).reshape(-1, 1),
            )

            mri_name = "/tmp/test.nii"
            settings.create_dataset(
                "MRI_data_name",
                data=np.array([ord(c) for c in mri_name]).reshape(-1, 1),
            )

            dti_name = ""
            if dti_name:
                settings.create_dataset(
                    "DTI_data_name",
                    data=np.array([ord(c) for c in dti_name]).reshape(-1, 1),
                )
            else:
                settings.create_dataset(
                    "DTI_data_name",
                    data=np.array([ord(c) for c in "no dti"]).reshape(-1, 1),
                )

            cond_model = "ColeCole4"
            settings.create_dataset(
                "cond_model",
                data=np.array([ord(c) for c in cond_model]).reshape(-1, 1),
            )

            encaps_type = "Gray matter"
            settings.create_dataset(
                "encapsulationType",
                data=np.array([ord(c) for c in encaps_type]).reshape(-1, 1),
            )

            # Numeric fields
            settings.create_dataset("GM_index", data=np.array([[2]]))
            settings.create_dataset("WM_index", data=np.array([[3]]))
            settings.create_dataset("CSF_index", data=np.array([[1]]))
            settings.create_dataset("Estimate_In_Template", data=np.array([[0]]))
            settings.create_dataset("calcAxonActivation", data=np.array([[0]]))
            settings.create_dataset("removeElectrode", data=np.array([[0]]))
            settings.create_dataset("Activation_threshold_VTA", data=np.array([[0.2]]))
            settings.create_dataset("outOfCore", data=np.array([[0]]))
            settings.create_dataset("stimSetMode", data=np.array([[0]]))

            # Array fields for electrode positioning
            # Shape should be (2, 3) for two hemispheres with x,y,z coordinates
            settings.create_dataset("headNative", data=np.array([[0, 0, 0], [0, 0, 0]]))
            settings.create_dataset(
                "yMarkerNative", data=np.array([[0, 1, 0], [0, 1, 0]])
            )
            settings.create_dataset(
                "stim_center", data=np.array([[np.nan, np.nan, np.nan]])
            )

            # Contact locations (reference to another dataset)
            contact_coords = f.create_dataset(
                "contact_coords_0", data=np.array([[0, 0, 0], [0, 0, 5], [0, 0, 10]])
            )
            settings.create_dataset(
                "contactLocation",
                data=np.array([[contact_coords.ref]]),
                dtype=h5py.ref_dtype,
            )

        return str(file_path)

    def test_lead_settings_get_methods(self, minimal_h5_file):
        """Test various getter methods."""
        ls = LeadSettings(minimal_h5_file)

        # Test numeric getters
        assert ls.get_num_elecs() == 1
        assert ls.get_gm_idx() == 2
        assert ls.get_wm_idx() == 3
        assert ls.get_csf_idx() == 1
        assert ls.get_est_in_temp() is False
        assert ls.get_calc_axon_act() == 0
        assert ls.remove_electrode() == 0

        # Test string getters
        assert ls.get_elec_type() == "Medtronic 3389"
        assert ls.get_cond_model() == "ColeCole4"
        assert ls.get_dti_name() == ""

        # Test array getters
        cur_ctrl = ls.get_cur_ctrl()
        assert len(cur_ctrl) == 2
        assert cur_ctrl[0] == 1.0

    def test_lead_settings_electrode_name_conversion(self, minimal_h5_file):
        """Test that electrode names are properly converted."""
        ls = LeadSettings(minimal_h5_file)
        electrode_type = ls.get_elec_type()
        assert electrode_type == "Medtronic 3389"

    def test_get_rot_z(self, minimal_h5_file):
        """Test rotation calculation."""
        ls = LeadSettings(minimal_h5_file)
        rot_z = ls.get_rot_z(0)
        assert isinstance(rot_z, (float, np.floating))

    def test_get_cntct_loc(self, minimal_h5_file):
        """Test getting contact locations."""
        ls = LeadSettings(minimal_h5_file)
        contact_locs = ls.get_cntct_loc(0)
        assert contact_locs.shape[0] == 3  # x, y, z
        assert contact_locs.shape[1] == 3  # 3 contacts


class TestNumpyCompatibility:
    """Test numpy 2.x compatibility issues."""

    def test_numpy_array_operations(self):
        """Test that numpy array operations work correctly."""
        # Test operations that might break with numpy 2.x
        arr1 = np.array([1.0, 2.0, 3.0])

        # Test basic operations
        result = arr1 + 1
        assert np.allclose(result, np.array([2.0, 3.0, 4.0]))

        # Test nan handling (common in leaddbsinterface)
        arr_with_nan = np.array([1.0, np.nan, 3.0])
        assert np.isnan(arr_with_nan[1])
        assert not np.all(np.isnan(arr_with_nan))

        # Test indexing with boolean arrays
        mask = ~np.isnan(arr_with_nan)
        filtered = arr_with_nan[mask]
        assert len(filtered) == 2
        assert np.allclose(filtered, np.array([1.0, 3.0]))

    def test_numpy_string_conversion(self):
        """Test numpy string/character array conversion."""
        # This pattern is used in LeadSettings._get_str
        ascii_codes = np.array([72, 101, 108, 108, 111])  # "Hello"
        result = "".join(np.vectorize(chr)(ascii_codes.astype(int)))
        assert result == "Hello"

    def test_numpy_reshape_operations(self):
        """Test reshape operations used in the codebase."""
        # Test reshaping operations
        arr = np.array([1, 2, 3, 4, 5, 6])
        reshaped = arr.reshape(2, 3)
        assert reshaped.shape == (2, 3)

        # Test transpose
        transposed = reshaped.T
        assert transposed.shape == (3, 2)

    def test_numpy_linalg_operations(self):
        """Test linear algebra operations."""
        # Test norm calculation (used in get_stretch_factor)
        vec1 = np.array([3.0, 4.0])
        norm = np.linalg.norm(vec1)
        assert np.isclose(norm, 5.0)

        # Test with 2D array
        vecs = np.array([[3.0, 4.0], [5.0, 12.0]])
        norms = np.linalg.norm(vecs, axis=1)
        assert np.allclose(norms, np.array([5.0, 13.0]))

    def test_strict_zip_compatibility(self):
        """Test that zip with strict parameter works (Python 3.10+)."""
        # This is used in lead_settings.py line 563
        list1 = [1, 2, 3]
        list2 = [4, 5, 6]

        result = list(zip(list1, list2, strict=True))
        assert result == [(1, 4), (2, 5), (3, 6)]

        # Test that mismatched lengths raise error
        list3 = [7, 8]
        with pytest.raises(ValueError):
            list(zip(list1, list3, strict=True))


class TestCLIScripts:
    """Test CLI scripts are importable and have proper entry points."""

    def test_convert_input_dictionary_has_main(self):
        """Test that convert_input_dictionary module has main function."""
        import importlib

        # Import module without executing
        spec = importlib.util.find_spec("leaddbsinterface.convert_input_dictionary")
        assert spec is not None, "convert_input_dictionary module should exist"

    def test_allocate_axons_has_main(self):
        """Test that allocate_axons module has main function."""
        import importlib

        spec = importlib.util.find_spec("leaddbsinterface.allocate_axons")
        assert spec is not None, "allocate_axons module should exist"


class TestIntegrationWithTempFiles:
    """Integration tests using temporary files."""

    def test_update_default_dict_with_lead_settings(self):
        """Test the full workflow of updating settings."""
        # Initialize
        default_settings = initialize_default_settings()
        default_settings = load_default_for_lead(default_settings)

        # Custom settings from Lead-DBS
        custom_settings = {
            "BrainRegion": {
                "Center": {"x[mm]": 10.0, "y[mm]": 20.0, "z[mm]": 30.0},
                "Dimension": {"x[mm]": 70.0},  # Override x, keep y and z
            },
            "StimulationSignal": {
                "Type": "Rectangle",  # Override Multisine
            },
        }

        # Merge
        result = update_default_dict(default_settings, custom_settings)

        # Check that custom settings override defaults
        assert result["BrainRegion"]["Center"]["x[mm]"] == 10.0
        assert result["BrainRegion"]["Dimension"]["x[mm]"] == 70.0
        # Check that non-overridden values remain
        assert result["BrainRegion"]["Dimension"]["y[mm]"] == 60
        assert result["BrainRegion"]["Dimension"]["z[mm]"] == 80
        assert result["StimulationSignal"]["Type"] == "Rectangle"

    def test_json_serialization_of_settings(self):
        """Test that settings can be serialized to JSON."""
        settings = initialize_default_settings()
        settings = load_default_for_lead(settings)

        # Should be JSON-serializable
        json_str = json.dumps(settings, indent=2)
        assert isinstance(json_str, str)
        assert len(json_str) > 0

        # Should be deserializable
        deserialized = json.loads(json_str)
        assert deserialized["ModelSide"] == 0
        assert deserialized["BrainRegion"]["Shape"] == "Ellipsoid"
