"""Tests for ossdbs.axon_processing.axon_models module.

Tests cover the AxonMorphology classes and AxonModels class.
"""

import json

import numpy as np
import pytest
# do not test if neuron is not installed
neuron = pytest.importorskip("neuron")

from ossdbs.axon_processing.axon_models import (
    AxonModels,
    AxonMorphologyMcNeal1976,
    AxonMorphologyMRG2002,
)


class TestAxonMorphologyMRG2002:
    """Tests for AxonMorphologyMRG2002 class."""

    def test_init_default_not_downsampled(self):
        """Test default initialization is not downsampled."""
        morph = AxonMorphologyMRG2002()
        assert morph.downsampled is False

    def test_init_downsampled(self):
        """Test initialization with downsampled=True."""
        morph = AxonMorphologyMRG2002(downsampled=True)
        assert morph.downsampled is True

    def test_fiber_diam_positive_validation(self):
        """Test that negative fiber diameter raises ValueError."""
        morph = AxonMorphologyMRG2002()
        with pytest.raises(ValueError, match="greater than zero"):
            morph.fiber_diam = -1.0

    def test_fiber_diam_setter_positive(self):
        """Test that positive fiber diameter is accepted."""
        morph = AxonMorphologyMRG2002()
        morph.fiber_diam = 5.7
        assert morph.fiber_diam == 5.7

    def test_n_ranvier_must_be_odd(self):
        """Test that n_Ranvier is adjusted to be odd."""
        morph = AxonMorphologyMRG2002()
        morph.update_axon_morphology(5.7, n_Ranvier=22)
        # 22 is even, should become 21
        assert morph.n_Ranvier == 21

    def test_n_ranvier_odd_stays_same(self):
        """Test that odd n_Ranvier stays the same."""
        morph = AxonMorphologyMRG2002()
        morph.update_axon_morphology(5.7, n_Ranvier=21)
        assert morph.n_Ranvier == 21

    def test_update_axon_morphology_large_fiber(self):
        """Test morphology update for large fiber (>= 5.7 um)."""
        morph = AxonMorphologyMRG2002()
        morph.update_axon_morphology(5.7, n_Ranvier=21)

        assert morph.fiber_diam == 5.7
        assert morph.n_Ranvier == 21
        # Large fiber has 11 compartments per internode
        assert morph.n_comp == 11
        # n_segments = (n_Ranvier - 1) * n_comp + 1
        assert morph.n_segments == (21 - 1) * 11 + 1  # 221

    def test_update_axon_morphology_small_fiber(self):
        """Test morphology update for small fiber (< 5.7 um)."""
        morph = AxonMorphologyMRG2002()
        morph.update_axon_morphology(2.0, n_Ranvier=21)

        assert morph.fiber_diam == 2.0
        assert morph.n_Ranvier == 21
        # Small fiber has 8 compartments per internode
        assert morph.n_comp == 8
        # n_segments = (n_Ranvier - 1) * n_comp + 1
        assert morph.n_segments == (21 - 1) * 8 + 1  # 161

    def test_update_axon_morphology_downsampled_large(self):
        """Test downsampled morphology for large fiber."""
        morph = AxonMorphologyMRG2002(downsampled=True)
        morph.update_axon_morphology(5.7, n_Ranvier=21)

        # Downsampled large fiber has 3 compartments per internode
        assert morph.n_comp == 3
        assert morph.n_segments == (21 - 1) * 3 + 1  # 61

    def test_update_axon_morphology_downsampled_small(self):
        """Test downsampled morphology for small fiber."""
        morph = AxonMorphologyMRG2002(downsampled=True)
        morph.update_axon_morphology(2.0, n_Ranvier=21)

        # Downsampled small fiber has 2 compartments per internode
        assert morph.n_comp == 2
        assert morph.n_segments == (21 - 1) * 2 + 1  # 41

    def test_get_n_comp_large_downsampled(self):
        """Test get_n_comp for large fiber downsampled."""
        morph = AxonMorphologyMRG2002()
        morph.update_axon_morphology(5.7, n_Ranvier=21)
        assert morph.get_n_comp(downsampled=True) == 3
        assert morph.get_n_comp(downsampled=False) == 11

    def test_get_n_comp_small_downsampled(self):
        """Test get_n_comp for small fiber downsampled."""
        morph = AxonMorphologyMRG2002()
        morph.update_axon_morphology(2.0, n_Ranvier=21)
        assert morph.get_n_comp(downsampled=True) == 2
        assert morph.get_n_comp(downsampled=False) == 8

    def test_get_n_segments_downsampled_vs_full(self):
        """Test get_n_segments returns different values for downsampled."""
        morph = AxonMorphologyMRG2002()
        morph.update_axon_morphology(5.7, n_Ranvier=21)

        n_seg_full = morph.get_n_segments(downsampled=False)
        n_seg_ds = morph.get_n_segments(downsampled=True)

        assert n_seg_full == 221
        assert n_seg_ds == 61
        assert n_seg_full > n_seg_ds


class TestGetLocalCompartmentCoords:
    """Tests for get_local_compartment_coords method."""

    def test_large_fiber_full_resolution(self):
        """Test local coords for large fiber at full resolution."""
        morph = AxonMorphologyMRG2002(downsampled=False)
        morph.update_axon_morphology(5.7, n_Ranvier=21)

        loc_coords = morph.get_local_compartment_coords()

        # Should have n_comp - 1 coordinates (10 for large fiber)
        assert loc_coords.shape == (10,)
        # All coordinates should be non-negative (distance from node)
        assert np.all(loc_coords >= 0)
        # Last coordinate should be the largest (farthest from node)
        assert loc_coords[-1] >= loc_coords[0]

    def test_large_fiber_downsampled(self):
        """Test local coords for large fiber downsampled."""
        morph = AxonMorphologyMRG2002(downsampled=True)
        morph.update_axon_morphology(5.7, n_Ranvier=21)

        loc_coords = morph.get_local_compartment_coords()

        # Should have n_comp - 1 = 2 coordinates
        assert loc_coords.shape == (2,)
        # At least the last coordinate should be positive
        assert loc_coords[-1] > 0
        assert loc_coords[1] > loc_coords[0]

    def test_small_fiber_full_resolution(self):
        """Test local coords for small fiber at full resolution."""
        morph = AxonMorphologyMRG2002(downsampled=False)
        morph.update_axon_morphology(2.0, n_Ranvier=21)

        loc_coords = morph.get_local_compartment_coords()

        # Should have n_comp - 1 = 7 coordinates
        assert loc_coords.shape == (7,)
        assert np.all(loc_coords >= 0)

    def test_small_fiber_downsampled(self):
        """Test local coords for small fiber downsampled."""
        morph = AxonMorphologyMRG2002(downsampled=True)
        morph.update_axon_morphology(2.0, n_Ranvier=21)

        loc_coords = morph.get_local_compartment_coords()

        # Should have n_comp - 1 = 1 coordinate
        assert loc_coords.shape == (1,)
        assert loc_coords[0] > 0


class TestAxonMorphologyMcNeal1976:
    """Tests for AxonMorphologyMcNeal1976 class."""

    def test_init_default_not_downsampled(self):
        """Test default initialization is not downsampled."""
        morph = AxonMorphologyMcNeal1976()
        assert morph.downsampled is False

    def test_downsampled_raises_not_implemented(self):
        """Test that setting downsampled=True raises NotImplementedError."""
        with pytest.raises(NotImplementedError, match="Downsampled McNeal1976"):
            AxonMorphologyMcNeal1976(downsampled=True)

    def test_n_comp_is_fixed_at_2(self):
        """Test that n_comp is always 2 for McNeal1976."""
        morph = AxonMorphologyMcNeal1976()
        assert morph.n_comp == 2

    def test_node_step_depends_on_fiber_diam(self):
        """Test that node_step is 0.2 * fiber_diam."""
        morph = AxonMorphologyMcNeal1976()
        morph.update_axon_morphology(10.0, n_Ranvier=21)
        # node_step = fiber_diam * 0.2
        assert morph.node_step == 10.0 * 0.2  # 2.0

    def test_update_axon_morphology(self):
        """Test morphology update for McNeal1976."""
        morph = AxonMorphologyMcNeal1976()
        morph.update_axon_morphology(10.0, n_Ranvier=21)

        assert morph.fiber_diam == 10.0
        assert morph.n_Ranvier == 21
        assert morph.n_comp == 2
        # n_segments = (n_Ranvier - 1) * n_comp + 1
        assert morph.n_segments == (21 - 1) * 2 + 1  # 41

    def test_get_local_compartment_coords(self):
        """Test local coords for McNeal1976."""
        morph = AxonMorphologyMcNeal1976()
        morph.update_axon_morphology(10.0, n_Ranvier=21)

        loc_coords = morph.get_local_compartment_coords()

        # Should have n_comp - 1 = 1 coordinate
        assert loc_coords.shape == (1,)
        # loc_coords[0] = node_step * 0.5 = 2.0 * 0.5 = 1.0
        assert loc_coords[0] == morph.node_step * 0.5

    def test_get_local_compartment_coords_downsampled_raises(self):
        """Test that get_local_compartment_coords with downsampled raises."""
        morph = AxonMorphologyMcNeal1976()
        morph.update_axon_morphology(10.0, n_Ranvier=21)
        # Manually set downsampled to trigger error
        morph._downsampled = True

        with pytest.raises(NotImplementedError, match="Downsampled McNeal1976"):
            morph.get_local_compartment_coords()


class TestAxonModelSetter:
    """Tests for AxonModels.axon_model property setter."""

    def test_axon_model_setter_valid_mrg2002(self, tmp_path):
        """Test that valid MRG2002 model is accepted."""
        # Create minimal json file - pathway_mat_file needs full path with separator
        config = {
            "pathway_mat_file": [str(tmp_path / "test.mat")],
            "axon_diams_all": [5.7],
            "axon_lengths_all": [20.0],
            "centering_coordinates": [[0, 0, 0]],
            "axon_model": "MRG2002",
            "combined_h5_file": str(tmp_path / "test.h5"),
        }
        config_file = tmp_path / "config.json"
        with open(config_file, "w") as f:
            json.dump(config, f)

        model = AxonModels(str(tmp_path), 0, str(config_file))
        assert model.axon_model == "MRG2002"

    def test_axon_model_setter_valid_mrg2002_ds(self, tmp_path):
        """Test that valid MRG2002_DS model is accepted."""
        config = {
            "pathway_mat_file": [str(tmp_path / "test.mat")],
            "axon_diams_all": [5.7],
            "axon_lengths_all": [20.0],
            "centering_coordinates": [[0, 0, 0]],
            "axon_model": "MRG2002_DS",
            "combined_h5_file": str(tmp_path / "test.h5"),
        }
        config_file = tmp_path / "config.json"
        with open(config_file, "w") as f:
            json.dump(config, f)

        model = AxonModels(str(tmp_path), 0, str(config_file))
        assert model.axon_model == "MRG2002_DS"

    def test_axon_model_setter_valid_mcneal1976(self, tmp_path):
        """Test that valid McNeal1976 model is accepted."""
        config = {
            "pathway_mat_file": [str(tmp_path / "test.mat")],
            "axon_diams_all": [5.7],
            "axon_lengths_all": [20.0],
            "centering_coordinates": [[0, 0, 0]],
            "axon_model": "McNeal1976",
            "combined_h5_file": str(tmp_path / "test.h5"),
        }
        config_file = tmp_path / "config.json"
        with open(config_file, "w") as f:
            json.dump(config, f)

        model = AxonModels(str(tmp_path), 0, str(config_file))
        assert model.axon_model == "McNeal1976"

    def test_axon_model_setter_invalid_model_raises(self, tmp_path):
        """Test that invalid model raises ValueError."""
        config = {
            "pathway_mat_file": [str(tmp_path / "test.mat")],
            "axon_diams_all": [5.7],
            "axon_lengths_all": [20.0],
            "centering_coordinates": [[0, 0, 0]],
            "axon_model": "InvalidModel",
            "combined_h5_file": str(tmp_path / "test.h5"),
        }
        config_file = tmp_path / "config.json"
        with open(config_file, "w") as f:
            json.dump(config, f)

        with pytest.raises(ValueError, match="NEURON model is not valid"):
            AxonModels(str(tmp_path), 0, str(config_file))


class TestAxonModelsInit:
    """Tests for AxonModels initialization."""

    def test_stim_dir_must_exist(self, tmp_path):
        """Test that stim_dir must be an existing directory."""
        with pytest.raises(ValueError, match="does not exist"):
            AxonModels("/nonexistent/path", 0, "config.json")

    def test_hemis_idx_must_be_0_or_1(self, tmp_path):
        """Test that hemis_idx must be 0 or 1."""
        config = {
            "pathway_mat_file": [str(tmp_path / "test.mat")],
            "axon_diams_all": [5.7],
            "axon_lengths_all": [20.0],
            "centering_coordinates": [[0, 0, 0]],
            "axon_model": "MRG2002",
            "combined_h5_file": str(tmp_path / "test.h5"),
        }
        config_file = tmp_path / "config.json"
        with open(config_file, "w") as f:
            json.dump(config, f)

        with pytest.raises(ValueError, match="0 or 1"):
            AxonModels(str(tmp_path), 2, str(config_file))

    def test_hemis_idx_0_sets_rh_folder(self, tmp_path):
        """Test that hemis_idx=0 sets right hemisphere folder."""
        config = {
            "pathway_mat_file": [str(tmp_path / "test.mat")],
            "axon_diams_all": [5.7],
            "axon_lengths_all": [20.0],
            "centering_coordinates": [[0, 0, 0]],
            "axon_model": "MRG2002",
            "combined_h5_file": str(tmp_path / "test.h5"),
        }
        config_file = tmp_path / "config.json"
        with open(config_file, "w") as f:
            json.dump(config, f)

        model = AxonModels(str(tmp_path), 0, str(config_file))
        assert model.oss_sim_folder == "OSS_sim_files_rh"

    def test_hemis_idx_1_sets_lh_folder(self, tmp_path):
        """Test that hemis_idx=1 sets left hemisphere folder."""
        config = {
            "pathway_mat_file": [str(tmp_path / "test.mat")],
            "axon_diams_all": [5.7],
            "axon_lengths_all": [20.0],
            "centering_coordinates": [[0, 0, 0]],
            "axon_model": "MRG2002",
            "combined_h5_file": str(tmp_path / "test.h5"),
        }
        config_file = tmp_path / "config.json"
        with open(config_file, "w") as f:
            json.dump(config, f)

        model = AxonModels(str(tmp_path), 1, str(config_file))
        assert model.oss_sim_folder == "OSS_sim_files_lh"

    def test_unsupported_file_format_raises(self, tmp_path):
        """Test that unsupported file format raises NotImplementedError."""
        txt_file = tmp_path / "config.txt"
        txt_file.write_text("invalid content")

        with pytest.raises(NotImplementedError, match="Unsupported input format"):
            AxonModels(str(tmp_path), 0, str(txt_file))


class TestCombinedH5File:
    """Tests for combined_h5_file property."""

    def test_h5_extension_added_if_missing(self, tmp_path):
        """Test that .h5 extension is added if missing."""
        config = {
            "pathway_mat_file": [str(tmp_path / "test.mat")],
            "axon_diams_all": [5.7],
            "axon_lengths_all": [20.0],
            "centering_coordinates": [[0, 0, 0]],
            "axon_model": "MRG2002",
            "combined_h5_file": str(tmp_path / "test"),  # No .h5
        }
        config_file = tmp_path / "config.json"
        with open(config_file, "w") as f:
            json.dump(config, f)

        model = AxonModels(str(tmp_path), 0, str(config_file))
        assert model.combined_h5_file.endswith(".h5")

    def test_h5_extension_not_duplicated(self, tmp_path):
        """Test that .h5 extension is not duplicated if present."""
        config = {
            "pathway_mat_file": [str(tmp_path / "test.mat")],
            "axon_diams_all": [5.7],
            "axon_lengths_all": [20.0],
            "centering_coordinates": [[0, 0, 0]],
            "axon_model": "MRG2002",
            "combined_h5_file": str(tmp_path / "test.h5"),
        }
        config_file = tmp_path / "config.json"
        with open(config_file, "w") as f:
            json.dump(config, f)

        model = AxonModels(str(tmp_path), 0, str(config_file))
        assert model.combined_h5_file == str(tmp_path / "test.h5")
        assert not model.combined_h5_file.endswith(".h5.h5")


class TestSelectAxonMorphologyModel:
    """Tests for _select_axon_morphology_model method."""

    def test_selects_mrg2002_for_mrg2002(self, tmp_path):
        """Test that MRG2002 model type selects AxonMorphologyMRG2002."""
        config = {
            "pathway_mat_file": [str(tmp_path / "test.mat")],
            "axon_diams_all": [5.7],
            "axon_lengths_all": [20.0],
            "centering_coordinates": [[0, 0, 0]],
            "axon_model": "MRG2002",
            "combined_h5_file": str(tmp_path / "test.h5"),
        }
        config_file = tmp_path / "config.json"
        with open(config_file, "w") as f:
            json.dump(config, f)

        model = AxonModels(str(tmp_path), 0, str(config_file))
        morph = model._select_axon_morphology_model()

        assert isinstance(morph, AxonMorphologyMRG2002)
        assert morph.downsampled is False

    def test_selects_mrg2002_ds_for_mrg2002_ds(self, tmp_path):
        """Test that MRG2002_DS model type selects downsampled MRG2002."""
        config = {
            "pathway_mat_file": [str(tmp_path / "test.mat")],
            "axon_diams_all": [5.7],
            "axon_lengths_all": [20.0],
            "centering_coordinates": [[0, 0, 0]],
            "axon_model": "MRG2002_DS",
            "combined_h5_file": str(tmp_path / "test.h5"),
        }
        config_file = tmp_path / "config.json"
        with open(config_file, "w") as f:
            json.dump(config, f)

        model = AxonModels(str(tmp_path), 0, str(config_file))
        morph = model._select_axon_morphology_model()

        assert isinstance(morph, AxonMorphologyMRG2002)
        assert morph.downsampled is True

    def test_selects_mcneal1976_for_mcneal1976(self, tmp_path):
        """Test that McNeal1976 model type selects AxonMorphologyMcNeal1976."""
        config = {
            "pathway_mat_file": [str(tmp_path / "test.mat")],
            "axon_diams_all": [5.7],
            "axon_lengths_all": [20.0],
            "centering_coordinates": [[0, 0, 0]],
            "axon_model": "McNeal1976",
            "combined_h5_file": str(tmp_path / "test.h5"),
        }
        config_file = tmp_path / "config.json"
        with open(config_file, "w") as f:
            json.dump(config, f)

        model = AxonModels(str(tmp_path), 0, str(config_file))
        morph = model._select_axon_morphology_model()

        assert isinstance(morph, AxonMorphologyMcNeal1976)
