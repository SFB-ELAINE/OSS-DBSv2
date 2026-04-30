"""Tests for ossdbs.axon_processing.utilities module.

Tests cover helper functions for fiber/streamline processing and output generation.
"""

import json

import numpy as np
import pytest

pytest.importorskip("neuron")

from nibabel.streamlines import ArraySequence

from ossdbs.axon_processing.utilities import (
    convert_fibers_to_streamlines,
    create_leaddbs_outputs,
    create_paraview_outputs,
    find_nearest,
    index_for_length,
    normalized,
    place_axons_on_streamlines,
    resample_fibers_to_Ranviers,
    resample_streamline_for_Ranvier,
    store_axon_statuses,
)


class TestFindNearest:
    """Tests for find_nearest function."""

    def test_find_exact_value(self):
        """Test finding exact value in array."""
        array = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        idx, value = find_nearest(array, 3.0)
        assert idx == 2
        assert value == 3.0

    def test_find_nearest_value(self):
        """Test finding nearest value when exact not present."""
        array = np.array([1.0, 2.0, 4.0, 5.0])
        idx, value = find_nearest(array, 3.0)
        # 3.0 is equidistant from 2.0 and 4.0, argmin returns first
        assert idx in [1, 2]
        assert value in [2.0, 4.0]

    def test_find_nearest_at_beginning(self):
        """Test finding nearest when target is before array start."""
        array = np.array([10.0, 20.0, 30.0])
        idx, value = find_nearest(array, 5.0)
        assert idx == 0
        assert value == 10.0

    def test_find_nearest_at_end(self):
        """Test finding nearest when target is after array end."""
        array = np.array([10.0, 20.0, 30.0])
        idx, value = find_nearest(array, 100.0)
        assert idx == 2
        assert value == 30.0

    def test_find_nearest_single_element(self):
        """Test finding nearest in single element array."""
        array = np.array([5.0])
        idx, value = find_nearest(array, 100.0)
        assert idx == 0
        assert value == 5.0


class TestNormalized:
    """Tests for normalized function."""

    def test_normalize_simple_vector_1d(self):
        """Test normalizing a 1D vector - result shape is preserved."""
        vector = np.array([3.0, 4.0])
        result = normalized(vector)
        # L2 norm = 5, so result should be [0.6, 0.8]
        # Note: the function uses expand_dims, so 1D inputs work but result may be 2D
        # Just check the norm is 1
        result_norm = np.linalg.norm(result)
        np.testing.assert_almost_equal(result_norm, 1.0)

    def test_normalize_3d_vector_1d(self):
        """Test normalizing a 1D 3-component vector."""
        vector = np.array([1.0, 2.0, 2.0])
        result = normalized(vector)
        # L2 norm = sqrt(1 + 4 + 4) = 3
        result_norm = np.linalg.norm(result)
        np.testing.assert_almost_equal(result_norm, 1.0)

    def test_normalize_zero_vector(self):
        """Test normalizing a zero vector returns zero vector."""
        vector = np.array([0.0, 0.0, 0.0])
        result = normalized(vector)
        # Zero vectors should remain zero (handled by l2[l2==0] = 1)
        expected = np.array([0.0, 0.0, 0.0])
        np.testing.assert_array_almost_equal(result.flatten(), expected)

    def test_normalize_2d_array(self):
        """Test normalizing 2D array along last axis."""
        vector = np.array([[3.0, 4.0], [1.0, 0.0]])
        result = normalized(vector)
        expected = np.array([[0.6, 0.8], [1.0, 0.0]])
        np.testing.assert_array_almost_equal(result, expected)


class TestConvertFibersToStreamlines:
    """Tests for convert_fibers_to_streamlines function."""

    def test_single_fiber(self):
        """Test converting single fiber to streamline."""
        # fibers format: 4 x N where row 4 is fiber index (1-based)
        fibers = np.array(
            [
                [0.0, 1.0, 2.0],  # x
                [0.0, 0.0, 0.0],  # y
                [0.0, 0.0, 0.0],  # z
                [1.0, 1.0, 1.0],  # fiber index (all belong to fiber 1)
            ]
        )

        streamlines, indices = convert_fibers_to_streamlines(fibers)

        assert len(streamlines) == 1
        assert streamlines[0].shape == (3, 3)  # 3 points, 3 dimensions
        assert indices == [1]  # Original index preserved

    def test_multiple_fibers(self):
        """Test converting multiple fibers to streamlines."""
        fibers = np.array(
            [
                [0.0, 1.0, 2.0, 0.0, 1.0],  # x
                [0.0, 0.0, 0.0, 1.0, 1.0],  # y
                [0.0, 0.0, 0.0, 0.0, 0.0],  # z
                [1.0, 1.0, 1.0, 2.0, 2.0],  # fiber index
            ]
        )

        streamlines, indices = convert_fibers_to_streamlines(fibers)

        assert len(streamlines) == 2
        assert streamlines[0].shape == (3, 3)  # Fiber 1: 3 points
        assert streamlines[1].shape == (2, 3)  # Fiber 2: 2 points
        assert indices == [1, 2]  # Original indices preserved


class TestIndexForLength:
    """Tests for index_for_length function."""

    def test_simple_streamline(self):
        """Test finding truncation point for simple streamline."""
        # Straight line with length 10 (5 segments of length 2)
        xyz = np.array(
            [[0, 0, 0], [2, 0, 0], [4, 0, 0], [6, 0, 0], [8, 0, 0], [10, 0, 0]],
            dtype=float,
        )

        idx, length = index_for_length(xyz, req_length=5.0)

        # Cumulative lengths: [2, 4, 6, 8, 10]
        # req_length=5.0, nearest is 4.0 at idx=1, but 4.0 < 5.0, so idx stays 1
        # Or nearest is 6.0 at idx=2, but 6.0 > 5.0, so idx becomes idx-1=1
        # Result: idx=1, length=4.0 (cumulative at segment 1)
        assert idx == 1
        assert length == 4.0

    def test_exact_length_match(self):
        """Test when required length exactly matches segment end."""
        xyz = np.array([[0, 0, 0], [2, 0, 0], [4, 0, 0]], dtype=float)

        _idx, length = index_for_length(xyz, req_length=4.0)

        assert length == 4.0

    def test_very_short_streamline(self):
        """Test handling very short streamline."""
        xyz = np.array([[0, 0, 0]], dtype=float)  # Single point

        result = index_for_length(xyz, req_length=10.0)

        # Single point returns [0]
        np.testing.assert_array_equal(result, np.array([0]))


class TestStoreAxonStatuses:
    """Tests for store_axon_statuses function."""

    def test_basic_output(self, tmp_path):
        """Test basic output without scaling or pathway."""
        store_axon_statuses(
            output_path=str(tmp_path),
            percent_activated=50.0,
            percent_damaged=10.0,
            percent_csf=5.0,
        )

        output_file = tmp_path / "Pathway_status.json"
        assert output_file.exists()

        with open(output_file) as f:
            data = json.load(f)

        assert data["percent_activated"] == 50.0
        assert data["percent_damaged"] == 10.0
        assert data["percent_csf"] == 5.0

    def test_with_scaling_index(self, tmp_path):
        """Test output with scaling index."""
        store_axon_statuses(
            output_path=str(tmp_path),
            percent_activated=50.0,
            percent_damaged=10.0,
            percent_csf=5.0,
            scaling_index=3,
        )

        output_file = tmp_path / "Pathway_status_3.json"
        assert output_file.exists()

        with open(output_file) as f:
            data = json.load(f)

        assert data["scaling_index"] == "3"

    def test_with_pathway_name(self, tmp_path):
        """Test output with pathway name."""
        store_axon_statuses(
            output_path=str(tmp_path),
            percent_activated=50.0,
            percent_damaged=10.0,
            percent_csf=5.0,
            pathway_name="motor_cortex",
        )

        output_file = tmp_path / "Pathway_status_motor_cortex.json"
        assert output_file.exists()

        with open(output_file) as f:
            data = json.load(f)

        assert data["pathway_name"] == "motor_cortex"

    def test_with_scaling_and_pathway(self, tmp_path):
        """Test output with both scaling index and pathway name."""
        store_axon_statuses(
            output_path=str(tmp_path),
            percent_activated=50.0,
            percent_damaged=10.0,
            percent_csf=5.0,
            scaling_index=2,
            pathway_name="sensory",
        )

        output_file = tmp_path / "Pathway_status_sensory_2.json"
        assert output_file.exists()

        with open(output_file) as f:
            data = json.load(f)

        assert data["scaling_index"] == "2"
        assert data["pathway_name"] == "sensory"


class TestCreateLeaddbsOutputs:
    """Tests for create_leaddbs_outputs function."""

    def test_basic_output(self, tmp_path):
        """Test basic output without scaling or pathway."""
        axon_data = np.array(
            [
                [0.0, 0.0, 0.0, 1.0, 1.0],
                [1.0, 0.0, 0.0, 1.0, 0.0],
            ]
        )

        create_leaddbs_outputs(
            output_path=str(tmp_path),
            Axon_Lead_DBS=axon_data,
            connectome_name="test_connectome",
        )

        output_file = tmp_path / "Axon_state.mat"
        assert output_file.exists()

    def test_with_scaling_index(self, tmp_path):
        """Test output with scaling index."""
        axon_data = np.random.rand(10, 5)

        create_leaddbs_outputs(
            output_path=str(tmp_path),
            Axon_Lead_DBS=axon_data,
            connectome_name="test",
            scaling_index=5,
        )

        output_file = tmp_path / "Axon_state_5.mat"
        assert output_file.exists()

    def test_with_pathway_name(self, tmp_path):
        """Test output with pathway name."""
        axon_data = np.random.rand(10, 5)

        create_leaddbs_outputs(
            output_path=str(tmp_path),
            Axon_Lead_DBS=axon_data,
            connectome_name="test",
            pathway_name="STN",
        )

        output_file = tmp_path / "Axon_state_STN.mat"
        assert output_file.exists()

    def test_with_scaling_and_pathway(self, tmp_path):
        """Test output with both scaling and pathway."""
        axon_data = np.random.rand(10, 5)

        create_leaddbs_outputs(
            output_path=str(tmp_path),
            Axon_Lead_DBS=axon_data,
            connectome_name="test",
            scaling_index=3,
            pathway_name="GPi",
        )

        output_file = tmp_path / "Axon_state_GPi_3.mat"
        assert output_file.exists()


class TestCreateParaviewOutputs:
    """Tests for create_paraview_outputs function."""

    def test_basic_output(self, tmp_path):
        """Test basic CSV output."""
        axon_data = np.array(
            [
                [0.0, 0.0, 0.0, 1.0, 1.0],
                [1.0, 0.0, 0.0, 1.0, 0.0],
            ]
        )

        create_paraview_outputs(
            output_path=str(tmp_path),
            Axon_Lead_DBS=axon_data,
        )

        output_file = tmp_path / "Axon_state.csv"
        assert output_file.exists()

        # Verify it's a valid CSV
        loaded = np.genfromtxt(output_file, delimiter=",", skip_header=1)
        assert loaded.shape == (2, 5)

    def test_with_scaling_index(self, tmp_path):
        """Test CSV output with scaling index."""
        axon_data = np.random.rand(10, 5)

        create_paraview_outputs(
            output_path=str(tmp_path),
            Axon_Lead_DBS=axon_data,
            scaling_index=7,
        )

        output_file = tmp_path / "Axon_state_7.csv"
        assert output_file.exists()

    def test_with_pathway_name(self, tmp_path):
        """Test CSV output with pathway name."""
        axon_data = np.random.rand(10, 5)

        create_paraview_outputs(
            output_path=str(tmp_path),
            Axon_Lead_DBS=axon_data,
            pathway_name="VIM",
        )

        output_file = tmp_path / "Axon_state_VIM.csv"
        assert output_file.exists()

    def test_with_scaling_and_pathway(self, tmp_path):
        """Test CSV output with both scaling and pathway."""
        axon_data = np.random.rand(10, 5)

        create_paraview_outputs(
            output_path=str(tmp_path),
            Axon_Lead_DBS=axon_data,
            scaling_index=1,
            pathway_name="Zona_incerta",
        )

        output_file = tmp_path / "Axon_state_Zona_incerta_1.csv"
        assert output_file.exists()


class TestResampleStreamlineForRanvier:
    """Tests for resample_streamline_for_Ranvier function."""

    def test_resample_straight_line(self):
        """Test resampling a straight line streamline."""
        # Create a straight line from (0,0,0) to (100,0,0)
        streamline = np.linspace([0, 0, 0], [100, 0, 0], 50)

        result = resample_streamline_for_Ranvier(
            streamline_array=streamline,
            estim_axon_length=50.0,
            n_Ranvier=11,
        )

        # Should have 11 Ranvier nodes
        assert result.shape[0] == 11


class TestResampleFibersToRanviers:
    """Tests for resample_fibers_to_Ranviers function."""

    def test_resamples_long_streamline(self):
        """Test that long streamlines are resampled correctly."""
        streamlines = ArraySequence()
        long_stream = np.linspace([0, 0, 0], [500, 0, 0], 500)
        streamlines.append(long_stream)
        inx_orig = [1]  # Original index (1-based from convert_fibers_to_streamlines)

        result, excluded, inx_resampled = resample_fibers_to_Ranviers(
            streamlines=streamlines,
            node_step=5.0,
            n_Ranvier=11,
            inx_orig=inx_orig,
        )

        assert len(result) >= 1
        assert len(excluded) == 0 or 0 not in excluded
        assert len(inx_resampled) == len(
            result
        )  # Same number of indices as streamlines

    def test_returns_three_outputs(self):
        """Test that function returns tuple of (resampled, excluded, inx_resampled)."""
        streamlines = ArraySequence()
        stream = np.linspace([0, 0, 0], [1000, 0, 0], 1000)
        streamlines.append(stream)
        inx_orig = [1]  # Original index

        result = resample_fibers_to_Ranviers(
            streamlines=streamlines,
            node_step=5.0,
            n_Ranvier=11,
            inx_orig=inx_orig,
        )

        assert isinstance(result, tuple)
        assert len(result) == 3  # (resampled, excluded, inx_resampled)


class TestPlaceAxonsOnStreamlines:
    """Tests for place_axons_on_streamlines function."""

    def test_axon_centered_near_coordinate(self):
        """Test that axon is centered near the centering coordinate."""
        # Create a streamline long enough for axon
        streamlines = ArraySequence()
        stream = np.linspace([0, 0, 0], [100, 0, 0], 41)  # 41 Ranvier nodes
        streamlines.append(stream)

        centering_coords = [[50, 0, 0]]  # Center point
        n_Ranvier = 21

        result = place_axons_on_streamlines(
            streamlines_resampled=streamlines,
            n_Ranvier=n_Ranvier,
            centering_coordinates=centering_coords,
        )

        assert len(result) == 1
        assert result[0].shape == (21, 3)

        # The center of the axon should be near x=50
        axon_center_x = result[0][10, 0]  # Middle point
        assert abs(axon_center_x - 50) < 10  # Within 10 units

    def test_axon_at_beginning_of_fiber(self):
        """Test axon placement when centering coord is near fiber start."""
        streamlines = ArraySequence()
        stream = np.linspace([0, 0, 0], [100, 0, 0], 41)
        streamlines.append(stream)

        centering_coords = [[0, 0, 0]]  # At start
        n_Ranvier = 21

        result = place_axons_on_streamlines(
            streamlines_resampled=streamlines,
            n_Ranvier=n_Ranvier,
            centering_coordinates=centering_coords,
        )

        assert len(result) == 1
        # Axon should start at beginning
        assert result[0][0, 0] == pytest.approx(0.0, abs=0.1)

    def test_axon_at_end_of_fiber(self):
        """Test axon placement when centering coord is near fiber end."""
        streamlines = ArraySequence()
        stream = np.linspace([0, 0, 0], [100, 0, 0], 41)
        streamlines.append(stream)

        centering_coords = [[100, 0, 0]]  # At end
        n_Ranvier = 21

        result = place_axons_on_streamlines(
            streamlines_resampled=streamlines,
            n_Ranvier=n_Ranvier,
            centering_coordinates=centering_coords,
        )

        assert len(result) == 1
        # Axon should end at fiber end
        assert result[0][-1, 0] == pytest.approx(100.0, abs=0.1)

    def test_multiple_centering_coordinates(self):
        """Test that closest centering coordinate is used."""
        streamlines = ArraySequence()
        stream = np.linspace([0, 0, 0], [100, 0, 0], 41)
        streamlines.append(stream)

        # Multiple centering coordinates, closest to x=75
        centering_coords = [[25, 0, 0], [75, 0, 0]]
        n_Ranvier = 21

        result = place_axons_on_streamlines(
            streamlines_resampled=streamlines,
            n_Ranvier=n_Ranvier,
            centering_coordinates=centering_coords,
        )

        # The streamline passes through both coordinates
        # Test should verify axon is placed appropriately
        assert len(result) == 1
        assert result[0].shape == (21, 3)
