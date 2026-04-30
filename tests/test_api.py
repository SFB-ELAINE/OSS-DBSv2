# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from ossdbs.api import (
    PAM_AVAILABLE,
    create_bounding_box,
    generate_point_models,
    generate_signal,
    prepare_stimulation_signal,
)
from ossdbs.stimulation_signals import (
    FrequencyDomainSignal,
    RectangleSignal,
    TrapezoidSignal,
    TriangleSignal,
)


class TestCreateBoundingBox:
    """Tests for create_bounding_box function."""

    def test_basic_bounding_box(self):
        """Test basic bounding box creation."""
        box_parameters = {
            "Dimension": {"x[mm]": 10.0, "y[mm]": 20.0, "z[mm]": 30.0},
            "Center": {"x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 0.0},
        }
        result = create_bounding_box(box_parameters)
        # Start should be -dimension/2, end should be +dimension/2
        assert result.start == (-5.0, -10.0, -15.0)
        assert result.end == (5.0, 10.0, 15.0)

    def test_offset_center(self):
        """Test bounding box with offset center."""
        box_parameters = {
            "Dimension": {"x[mm]": 10.0, "y[mm]": 10.0, "z[mm]": 10.0},
            "Center": {"x[mm]": 5.0, "y[mm]": 5.0, "z[mm]": 5.0},
        }
        result = create_bounding_box(box_parameters)
        assert result.start == (0.0, 0.0, 0.0)
        assert result.end == (10.0, 10.0, 10.0)

    def test_asymmetric_dimensions(self):
        """Test bounding box with asymmetric dimensions."""
        box_parameters = {
            "Dimension": {"x[mm]": 40.0, "y[mm]": 20.0, "z[mm]": 60.0},
            "Center": {"x[mm]": -10.0, "y[mm]": 5.0, "z[mm]": 0.0},
        }
        result = create_bounding_box(box_parameters)
        expected_start = (-10.0 - 20.0, 5.0 - 10.0, 0.0 - 30.0)
        expected_end = (-10.0 + 20.0, 5.0 + 10.0, 0.0 + 30.0)
        np.testing.assert_array_almost_equal(result.start, expected_start)
        np.testing.assert_array_almost_equal(result.end, expected_end)


class TestGenerateSignal:
    """Tests for generate_signal function."""

    @pytest.fixture
    def base_signal_settings(self):
        """Base settings for signal generation."""
        return {
            "Frequency[Hz]": 130.0,
            "PulseWidth[us]": 60.0,
            "InterPulseWidth[us]": 60.0,
            "CounterPulseWidth[us]": 120.0,
            "CounterAmplitude": 0.5,
            "CutoffFrequency": 100000.0,
        }

    @patch("ossdbs.stimulation_signals.signal.plt")
    def test_rectangle_signal(self, mock_plt, base_signal_settings, tmp_path):
        """Test Rectangle signal generation."""
        settings = {
            "StimulationSignal": {**base_signal_settings, "Type": "Rectangle"},
            "OutputPath": str(tmp_path),
        }
        signal = generate_signal(settings)
        assert isinstance(signal, RectangleSignal)
        assert signal.frequency == 130.0

    @patch("ossdbs.stimulation_signals.signal.plt")
    def test_triangle_signal(self, mock_plt, base_signal_settings, tmp_path):
        """Test Triangle signal generation."""
        settings = {
            "StimulationSignal": {**base_signal_settings, "Type": "Triangle"},
            "OutputPath": str(tmp_path),
        }
        signal = generate_signal(settings)
        assert isinstance(signal, TriangleSignal)
        assert signal.frequency == 130.0

    @patch("ossdbs.stimulation_signals.signal.plt")
    def test_trapezoid_signal(self, mock_plt, base_signal_settings, tmp_path):
        """Test Trapezoid signal generation."""
        signal_settings = {
            **base_signal_settings,
            "Type": "Trapezoid",
            "PulseTopWidth[us]": 30.0,
        }
        settings = {
            "StimulationSignal": signal_settings,
            "OutputPath": str(tmp_path),
        }
        signal = generate_signal(settings)
        assert isinstance(signal, TrapezoidSignal)
        assert signal.frequency == 130.0


class TestPrepareStimulationSignal:
    """Tests for prepare_stimulation_signal function."""

    def test_multisine_signal(self):
        """Test Multisine signal preparation."""
        settings = {
            "StimulationSignal": {
                "Type": "Multisine",
                "CurrentControlled": True,
                "ListOfFrequencies": [1000.0, 2000.0, 3000.0],
            }
        }
        result = prepare_stimulation_signal(settings)
        assert isinstance(result, FrequencyDomainSignal)
        assert result.current_controlled is True
        assert len(result.frequencies) == 3

    @patch("ossdbs.stimulation_signals.signal.plt")
    def test_rectangle_signal_preparation(self, mock_plt, tmp_path):
        """Test Rectangle signal preparation to frequency domain."""
        settings = {
            "StimulationSignal": {
                "Type": "Rectangle",
                "CurrentControlled": False,
                "Frequency[Hz]": 130.0,
                "PulseWidth[us]": 60.0,
                "InterPulseWidth[us]": 60.0,
                "CounterPulseWidth[us]": 120.0,
                "CounterAmplitude": 0.5,
                "CutoffFrequency": 10000.0,
                "SpectrumMode": "Full",
            },
            "OutputPath": str(tmp_path),
        }
        result = prepare_stimulation_signal(settings)
        assert isinstance(result, FrequencyDomainSignal)
        assert result.current_controlled is False
        assert result.base_frequency == 130.0

    @patch("ossdbs.stimulation_signals.signal.plt")
    def test_octave_band_approximation(self, mock_plt, tmp_path):
        """Test OctaveBand spectrum mode."""
        settings = {
            "StimulationSignal": {
                "Type": "Rectangle",
                "CurrentControlled": True,
                "Frequency[Hz]": 130.0,
                "PulseWidth[us]": 60.0,
                "InterPulseWidth[us]": 60.0,
                "CounterPulseWidth[us]": 120.0,
                "CounterAmplitude": 0.5,
                "CutoffFrequency": 10000.0,
                "SpectrumMode": "OctaveBand",
            },
            "OutputPath": str(tmp_path),
        }
        result = prepare_stimulation_signal(settings)
        assert isinstance(result, FrequencyDomainSignal)
        assert result.octave_band_approximation is True


class TestGeneratePointModels:
    """Tests for generate_point_models function."""

    def test_no_active_point_models(self):
        """Test with no active point models."""
        settings = {
            "PointModel": {
                "Pathway": {"Active": False},
                "Lattice": {"Active": False},
                "VoxelLattice": {"Active": False},
            }
        }
        result = generate_point_models(settings)
        assert result == []

    def test_lattice_point_model(self):
        """Test Lattice point model generation."""
        settings = {
            "PointModel": {
                "Pathway": {"Active": False},
                "Lattice": {
                    "Active": True,
                    "Shape": {"x": 10, "y": 10, "z": 10},
                    "Center": {"x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 0.0},
                    "Direction": {"x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 1.0},
                    "PointDistance[mm]": 0.5,
                    "CollapseVTA": False,
                    "ExportField": False,
                },
                "VoxelLattice": {"Active": False},
            }
        }
        result = generate_point_models(settings)
        assert len(result) == 1
        from ossdbs.point_analysis import Lattice

        assert isinstance(result[0], Lattice)

    @patch("ossdbs.api.Pathway")
    def test_pathway_point_model(self, mock_pathway):
        """Test Pathway point model generation with mock."""
        mock_pathway.return_value = MagicMock()
        settings = {
            "PointModel": {
                "Pathway": {
                    "Active": True,
                    "FileName": "test_pathway.h5",
                    "ExportField": True,
                },
                "Lattice": {"Active": False},
                "VoxelLattice": {"Active": False},
            }
        }
        result = generate_point_models(settings)
        assert len(result) == 1
        mock_pathway.assert_called_once_with("test_pathway.h5", export_field=True)

    def test_multiple_point_models(self):
        """Test multiple point models active."""
        settings = {
            "PointModel": {
                "Pathway": {"Active": False},
                "Lattice": {
                    "Active": True,
                    "Shape": {"x": 5, "y": 5, "z": 5},
                    "Center": {"x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 0.0},
                    "Direction": {"x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 1.0},
                    "PointDistance[mm]": 1.0,
                    "CollapseVTA": False,
                    "ExportField": False,
                },
                "VoxelLattice": {"Active": False},
            }
        }
        result = generate_point_models(settings)
        # Only Lattice is active and doesn't require external files
        assert len(result) == 1


class TestRunVolumeConductorModel:
    """Tests for run_volume_conductor_model function."""

    @pytest.fixture
    def base_settings(self):
        """Base settings for volume conductor tests."""
        return {
            "OutOfCore": False,
            "ComputeImpedance": False,
            "ExportVTK": False,
            "ExportFrequency": None,
            "ActivationThresholdVTA[V-per-m]": 0.2,
            "Mesh": {
                "AdaptiveMeshRefinement": {"Active": False},
                "MaterialRefinementSteps": 0,
            },
            "PointModel": {
                "Pathway": {"Active": False},
                "Lattice": {"Active": False},
                "VoxelLattice": {"Active": False},
            },
        }

    def test_basic_run(self, base_settings):
        """Test basic volume conductor model run."""
        from ossdbs.api import run_volume_conductor_model

        # Mock volume conductor
        mock_vc = MagicMock()
        mock_vc.run_full_analysis.return_value = {"solve_time": 1.0}

        # Mock frequency domain signal
        mock_signal = MagicMock()

        result = run_volume_conductor_model(base_settings, mock_vc, mock_signal)

        mock_vc.run_full_analysis.assert_called_once()
        assert result == {"solve_time": 1.0}

    def test_with_compute_impedance(self, base_settings):
        """Test with impedance computation enabled."""
        from ossdbs.api import run_volume_conductor_model

        base_settings["ComputeImpedance"] = True

        mock_vc = MagicMock()
        mock_vc.run_full_analysis.return_value = {}
        mock_signal = MagicMock()

        run_volume_conductor_model(base_settings, mock_vc, mock_signal)

        call_args = mock_vc.run_full_analysis.call_args
        # Second positional arg is compute_impedance
        assert call_args[0][1] is True

    def test_with_export_vtk(self, base_settings):
        """Test with VTK export enabled."""
        from ossdbs.api import run_volume_conductor_model

        base_settings["ExportVTK"] = True

        mock_vc = MagicMock()
        mock_vc.run_full_analysis.return_value = {}
        mock_signal = MagicMock()

        run_volume_conductor_model(base_settings, mock_vc, mock_signal)

        call_args = mock_vc.run_full_analysis.call_args
        # Third positional arg is export_vtk
        assert call_args[0][2] is True

    def test_with_truncation_time(self, base_settings):
        """Test with truncation time set."""
        from ossdbs.api import run_volume_conductor_model

        mock_vc = MagicMock()
        mock_vc.run_full_analysis.return_value = {}
        mock_signal = MagicMock()

        run_volume_conductor_model(
            base_settings, mock_vc, mock_signal, truncation_time=0.005
        )

        call_kwargs = mock_vc.run_full_analysis.call_args[1]
        assert call_kwargs["truncation_time"] == 0.005

    def test_with_export_frequency(self, base_settings):
        """Test with custom export frequency."""
        from ossdbs.api import run_volume_conductor_model

        base_settings["ExportFrequency"] = 1000.0

        mock_vc = MagicMock()
        mock_vc.run_full_analysis.return_value = {}
        mock_signal = MagicMock()

        run_volume_conductor_model(base_settings, mock_vc, mock_signal)

        call_kwargs = mock_vc.run_full_analysis.call_args[1]
        assert call_kwargs["export_frequency"] == 1000.0


class TestRunStimSets:
    """Tests for run_stim_sets function."""

    @pytest.fixture
    def base_settings(self, tmp_path):
        """Base settings for StimSets tests."""
        return {
            "OutOfCore": False,
            "ExportVTK": False,
            "ActivationThresholdVTA[V-per-m]": 0.2,
            "OutputPath": str(tmp_path),
            "Mesh": {
                "AdaptiveMeshRefinement": {"Active": False},
                "MaterialRefinementSteps": 0,
            },
            "PointModel": {
                "Pathway": {"Active": False},
                "Lattice": {"Active": False},
                "VoxelLattice": {"Active": False},
            },
            "FEMOrder": 2,
        }

    def test_no_ground_contact_raises_error(self, base_settings):
        """Test that missing ground contact raises ValueError."""
        from ossdbs.api import run_stim_sets

        # Mock geometry with no ground contact
        mock_contact = MagicMock()
        mock_contact.current = 1.0  # Not -1, so not ground
        mock_contact.active = True

        mock_geometry = MagicMock()
        mock_geometry.contacts = [mock_contact]

        mock_conductivity = MagicMock()
        mock_solver = MagicMock()
        mock_signal = MagicMock()
        mock_signal.current_controlled = True

        with pytest.raises(ValueError, match="No ground contact"):
            run_stim_sets(
                base_settings,
                mock_geometry,
                mock_conductivity,
                mock_solver,
                mock_signal,
            )

    @patch("ossdbs.api.prepare_volume_conductor_model")
    def test_with_ground_contact(self, mock_prepare_vcm, base_settings):
        """Test StimSets with proper ground contact."""
        from ossdbs.api import run_stim_sets

        # Mock ground contact
        ground_contact = MagicMock()
        ground_contact.current = -1.0
        ground_contact.active = True
        ground_contact.name = "Contact_Ground"

        # Mock active contact
        active_contact = MagicMock()
        active_contact.current = 1.0
        active_contact.active = True
        active_contact.name = "Contact_1"

        mock_geometry = MagicMock()
        mock_geometry.contacts = [ground_contact, active_contact]
        mock_geometry.get_contact_index.return_value = 0

        mock_conductivity = MagicMock()
        mock_solver = MagicMock()
        mock_signal = MagicMock()
        mock_signal.current_controlled = True

        # Mock volume conductor
        mock_vc = MagicMock()
        mock_vc.run_full_analysis.return_value = {}
        mock_prepare_vcm.return_value = mock_vc

        run_stim_sets(
            base_settings,
            mock_geometry,
            mock_conductivity,
            mock_solver,
            mock_signal,
        )

        # Should have run analysis for the active contact (not ground)
        mock_vc.run_full_analysis.assert_called()


class TestRunPAM:
    """Tests for run_PAM function."""

    def test_pam_not_available(self):
        """Test that run_PAM raises error when NEURON not installed."""
        from ossdbs.api import run_PAM

        if not PAM_AVAILABLE:
            with pytest.raises(RuntimeError, match="PAM not available"):
                run_PAM({})

    @pytest.mark.skipif(not PAM_AVAILABLE)
    @patch("ossdbs.axon_processing.get_neuron_model")
    def test_pam_basic_run(self, mock_get_neuron_model, tmp_path):
        """Test basic PAM run with mocked components."""
        from ossdbs.api import run_PAM

        # Create pathway file
        pathway_file = tmp_path / "pathway.json"
        pathway_file.write_text('{"Axon_Model_Type": "MRG2002"}')

        # Create time domain solution file
        output_path = tmp_path / "output"
        output_path.mkdir()
        solution_file = output_path / "oss_time_result_PAM.h5"
        solution_file.touch()

        settings = {
            "PathwayFile": str(pathway_file),
            "OutputPath": str(output_path),
            "StimSets": {"Active": False},
            "Scaling": 1.0,
            "ScalingIndex": 0,
        }

        # Mock neuron model
        mock_model = MagicMock()
        mock_model.load_solution.return_value = MagicMock()
        mock_get_neuron_model.return_value = mock_model

        run_PAM(settings)

        mock_get_neuron_model.assert_called_once()
        mock_model.load_solution.assert_called_once()
        mock_model.process_pathways.assert_called_once()

    @pytest.mark.skipif(not PAM_AVAILABLE)
    @patch("ossdbs.axon_processing.get_neuron_model")
    def test_pam_with_current_vector(self, mock_get_neuron_model, tmp_path):
        """Test PAM with CurrentVector (StimSets mode)."""
        from ossdbs.api import run_PAM

        # Create pathway file
        pathway_file = tmp_path / "pathway.json"
        pathway_file.write_text('{"Axon_Model_Type": "MRG2002"}')

        # Create output directory
        output_path = tmp_path / "output"
        output_path.mkdir()

        # Create unit solution directories
        for i in range(3):
            contact_dir = tmp_path / f"outputE1C{i + 1}"
            contact_dir.mkdir()
            (contact_dir / "oss_time_result_PAM.h5").touch()

        settings = {
            "PathwayFile": str(pathway_file),
            "OutputPath": str(output_path),
            "StimSets": {"Active": True, "StimSetsFile": None},
            "CurrentVector": [1.0, 0.5, 0.0],
            "Scaling": 1.0,
        }

        # Mock neuron model
        mock_model = MagicMock()
        mock_model.load_unit_solutions.return_value = [MagicMock()]
        mock_model.superimpose_unit_solutions.return_value = MagicMock()
        mock_get_neuron_model.return_value = mock_model

        run_PAM(settings)

        mock_model.load_unit_solutions.assert_called_once()
        mock_model.superimpose_unit_solutions.assert_called_once()
        mock_model.process_pathways.assert_called_once()

    @pytest.mark.skipif(not PAM_AVAILABLE)
    @patch("ossdbs.axon_processing.get_neuron_model")
    def test_pam_stimsets_no_current_vector_raises_error(
        self, mock_get_neuron_model, tmp_path
    ):
        """Test PAM raises error when StimSets active but no CurrentVector."""
        from ossdbs.api import run_PAM

        # Create pathway file
        pathway_file = tmp_path / "pathway.json"
        pathway_file.write_text('{"Axon_Model_Type": "MRG2002"}')

        output_path = tmp_path / "output"
        output_path.mkdir()

        settings = {
            "PathwayFile": str(pathway_file),
            "OutputPath": str(output_path),
            "StimSets": {"Active": True, "StimSetsFile": None},
            # No CurrentVector provided
        }

        mock_get_neuron_model.return_value = MagicMock()

        with pytest.raises(ValueError, match="StimSetsFile or a CurrentVector"):
            run_PAM(settings)
