# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import sys
from unittest.mock import MagicMock, patch

import pytest


class TestMainCLI:
    """Tests for main() CLI entry point."""

    def test_cli_argument_parsing(self, tmp_path):
        """Test CLI argument parsing."""
        # Create a minimal input file
        input_file = tmp_path / "input.json"
        input_file.write_text('{"test": "data"}')

        with patch("ossdbs.main.main_run") as mock_main_run:
            with patch.object(sys, "argv", ["ossdbs", str(input_file)]):
                from ossdbs.main import main

                main()

                # Verify main_run was called with the loaded settings
                mock_main_run.assert_called_once()
                call_args = mock_main_run.call_args[0][0]
                assert call_args["test"] == "data"
                assert "StimulationFolder" in call_args

    def test_cli_with_loglevel(self, tmp_path):
        """Test CLI with custom log level."""
        input_file = tmp_path / "input.json"
        input_file.write_text('{"test": "data"}')

        with patch("ossdbs.main.main_run") as mock_main_run:
            with patch("ossdbs.main.set_logger") as mock_set_logger:
                with patch.object(
                    sys, "argv", ["ossdbs", "--loglevel", "10", str(input_file)]
                ):
                    from ossdbs.main import main

                    main()

                    mock_set_logger.assert_called_once_with(level=10)
                    mock_main_run.assert_called_once()

    def test_cli_stimulation_folder_set(self, tmp_path):
        """Test that StimulationFolder is set correctly."""
        input_file = tmp_path / "input.json"
        input_file.write_text('{"key": "value"}')

        with patch("ossdbs.main.main_run") as mock_main_run:
            with patch.object(sys, "argv", ["ossdbs", str(input_file)]):
                from ossdbs.main import main

                main()

                call_args = mock_main_run.call_args[0][0]
                assert call_args["StimulationFolder"] == str(tmp_path)


class TestMainRun:
    """Tests for main_run() function."""

    @pytest.fixture
    def mock_settings(self, tmp_path):
        """Create mock settings for testing."""
        output_path = tmp_path / "output"
        output_path.mkdir()

        return {
            "OutputPath": str(output_path),
            "StimulationFolder": str(tmp_path),
            "FailFlag": "test_run",
            "BrainRegion": {
                "Shape": "Ellipsoid",
                "Center": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 0},
                "Dimension": {"x[mm]": 40, "y[mm]": 40, "z[mm]": 40},
            },
            "Electrodes": [
                {
                    "Name": "BostonScientificVercise",
                    "Direction": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 1},
                    "TipPosition": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 0},
                    "Rotation[Degrees]": 0,
                    "Contacts": [
                        {"Contact_ID": 1, "Active": True, "Voltage[V]": 1.0},
                        {"Contact_ID": 2, "Active": True, "Voltage[V]": 0.0},
                    ],
                }
            ],
            "MaterialDistribution": {
                "MRIPath": "test.nii.gz",
                "MRIMapping": {"CSF": 1, "Gray matter": 2},
                "DiffusionTensorActive": False,
                "WMMasking": False,
            },
            "DielectricModel": {"Type": "ColeCole4", "CustomParameters": None},
            "Mesh": {
                "LoadMesh": False,
                "SaveMesh": False,
                "MeshingHypothesis": {"Type": "Default"},
            },
            "Solver": {
                "Type": "CG",
                "Preconditioner": "bddc",
                "PreconditionerKwargs": {},
                "MaximumSteps": 100,
                "Precision": 1e-8,
            },
            "StimulationSignal": {
                "Type": "Multisine",
                "CurrentControlled": False,
                "ListOfFrequencies": [1000.0],
            },
            "FEMOrder": 2,
            "EQSMode": False,
            "StimSets": {"Active": False},
            "PathwayFile": None,
            "PointModel": {
                "Pathway": {"Active": False},
                "Lattice": {"Active": False},
                "VoxelLattice": {"Active": False},
            },
            "ActivationThresholdVTA[V-per-m]": 0.2,
            "OutOfCore": False,
            "ExportVTK": False,
            "ExportFrequency": None,
            "ComputeImpedance": False,
        }

    @patch("ossdbs.main.ngsolve.TaskManager")
    @patch("ossdbs.main.run_volume_conductor_model")
    @patch("ossdbs.main.prepare_stimulation_signal")
    @patch("ossdbs.main.prepare_volume_conductor_model")
    @patch("ossdbs.main.prepare_solver")
    @patch("ossdbs.main.ConductivityCF")
    @patch("ossdbs.main.prepare_dielectric_properties")
    @patch("ossdbs.main.set_contact_and_encapsulation_layer_properties")
    @patch("ossdbs.main.ModelGeometry")
    @patch("ossdbs.main.build_brain_model")
    @patch("ossdbs.main.generate_electrodes")
    @patch("ossdbs.main.load_images")
    @patch("ossdbs.main.Settings")
    @patch("ossdbs.main.TypeChecker")
    @patch("ossdbs.main.log_to_file")
    def test_main_run_basic_flow(
        self,
        mock_log_to_file,
        mock_type_checker,
        mock_settings_class,
        mock_load_images,
        mock_generate_electrodes,
        mock_build_brain_model,
        mock_model_geometry,
        mock_set_contact,
        mock_prepare_dielectric,
        mock_conductivity_cf,
        mock_prepare_solver,
        mock_prepare_vcm,
        mock_prepare_signal,
        mock_run_vcm,
        mock_task_manager,
        mock_settings,
        tmp_path,
    ):
        """Test basic main_run flow with all components mocked."""
        from ossdbs.main import main_run

        # Configure mocks
        mock_settings_class.return_value.complete_settings.return_value = mock_settings
        mock_load_images.return_value = (MagicMock(), None)
        mock_generate_electrodes.return_value = [MagicMock()]
        mock_build_brain_model.return_value = MagicMock()
        mock_model_geometry.return_value = MagicMock()
        mock_prepare_dielectric.return_value = {}
        mock_conductivity_cf.return_value = MagicMock()
        mock_prepare_solver.return_value = MagicMock()
        mock_prepare_vcm.return_value = MagicMock()
        mock_prepare_signal.return_value = MagicMock()
        mock_run_vcm.return_value = {"solve": 1.0}
        mock_task_manager.return_value.__enter__ = MagicMock()
        mock_task_manager.return_value.__exit__ = MagicMock()

        # Run
        main_run(mock_settings)

        # Verify key functions were called
        mock_load_images.assert_called_once()
        mock_generate_electrodes.assert_called_once()
        mock_build_brain_model.assert_called_once()
        mock_model_geometry.assert_called_once()
        mock_prepare_dielectric.assert_called_once()
        mock_prepare_solver.assert_called_once()
        mock_prepare_vcm.assert_called_once()
        mock_prepare_signal.assert_called_once()
        mock_run_vcm.assert_called_once()

        # Verify success file created and fail file removed
        success_file = tmp_path / f"success_{mock_settings['FailFlag']}.txt"
        assert success_file.exists()

    @patch("ossdbs.main.ngsolve.TaskManager")
    @patch("ossdbs.main.run_stim_sets")
    @patch("ossdbs.main.prepare_stimulation_signal")
    @patch("ossdbs.main.prepare_volume_conductor_model")
    @patch("ossdbs.main.prepare_solver")
    @patch("ossdbs.main.ConductivityCF")
    @patch("ossdbs.main.prepare_dielectric_properties")
    @patch("ossdbs.main.set_contact_and_encapsulation_layer_properties")
    @patch("ossdbs.main.ModelGeometry")
    @patch("ossdbs.main.build_brain_model")
    @patch("ossdbs.main.generate_electrodes")
    @patch("ossdbs.main.load_images")
    @patch("ossdbs.main.Settings")
    @patch("ossdbs.main.TypeChecker")
    @patch("ossdbs.main.log_to_file")
    def test_main_run_with_stimsets(
        self,
        mock_log_to_file,
        mock_type_checker,
        mock_settings_class,
        mock_load_images,
        mock_generate_electrodes,
        mock_build_brain_model,
        mock_model_geometry,
        mock_set_contact,
        mock_prepare_dielectric,
        mock_conductivity_cf,
        mock_prepare_solver,
        mock_prepare_vcm,
        mock_prepare_signal,
        mock_run_stim_sets,
        mock_task_manager,
        mock_settings,
        tmp_path,
    ):
        """Test main_run with StimSets active."""
        from ossdbs.main import main_run

        # Enable StimSets
        mock_settings["StimSets"]["Active"] = True

        # Configure mocks
        mock_settings_class.return_value.complete_settings.return_value = mock_settings
        mock_load_images.return_value = (MagicMock(), None)
        mock_generate_electrodes.return_value = [MagicMock()]
        mock_build_brain_model.return_value = MagicMock()
        mock_model_geometry.return_value = MagicMock()
        mock_prepare_dielectric.return_value = {}
        mock_conductivity_cf.return_value = MagicMock()
        mock_prepare_solver.return_value = MagicMock()
        mock_prepare_vcm.return_value = MagicMock()
        mock_prepare_signal.return_value = MagicMock()
        mock_task_manager.return_value.__enter__ = MagicMock()
        mock_task_manager.return_value.__exit__ = MagicMock()

        # Run
        main_run(mock_settings)

        # Verify run_stim_sets was called instead of run_volume_conductor_model
        mock_run_stim_sets.assert_called_once()

    @patch("ossdbs.main.ngsolve.TaskManager")
    @patch("ossdbs.main.run_volume_conductor_model")
    @patch("ossdbs.main.prepare_stimulation_signal")
    @patch("ossdbs.main.prepare_volume_conductor_model")
    @patch("ossdbs.main.prepare_solver")
    @patch("ossdbs.main.ConductivityCF")
    @patch("ossdbs.main.prepare_dielectric_properties")
    @patch("ossdbs.main.set_contact_and_encapsulation_layer_properties")
    @patch("ossdbs.main.ModelGeometry")
    @patch("ossdbs.main.build_brain_model")
    @patch("ossdbs.main.generate_electrodes")
    @patch("ossdbs.main.load_images")
    @patch("ossdbs.main.Settings")
    @patch("ossdbs.main.TypeChecker")
    @patch("ossdbs.main.log_to_file")
    @patch("ossdbs.main.generate_signal")
    def test_main_run_with_truncation(
        self,
        mock_generate_signal,
        mock_log_to_file,
        mock_type_checker,
        mock_settings_class,
        mock_load_images,
        mock_generate_electrodes,
        mock_build_brain_model,
        mock_model_geometry,
        mock_set_contact,
        mock_prepare_dielectric,
        mock_conductivity_cf,
        mock_prepare_solver,
        mock_prepare_vcm,
        mock_prepare_signal,
        mock_run_vcm,
        mock_task_manager,
        mock_settings,
        tmp_path,
    ):
        """Test main_run with truncation ratio set."""
        from ossdbs.main import main_run

        # Set truncation ratio
        mock_settings["TruncateAfterActivePartRatio"] = 20.0

        # Configure mocks
        mock_settings_class.return_value.complete_settings.return_value = mock_settings
        mock_load_images.return_value = (MagicMock(), None)
        mock_generate_electrodes.return_value = [MagicMock()]
        mock_build_brain_model.return_value = MagicMock()
        mock_model_geometry.return_value = MagicMock()
        mock_prepare_dielectric.return_value = {}
        mock_conductivity_cf.return_value = MagicMock()
        mock_prepare_solver.return_value = MagicMock()
        mock_prepare_vcm.return_value = MagicMock()
        mock_prepare_signal.return_value = MagicMock()
        mock_run_vcm.return_value = {"solve": 1.0}
        mock_task_manager.return_value.__enter__ = MagicMock()
        mock_task_manager.return_value.__exit__ = MagicMock()

        # Mock signal with get_active_time
        mock_signal = MagicMock()
        mock_signal.get_active_time.return_value = 0.001
        mock_generate_signal.return_value = mock_signal

        # Run
        main_run(mock_settings)

        # Verify generate_signal was called for truncation calculation
        mock_generate_signal.assert_called_once()
        mock_signal.get_active_time.assert_called_once()

        # Verify run_volume_conductor_model was called with truncation_time
        call_kwargs = mock_run_vcm.call_args[1]
        assert call_kwargs["truncation_time"] == 20.0 * 0.001

    @patch("ossdbs.main.ConductivityCF")
    @patch("ossdbs.main.prepare_dielectric_properties")
    @patch("ossdbs.main.set_contact_and_encapsulation_layer_properties")
    @patch("ossdbs.main.ModelGeometry")
    @patch("ossdbs.main.build_brain_model")
    @patch("ossdbs.main.generate_electrodes")
    @patch("ossdbs.main.load_images")
    @patch("ossdbs.main.Settings")
    @patch("ossdbs.main.TypeChecker")
    @patch("ossdbs.main.log_to_file")
    def test_main_run_invalid_truncation_ratio(
        self,
        mock_log_to_file,
        mock_type_checker,
        mock_settings_class,
        mock_load_images,
        mock_generate_electrodes,
        mock_build_brain_model,
        mock_model_geometry,
        mock_set_contact,
        mock_prepare_dielectric,
        mock_conductivity_cf,
        mock_settings,
        tmp_path,
    ):
        """Test main_run raises error for invalid truncation ratio."""
        from ossdbs.main import main_run

        # Set invalid truncation ratio (not a float)
        mock_settings["TruncateAfterActivePartRatio"] = "invalid"

        mock_settings_class.return_value.complete_settings.return_value = mock_settings
        mock_load_images.return_value = (MagicMock(), None)
        mock_generate_electrodes.return_value = [MagicMock()]
        mock_build_brain_model.return_value = MagicMock()
        mock_model_geometry.return_value = MagicMock()
        mock_prepare_dielectric.return_value = {}
        mock_conductivity_cf.return_value = MagicMock()

        with pytest.raises(ValueError, match="floating-point"):
            main_run(mock_settings)

    @patch("ossdbs.main.ngsolve.TaskManager")
    @patch("ossdbs.main.run_volume_conductor_model")
    @patch("ossdbs.main.prepare_stimulation_signal")
    @patch("ossdbs.main.prepare_volume_conductor_model")
    @patch("ossdbs.main.prepare_solver")
    @patch("ossdbs.main.ConductivityCF")
    @patch("ossdbs.main.prepare_dielectric_properties")
    @patch("ossdbs.main.set_contact_and_encapsulation_layer_properties")
    @patch("ossdbs.main.ModelGeometry")
    @patch("ossdbs.main.build_brain_model")
    @patch("ossdbs.main.generate_electrodes")
    @patch("ossdbs.main.load_images")
    @patch("ossdbs.main.Settings")
    @patch("ossdbs.main.TypeChecker")
    @patch("ossdbs.main.log_to_file")
    def test_main_run_geometry_retry(
        self,
        mock_log_to_file,
        mock_type_checker,
        mock_settings_class,
        mock_load_images,
        mock_generate_electrodes,
        mock_build_brain_model,
        mock_model_geometry,
        mock_set_contact,
        mock_prepare_dielectric,
        mock_conductivity_cf,
        mock_prepare_solver,
        mock_prepare_vcm,
        mock_prepare_signal,
        mock_run_vcm,
        mock_task_manager,
        mock_settings,
        tmp_path,
    ):
        """Test main_run retries geometry with rotation on failure."""
        from ossdbs.main import main_run

        # Configure mocks
        mock_settings_class.return_value.complete_settings.return_value = mock_settings
        mock_load_images.return_value = (MagicMock(), None)
        mock_generate_electrodes.return_value = [MagicMock()]
        mock_build_brain_model.return_value = MagicMock()

        # First call to ModelGeometry raises RuntimeError, second succeeds
        mock_model_geometry.side_effect = [RuntimeError("Geometry failed"), MagicMock()]

        mock_prepare_dielectric.return_value = {}
        mock_conductivity_cf.return_value = MagicMock()
        mock_prepare_solver.return_value = MagicMock()
        mock_prepare_vcm.return_value = MagicMock()
        mock_prepare_signal.return_value = MagicMock()
        mock_run_vcm.return_value = {"solve": 1.0}
        mock_task_manager.return_value.__enter__ = MagicMock()
        mock_task_manager.return_value.__exit__ = MagicMock()

        # Run
        main_run(mock_settings)

        # Verify build_brain_model was called twice (once with rotation)
        assert mock_build_brain_model.call_count == 2
        # Second call should have rotate_initial_geo=True
        second_call_kwargs = mock_build_brain_model.call_args_list[1][1]
        assert second_call_kwargs.get("rotate_initial_geo") is True

    @patch("ossdbs.main.ngsolve.TaskManager")
    @patch("ossdbs.main.run_volume_conductor_model")
    @patch("ossdbs.main.prepare_stimulation_signal")
    @patch("ossdbs.main.prepare_volume_conductor_model")
    @patch("ossdbs.main.prepare_solver")
    @patch("ossdbs.main.ConductivityCF")
    @patch("ossdbs.main.prepare_dielectric_properties")
    @patch("ossdbs.main.set_contact_and_encapsulation_layer_properties")
    @patch("ossdbs.main.ModelGeometry")
    @patch("ossdbs.main.build_brain_model")
    @patch("ossdbs.main.generate_electrodes")
    @patch("ossdbs.main.load_images")
    @patch("ossdbs.main.Settings")
    @patch("ossdbs.main.TypeChecker")
    @patch("ossdbs.main.log_to_file")
    def test_main_run_with_pathway_file(
        self,
        mock_log_to_file,
        mock_type_checker,
        mock_settings_class,
        mock_load_images,
        mock_generate_electrodes,
        mock_build_brain_model,
        mock_model_geometry,
        mock_set_contact,
        mock_prepare_dielectric,
        mock_conductivity_cf,
        mock_prepare_solver,
        mock_prepare_vcm,
        mock_prepare_signal,
        mock_run_vcm,
        mock_task_manager,
        mock_settings,
        tmp_path,
    ):
        """Test main_run with PathwayFile set (logs message)."""
        from ossdbs.main import main_run

        # Set pathway file
        mock_settings["PathwayFile"] = "pathway.h5"

        # Configure mocks
        mock_settings_class.return_value.complete_settings.return_value = mock_settings
        mock_load_images.return_value = (MagicMock(), None)
        mock_generate_electrodes.return_value = [MagicMock()]
        mock_build_brain_model.return_value = MagicMock()
        mock_model_geometry.return_value = MagicMock()
        mock_prepare_dielectric.return_value = {}
        mock_conductivity_cf.return_value = MagicMock()
        mock_prepare_solver.return_value = MagicMock()
        mock_prepare_vcm.return_value = MagicMock()
        mock_prepare_signal.return_value = MagicMock()
        mock_run_vcm.return_value = {"solve": 1.0}
        mock_task_manager.return_value.__enter__ = MagicMock()
        mock_task_manager.return_value.__exit__ = MagicMock()

        # Run - should complete without error
        main_run(mock_settings)

        # Verify run completed (success file exists)
        success_file = tmp_path / f"success_{mock_settings['FailFlag']}.txt"
        assert success_file.exists()
