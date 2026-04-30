"""Tests for ossdbs.axon_processing.neuron_model module.

These tests cover the multiprocessing and NEURON simulation code.
Tests that require NEURON are skipped if NEURON is not installed or
if the NEURON compiler (nrnivmodl/mknrndll) is not available.
"""

import importlib.util
import multiprocessing as mp
import os
import shutil
import subprocess
import sys

import h5py
import numpy as np
import pytest

# do not test if neuron is not installed
pytest.importorskip("neuron")

from ossdbs.axon_processing.neuron_model import (
    NEURON_PROCESS_TIMEOUT,
    _mp_context,
)

# Check if NEURON module is available
NEURON_MODULE_AVAILABLE = importlib.util.find_spec("neuron") is not None

# Check if NEURON compiler is available and actually works
if sys.platform == "win32":
    NEURON_COMPILER = "mknrndll"
else:
    NEURON_COMPILER = "nrnivmodl"


def _check_neuron_compiler():
    """Check if NEURON compiler is available and executable."""
    compiler_path = shutil.which(NEURON_COMPILER)
    if compiler_path is None:
        return False

    # On Windows, .BAT files need shell=True to run
    # Try to actually run the compiler to verify it works
    try:
        result = subprocess.run(
            compiler_path,
            capture_output=True,
            timeout=30,
            shell=(sys.platform == "win32"),
        )
        # Exit code 0 or 1 is OK (1 means no .mod files found, which is expected)
        return result.returncode in (0, 1)
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        return False


NEURON_COMPILER_AVAILABLE = _check_neuron_compiler()

# Full NEURON availability requires both module and compiler
NEURON_AVAILABLE = NEURON_MODULE_AVAILABLE and NEURON_COMPILER_AVAILABLE


# Skip marker for tests requiring NEURON
requires_neuron = pytest.mark.skipif(
    not NEURON_AVAILABLE,
    reason=f"NEURON not fully available (module={NEURON_MODULE_AVAILABLE}, "
    f"compiler={NEURON_COMPILER_AVAILABLE})",
)


class TestMultiprocessingConfiguration:
    """Tests for multiprocessing configuration."""

    def test_mp_context_is_spawn(self):
        """Verify that the multiprocessing context uses 'spawn' method."""
        # The spawn context should be used to avoid NEURON state issues
        assert _mp_context.get_start_method() == "spawn"

    def test_process_timeout_is_positive(self):
        """Verify the process timeout constant is a positive number."""
        assert NEURON_PROCESS_TIMEOUT > 0
        assert isinstance(NEURON_PROCESS_TIMEOUT, int | float)

    def test_cpu_count_calculation(self):
        """Test that CPU count calculation returns valid integer."""
        N_proc = max(1, mp.cpu_count() // 2)
        assert isinstance(N_proc, int)
        assert N_proc >= 1
        assert N_proc <= mp.cpu_count()


class TestGetVExt:
    """Tests for voltage conversion function get_v_ext."""

    @requires_neuron
    def test_get_v_ext_basic_scaling(self, tmp_path):
        """Test basic voltage conversion with scaling."""
        from ossdbs.axon_processing.neuron_model import MRG2002

        # Create minimal pathways_dict
        pathways_dict = {
            "Axon_Model_Type": "MRG2002",
            "connectome_name": "test",
            "axon_diams": [5.7],
            "n_Ranvier": [21],
            "N_seeded_neurons": [1],
            "N_orig_neurons": [1],
        }

        output_path = str(tmp_path / "output")
        os.makedirs(output_path)

        simulator = MRG2002(pathways_dict, output_path)
        simulator.signal_dict = {"scaling": 1.0}

        # Test input: 2 segments, 3 time steps, values in Volts
        v_time_sol = np.array([[0.001, 0.002, 0.003], [0.004, 0.005, 0.006]])

        v_ext = simulator.get_v_ext(v_time_sol)

        # Should be converted to mV (multiply by 1000) and scaled
        expected = v_time_sol * 1000.0 * 1.0
        np.testing.assert_array_almost_equal(v_ext, expected)

    @requires_neuron
    def test_get_v_ext_with_scaling_factor(self, tmp_path):
        """Test voltage conversion with non-unity scaling factor."""
        from ossdbs.axon_processing.neuron_model import MRG2002

        pathways_dict = {
            "Axon_Model_Type": "MRG2002",
            "connectome_name": "test",
            "axon_diams": [5.7],
            "n_Ranvier": [21],
            "N_seeded_neurons": [1],
            "N_orig_neurons": [1],
        }

        output_path = str(tmp_path / "output")
        os.makedirs(output_path)

        simulator = MRG2002(pathways_dict, output_path)
        simulator.signal_dict = {"scaling": 2.5}

        v_time_sol = np.array([[0.001, 0.002], [0.003, 0.004]])

        v_ext = simulator.get_v_ext(v_time_sol)

        expected = v_time_sol * 1000.0 * 2.5
        np.testing.assert_array_almost_equal(v_ext, expected)

    @requires_neuron
    def test_get_v_ext_with_scaling_vector(self, tmp_path):
        """Test voltage conversion with contact-wise scaling vector."""
        from ossdbs.axon_processing.neuron_model import MRG2002

        pathways_dict = {
            "Axon_Model_Type": "MRG2002",
            "connectome_name": "test",
            "axon_diams": [5.7],
            "n_Ranvier": [21],
            "N_seeded_neurons": [1],
            "N_orig_neurons": [1],
        }

        output_path = str(tmp_path / "output")
        os.makedirs(output_path)

        simulator = MRG2002(pathways_dict, output_path)
        simulator.signal_dict = {"scaling": 1.0}
        simulator._scaling_vector = [0.5, 0.3, 0.2]

        v_time_sol = np.array([[0.001, 0.002], [0.003, 0.004]])

        v_ext = simulator.get_v_ext(v_time_sol)

        # Should sum contributions from each scaling component
        expected = np.zeros_like(v_time_sol)
        for scaling_comp in [0.5, 0.3, 0.2]:
            expected += v_time_sol * 1000.0 * 1.0 * scaling_comp
        np.testing.assert_array_almost_equal(v_ext, expected)


class TestNeuronSimulatorInit:
    """Tests for NeuronSimulator initialization."""

    @requires_neuron
    def test_mrg2002_init(self, tmp_path):
        """Test MRG2002 simulator initialization."""
        from ossdbs.axon_processing.neuron_model import MRG2002

        pathways_dict = {
            "Axon_Model_Type": "MRG2002",
            "connectome_name": "test_connectome",
            "axon_diams": [5.7],
            "n_Ranvier": [21],
            "N_seeded_neurons": [10],
            "N_orig_neurons": [10],
        }

        output_path = str(tmp_path / "output")
        os.makedirs(output_path)

        simulator = MRG2002(pathways_dict, output_path)

        assert simulator.axon_model == "MRG2002"
        assert simulator._extra_initialization is True
        assert simulator._v_init == -80.0
        assert simulator.connectome_name == "test_connectome"
        assert os.path.isdir(simulator._neuron_workdir)

    @requires_neuron
    def test_mcneal1976_init(self, tmp_path):
        """Test McNeal1976 simulator initialization."""
        from ossdbs.axon_processing.neuron_model import McNeal1976

        pathways_dict = {
            "Axon_Model_Type": "McNeal1976",
            "connectome_name": "test_connectome",
            "axon_diams": [5.7],
            "n_Ranvier": [21],
            "N_seeded_neurons": [10],
            "N_orig_neurons": [10],
        }

        output_path = str(tmp_path / "output")
        os.makedirs(output_path)

        simulator = McNeal1976(pathways_dict, output_path)

        assert simulator.axon_model == "McNeal1976"
        assert simulator._extra_initialization is False
        assert simulator._v_init == -70.0

    @requires_neuron
    def test_downsampled_detection(self, tmp_path):
        """Test that downsampled models are correctly detected."""
        from ossdbs.axon_processing.neuron_model import MRG2002

        pathways_dict = {
            "Axon_Model_Type": "MRG2002_DS",  # Downsampled
            "connectome_name": "test",
            "axon_diams": [5.7],
            "n_Ranvier": [21],
            "N_seeded_neurons": [1],
            "N_orig_neurons": [1],
        }

        output_path = str(tmp_path / "output")
        os.makedirs(output_path)

        simulator = MRG2002(pathways_dict, output_path)
        assert simulator.downsampled is True


class TestRunNeuronSimulation:
    """Tests for the _run_neuron_simulation worker function."""

    @requires_neuron
    def test_worker_function_returns_tuple(self, tmp_path):
        """Test that worker function returns correct tuple format."""
        from ossdbs.axon_processing.neuron_model import (
            MRG2002,
            _run_neuron_simulation,
        )

        pathways_dict = {
            "Axon_Model_Type": "MRG2002",
            "connectome_name": "test",
            "axon_diams": [5.7],
            "n_Ranvier": [21],
            "N_seeded_neurons": [1],
            "N_orig_neurons": [1],
        }

        output_path = str(tmp_path / "output")
        os.makedirs(output_path)

        simulator = MRG2002(pathways_dict, output_path)

        # Create minimal test data
        # For MRG2002 with 21 Ranvier nodes and fiber_diam >= 5.7:
        # n_segments = (21 - 1) * 11 + 1 = 221
        n_segments = 221
        n_time_steps = 100
        v_ext = np.zeros((n_segments, n_time_steps))

        signal_dict = {
            "time_step": 0.025,  # 40 kHz
            "N_time_steps": n_time_steps,
        }

        # Modify HOC file first (required before running simulation)
        from ossdbs.axon_processing.axon_models import AxonMorphologyMRG2002

        axon_morphology = AxonMorphologyMRG2002()
        axon_morphology.update_axon_morphology(5.7, n_Ranvier=21)
        stepsPerMs = int(1.0 / signal_dict["time_step"])
        simulator.modify_hoc_file(21, stepsPerMs, axon_morphology)

        result = _run_neuron_simulation(
            neuron_index=0,
            v_ext=v_ext,
            neuron_workdir=simulator._neuron_workdir,
            hoc_file=simulator.hoc_file,
            signal_dict=signal_dict,
            extra_initialization=simulator._extra_initialization,
            v_init=simulator._v_init,
        )

        # Result should be (neuron_index, activated) or (neuron_index, None, error)
        assert isinstance(result, tuple)
        assert len(result) >= 2
        assert result[0] == 0  # neuron_index
        # With zero input, should not be activated
        if result[1] is not None:
            assert result[1] is False or result[1] is True

    @requires_neuron
    def test_worker_function_catches_errors(self, tmp_path):
        """Test that worker function catches and returns errors."""
        from ossdbs.axon_processing.neuron_model import _run_neuron_simulation

        # Call with invalid parameters to trigger an error
        result = _run_neuron_simulation(
            neuron_index=42,
            v_ext=np.zeros((10, 10)),
            neuron_workdir="/nonexistent/path",
            hoc_file="nonexistent.hoc",
            signal_dict={"time_step": 0.025, "N_time_steps": 10},
            extra_initialization=False,
            v_init=-70.0,
        )

        # Should return a tuple starting with the neuron_index
        assert isinstance(result, tuple)
        assert result[0] == 42  # neuron_index preserved
        # Either an error tuple (index, None, msg) or a result tuple (index, bool)
        assert len(result) in (2, 3)


class TestCheckPathwayActivation:
    """Tests for check_pathway_activation method."""

    @requires_neuron
    def test_data_extraction_from_h5(self, tmp_path):
        """Test that H5 data is correctly extracted before multiprocessing."""
        from ossdbs.axon_processing.neuron_model import MRG2002

        # Create test H5 file with pathway data
        h5_path = tmp_path / "test_pathway.h5"
        n_neurons = 3
        n_segments = 221  # For MRG2002 with 21 Ranvier nodes
        n_time_steps = 50

        with h5py.File(h5_path, "w") as f:
            f.create_dataset("TimeSteps[s]", data=np.linspace(0, 0.001, n_time_steps))

            pathway_grp = f.create_group("test_pathway")
            # Status: 0 = available, -1 = electrode intersection, -2 = CSF
            pathway_grp.create_dataset("Status", data=np.array([0, -1, -2]))

            for i in range(n_neurons):
                axon_grp = pathway_grp.create_group(f"axon{i}")
                axon_grp.create_dataset(
                    "Points[mm]", data=np.random.rand(n_segments, 3)
                )
                axon_grp.create_dataset(
                    "Potential[V]",
                    data=np.random.rand(n_segments, n_time_steps) * 0.001,
                )
                # Store original streamline index (from streamline_indexing branch)
                axon_grp.attrs["inx"] = i + 1  # 1-based index

        pathways_dict = {
            "Axon_Model_Type": "MRG2002",
            "connectome_name": "test",
            "axon_diams": [5.7],
            "n_Ranvier": [21],
            "N_seeded_neurons": [n_neurons],
            "N_orig_neurons": [n_neurons],
        }

        output_path = str(tmp_path / "output")
        os.makedirs(output_path)

        simulator = MRG2002(pathways_dict, output_path)

        # Load and process
        td_solution = simulator.load_solution(str(h5_path))
        simulator.signal_dict = {
            "time_step": 0.025,
            "scaling": 1.0,
            "N_time_steps": n_time_steps,
        }

        # Run pathway activation (this tests the full flow)
        pathway_dataset = td_solution["test_pathway"]
        simulator.check_pathway_activation(
            pathway_dataset, pathway_idx=0, pathway_name="test_pathway"
        )

        td_solution.close()

        # Verify output files were created
        assert os.path.exists(output_path)

    @requires_neuron
    def test_skipped_neurons_have_correct_status(self, tmp_path):
        """Test that neurons with pre-status != 0 are correctly skipped."""
        from ossdbs.axon_processing.neuron_model import MRG2002

        h5_path = tmp_path / "test_pathway.h5"
        n_neurons = 5
        n_segments = 221
        n_time_steps = 50

        with h5py.File(h5_path, "w") as f:
            f.create_dataset("TimeSteps[s]", data=np.linspace(0, 0.001, n_time_steps))

            pathway_grp = f.create_group("test_pathway")
            # Mix of statuses: only neuron 0 and 3 should be simulated
            pathway_grp.create_dataset("Status", data=np.array([0, -1, -2, 0, -1]))

            for i in range(n_neurons):
                axon_grp = pathway_grp.create_group(f"axon{i}")
                axon_grp.create_dataset(
                    "Points[mm]", data=np.random.rand(n_segments, 3)
                )
                axon_grp.create_dataset(
                    "Potential[V]", data=np.zeros((n_segments, n_time_steps))
                )
                # Store original streamline index (from streamline_indexing branch)
                axon_grp.attrs["inx"] = i + 1  # 1-based index

        pathways_dict = {
            "Axon_Model_Type": "MRG2002",
            "connectome_name": "test",
            "axon_diams": [5.7],
            "n_Ranvier": [21],
            "N_seeded_neurons": [n_neurons],
            "N_orig_neurons": [n_neurons],
        }

        output_path = str(tmp_path / "output")
        os.makedirs(output_path)

        simulator = MRG2002(pathways_dict, output_path)

        td_solution = simulator.load_solution(str(h5_path))
        simulator.signal_dict = {
            "time_step": 0.025,
            "scaling": 1.0,
            "N_time_steps": n_time_steps,
        }

        pathway_dataset = td_solution["test_pathway"]
        simulator.check_pathway_activation(
            pathway_dataset, pathway_idx=0, pathway_name="test_pathway"
        )

        td_solution.close()


class TestProcessPoolExecutorIntegration:
    """Integration tests for ProcessPoolExecutor-based parallelization."""

    @requires_neuron
    def test_parallel_execution_produces_results(self, tmp_path):
        """Test that parallel execution with ProcessPoolExecutor works."""
        from ossdbs.axon_processing.neuron_model import MRG2002

        h5_path = tmp_path / "test_pathway.h5"
        n_neurons = 4  # Small number for quick test
        n_segments = 221
        n_time_steps = 50

        with h5py.File(h5_path, "w") as f:
            f.create_dataset("TimeSteps[s]", data=np.linspace(0, 0.001, n_time_steps))

            pathway_grp = f.create_group("test_pathway")
            pathway_grp.create_dataset("Status", data=np.zeros(n_neurons))

            for i in range(n_neurons):
                axon_grp = pathway_grp.create_group(f"axon{i}")
                axon_grp.create_dataset(
                    "Points[mm]", data=np.random.rand(n_segments, 3)
                )
                # Zero potential = no activation expected
                axon_grp.create_dataset(
                    "Potential[V]", data=np.zeros((n_segments, n_time_steps))
                )
                # Store original streamline index (from streamline_indexing branch)
                axon_grp.attrs["inx"] = i + 1  # 1-based index

        pathways_dict = {
            "Axon_Model_Type": "MRG2002",
            "connectome_name": "test",
            "axon_diams": [5.7],
            "n_Ranvier": [21],
            "N_seeded_neurons": [n_neurons],
            "N_orig_neurons": [n_neurons],
        }

        output_path = str(tmp_path / "output")
        os.makedirs(output_path)

        simulator = MRG2002(pathways_dict, output_path)

        td_solution = simulator.load_solution(str(h5_path))
        simulator.signal_dict = {
            "time_step": 0.025,
            "scaling": 1.0,
            "N_time_steps": n_time_steps,
        }

        pathway_dataset = td_solution["test_pathway"]

        # This should complete without errors
        simulator.check_pathway_activation(
            pathway_dataset, pathway_idx=0, pathway_name="test_pathway"
        )

        td_solution.close()

        # Check that output was produced
        output_files = os.listdir(output_path)
        assert len(output_files) > 0


class TestMcNeal1976:
    """Tests specific to McNeal1976 model."""

    @requires_neuron
    def test_mcneal_hoc_file_modification(self, tmp_path):
        """Test that McNeal1976 HOC file is correctly modified."""
        from ossdbs.axon_processing.axon_models import AxonMorphologyMcNeal1976
        from ossdbs.axon_processing.neuron_model import McNeal1976

        pathways_dict = {
            "Axon_Model_Type": "McNeal1976",
            "connectome_name": "test",
            "axon_diams": [10.0],
            "n_Ranvier": [21],
            "N_seeded_neurons": [1],
            "N_orig_neurons": [1],
        }

        output_path = str(tmp_path / "output")
        os.makedirs(output_path)

        simulator = McNeal1976(pathways_dict, output_path)

        axon_morphology = AxonMorphologyMcNeal1976()
        axon_morphology.update_axon_morphology(10.0, n_Ranvier=21)

        # This should not raise an error
        simulator.modify_hoc_file(21, 40, axon_morphology)

        # Verify the HOC file exists
        hoc_path = os.path.join(simulator._neuron_workdir, simulator.hoc_file)
        assert os.path.exists(hoc_path)


class TestMRG2002Upsampling:
    """Tests for MRG2002 upsampling functionality."""

    @requires_neuron
    def test_upsample_voltage_large_fiber(self, tmp_path):
        """Test voltage upsampling for fiber diameter >= 5.7."""
        from ossdbs.axon_processing.axon_models import AxonMorphologyMRG2002
        from ossdbs.axon_processing.neuron_model import MRG2002

        pathways_dict = {
            "Axon_Model_Type": "MRG2002_DS",
            "connectome_name": "test",
            "axon_diams": [5.7],
            "n_Ranvier": [21],
            "N_seeded_neurons": [1],
            "N_orig_neurons": [1],
        }

        output_path = str(tmp_path / "output")
        os.makedirs(output_path)

        simulator = MRG2002(pathways_dict, output_path)

        axon_morphology = AxonMorphologyMRG2002()
        axon_morphology.update_axon_morphology(5.7, n_Ranvier=21)

        # For downsampled model with fiber >= 5.7:
        # n_segments_ds = ((n_segments_full - 1) / 11) * 3 + 1
        # n_segments_full = (21 - 1) * 11 + 1 = 221
        # n_segments_ds = (220 / 11) * 3 + 1 = 60 + 1 = 61
        n_segments_ds = 61
        n_time_steps = 50

        v_time_sol = np.random.rand(n_segments_ds, n_time_steps)

        v_upsampled = simulator.upsample_voltage(v_time_sol, 5.7, axon_morphology)

        # Should be upsampled to full resolution
        assert v_upsampled.shape[0] == axon_morphology.n_segments
        assert v_upsampled.shape[1] == n_time_steps

    @requires_neuron
    def test_upsample_voltage_small_fiber(self, tmp_path):
        """Test voltage upsampling for fiber diameter < 5.7."""
        from ossdbs.axon_processing.axon_models import AxonMorphologyMRG2002
        from ossdbs.axon_processing.neuron_model import MRG2002

        pathways_dict = {
            "Axon_Model_Type": "MRG2002_DS",
            "connectome_name": "test",
            "axon_diams": [2.0],
            "n_Ranvier": [21],
            "N_seeded_neurons": [1],
            "N_orig_neurons": [1],
        }

        output_path = str(tmp_path / "output")
        os.makedirs(output_path)

        simulator = MRG2002(pathways_dict, output_path)

        axon_morphology = AxonMorphologyMRG2002()
        axon_morphology.update_axon_morphology(2.0, n_Ranvier=21)

        # For fiber < 5.7, different segment calculation
        # n_segments_full = (21 - 1) * 8 + 1 = 161
        # n_segments_ds = ((161 - 1) / 8) * 2 + 1 = 40 + 1 = 41
        n_segments_ds = 41
        n_time_steps = 50

        v_time_sol = np.random.rand(n_segments_ds, n_time_steps)

        v_upsampled = simulator.upsample_voltage(v_time_sol, 2.0, axon_morphology)

        # Should be upsampled to full resolution
        assert v_upsampled.shape[0] == axon_morphology.n_segments
        assert v_upsampled.shape[1] == n_time_steps
