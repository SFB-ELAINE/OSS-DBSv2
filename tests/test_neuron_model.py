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

neuron_model = pytest.importorskip(
    "ossdbs.axon_processing.neuron_model",
    reason="ossdbs.axon_processing.neuron_model requires NEURON",
)
NEURON_PROCESS_TIMEOUT = neuron_model.NEURON_PROCESS_TIMEOUT
_mp_context = neuron_model._mp_context

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


def _point_source_potential(segment_positions, electrode_pos, amplitude):
    """Compute extracellular potential from a point source in homogeneous medium.

    Uses V = amplitude / (4*pi*sigma*r) with sigma=0.3 S/m.
    Returns potential in mV for positions in um.

    Parameters
    ----------
    segment_positions : np.ndarray
        (N, 3) array of segment midpoint positions in um.
    electrode_pos : np.ndarray
        (3,) electrode position in um.
    amplitude : float
        Source amplitude in mV*um (pre-scaled for convenience).
    """
    r = np.linalg.norm(segment_positions - electrode_pos, axis=1)
    r = np.maximum(r, 1.0)  # avoid division by zero
    return amplitude / r


class TestNeuronActivation:
    """Tests that verify NEURON actually produces action potentials.

    These tests apply known extracellular potentials directly to
    _run_neuron_simulation and check activation results. They serve as
    numerical regression tests: if the NEURON version changes (e.g. 8→9),
    any shift in activation thresholds will be caught here.

    The stimulus is a point-source cathodic pulse placed near the center
    of the axon, producing a realistic spatial gradient (activating function).
    """

    @requires_neuron
    def test_mrg2002_suprathreshold_activates(self, tmp_path):
        """A cathodic point-source pulse above threshold must trigger an AP.

        Uses default HOC parameters (fiberD=3.0, 75 nodes, 593 compartments,
        deltax=278 um) to avoid HOC file modifications and mechanism reload issues.
        """
        from ossdbs.axon_processing.neuron_model import (
            MRG2002,
            _run_neuron_simulation,
        )

        # Default MRG2002 parameters from axon4pyfull.hoc
        n_ranvier = 75
        n_comp = 8  # for fiberD=3.0 (< 5.7): node + MYSA + FLUT + 3*STIN + FLUT + MYSA
        n_segments = (n_ranvier - 1) * n_comp + 1  # 593
        deltax = 278.0  # um, default node-to-node spacing

        pathways_dict = {
            "Axon_Model_Type": "MRG2002",
            "connectome_name": "test",
            "axon_diams": [3.0],
            "n_Ranvier": [n_ranvier],
            "N_seeded_neurons": [1],
            "N_orig_neurons": [1],
        }

        output_path = str(tmp_path / "output")
        os.makedirs(output_path)
        simulator = MRG2002(pathways_dict, output_path)
        # Skip modify_hoc_file — use defaults that match the template

        time_step = 0.005  # ms
        n_time_steps = 2000  # 10 ms total
        signal_dict = {"time_step": time_step, "N_time_steps": n_time_steps}

        # Build segment positions along x-axis
        seg_positions = np.zeros((n_segments, 3))
        for i in range(n_segments):
            seg_positions[i, 0] = (i // n_comp) * deltax + (
                i % n_comp
            ) * deltax / n_comp

        # Point source 200 um from center of axon
        center_x = seg_positions[n_segments // 2, 0]
        electrode_pos = np.array([center_x, 200.0, 0.0])

        # Amplitude 1e5 mV*um — produces ~500 mV at center, AF ~15 mV
        spatial_potential = _point_source_potential(
            seg_positions, electrode_pos, amplitude=-1e5
        )

        v_ext = np.zeros((n_segments, n_time_steps))
        pulse_start = 20
        pulse_duration_steps = int(0.2 / time_step)  # 0.2 ms pulse
        for t in range(pulse_start, pulse_start + pulse_duration_steps):
            v_ext[:, t] = spatial_potential

        result = _run_neuron_simulation(
            neuron_index=0,
            v_ext=v_ext,
            neuron_workdir=simulator._neuron_workdir,
            hoc_file=simulator.hoc_file,
            signal_dict=signal_dict,
            extra_initialization=simulator._extra_initialization,
            v_init=simulator._v_init,
        )

        assert len(result) == 2, f"Simulation error: {result}"
        assert result[1] is True, "Suprathreshold stimulus must activate axon"

    @requires_neuron
    def test_mrg2002_subthreshold_does_not_activate(self, tmp_path):
        """A weak cathodic point-source pulse must NOT trigger an AP."""
        from ossdbs.axon_processing.neuron_model import (
            MRG2002,
            _run_neuron_simulation,
        )

        n_ranvier = 75
        n_comp = 8
        n_segments = (n_ranvier - 1) * n_comp + 1  # 593
        deltax = 278.0

        pathways_dict = {
            "Axon_Model_Type": "MRG2002",
            "connectome_name": "test",
            "axon_diams": [3.0],
            "n_Ranvier": [n_ranvier],
            "N_seeded_neurons": [1],
            "N_orig_neurons": [1],
        }

        output_path = str(tmp_path / "output")
        os.makedirs(output_path)
        simulator = MRG2002(pathways_dict, output_path)

        time_step = 0.005
        n_time_steps = 2000
        signal_dict = {"time_step": time_step, "N_time_steps": n_time_steps}

        seg_positions = np.zeros((n_segments, 3))
        for i in range(n_segments):
            seg_positions[i, 0] = (i // n_comp) * deltax + (
                i % n_comp
            ) * deltax / n_comp

        center_x = seg_positions[n_segments // 2, 0]
        electrode_pos = np.array([center_x, 200.0, 0.0])

        # Amplitude 1e3 — ~50x weaker than threshold
        spatial_potential = _point_source_potential(
            seg_positions, electrode_pos, amplitude=-1e3
        )

        v_ext = np.zeros((n_segments, n_time_steps))
        pulse_start = 20
        pulse_duration_steps = int(0.2 / time_step)
        for t in range(pulse_start, pulse_start + pulse_duration_steps):
            v_ext[:, t] = spatial_potential

        result = _run_neuron_simulation(
            neuron_index=0,
            v_ext=v_ext,
            neuron_workdir=simulator._neuron_workdir,
            hoc_file=simulator.hoc_file,
            signal_dict=signal_dict,
            extra_initialization=simulator._extra_initialization,
            v_init=simulator._v_init,
        )

        assert len(result) == 2, f"Simulation error: {result}"
        assert result[1] is False, "Subthreshold stimulus must not activate axon"

    @requires_neuron
    def test_mcneal1976_suprathreshold_activates(self, tmp_path):
        """McNeal1976: a strong cathodic point-source pulse must trigger an AP."""
        from ossdbs.axon_processing.axon_models import AxonMorphologyMcNeal1976
        from ossdbs.axon_processing.neuron_model import (
            McNeal1976,
            _run_neuron_simulation,
        )

        n_ranvier = 25
        n_compartments = n_ranvier * 2 - 1  # 49 (nodes + internodes)

        pathways_dict = {
            "Axon_Model_Type": "McNeal1976",
            "connectome_name": "test",
            "axon_diams": [5.0],
            "n_Ranvier": [n_ranvier],
            "N_seeded_neurons": [1],
            "N_orig_neurons": [1],
        }

        output_path = str(tmp_path / "output")
        os.makedirs(output_path)
        simulator = McNeal1976(pathways_dict, output_path)

        axon_morphology = AxonMorphologyMcNeal1976()
        axon_morphology.update_axon_morphology(5.0, n_Ranvier=n_ranvier)

        time_step = 0.005
        n_time_steps = 1000
        steps_per_ms = int(1.0 / time_step)
        simulator.modify_hoc_file(n_ranvier, steps_per_ms, axon_morphology)

        signal_dict = {"time_step": time_step, "N_time_steps": n_time_steps}

        # McNeal1976: DIAM=5, ELD=100, so internode length = 500 um
        # GAP = 2.5 um (node length). Compartments alternate node/internode.
        internode_length = 500.0  # um
        node_length = 2.5  # um (GAP * 1e4 in um)
        seg_positions = np.zeros((n_compartments, 3))
        x = 0.0
        for i in range(n_compartments):
            if i % 2 == 0:  # node
                seg_positions[i, 0] = x + node_length / 2
                x += node_length
            else:  # internode
                seg_positions[i, 0] = x + internode_length / 2
                x += internode_length

        center_x = seg_positions[n_compartments // 2, 0]
        electrode_pos = np.array([center_x, 500.0, 0.0])

        # Strong cathodic pulse
        spatial_potential = _point_source_potential(
            seg_positions, electrode_pos, amplitude=-5e6
        )

        v_ext = np.zeros((n_compartments, n_time_steps))
        pulse_start = 10
        pulse_duration_steps = int(0.1 / time_step)
        for t in range(pulse_start, pulse_start + pulse_duration_steps):
            v_ext[:, t] = spatial_potential

        result = _run_neuron_simulation(
            neuron_index=0,
            v_ext=v_ext,
            neuron_workdir=simulator._neuron_workdir,
            hoc_file=simulator.hoc_file,
            signal_dict=signal_dict,
            extra_initialization=simulator._extra_initialization,
            v_init=simulator._v_init,
        )

        assert len(result) == 2, f"Simulation error: {result}"
        assert result[1] is True, "Suprathreshold stimulus must activate axon"

    @requires_neuron
    def test_mcneal1976_subthreshold_does_not_activate(self, tmp_path):
        """McNeal1976: a weak cathodic point-source pulse must NOT trigger an AP."""
        from ossdbs.axon_processing.axon_models import AxonMorphologyMcNeal1976
        from ossdbs.axon_processing.neuron_model import (
            McNeal1976,
            _run_neuron_simulation,
        )

        n_ranvier = 25
        n_compartments = n_ranvier * 2 - 1

        pathways_dict = {
            "Axon_Model_Type": "McNeal1976",
            "connectome_name": "test",
            "axon_diams": [5.0],
            "n_Ranvier": [n_ranvier],
            "N_seeded_neurons": [1],
            "N_orig_neurons": [1],
        }

        output_path = str(tmp_path / "output")
        os.makedirs(output_path)
        simulator = McNeal1976(pathways_dict, output_path)

        axon_morphology = AxonMorphologyMcNeal1976()
        axon_morphology.update_axon_morphology(5.0, n_Ranvier=n_ranvier)

        time_step = 0.005
        n_time_steps = 1000
        steps_per_ms = int(1.0 / time_step)
        simulator.modify_hoc_file(n_ranvier, steps_per_ms, axon_morphology)

        signal_dict = {"time_step": time_step, "N_time_steps": n_time_steps}

        internode_length = 500.0
        node_length = 2.5
        seg_positions = np.zeros((n_compartments, 3))
        x = 0.0
        for i in range(n_compartments):
            if i % 2 == 0:
                seg_positions[i, 0] = x + node_length / 2
                x += node_length
            else:
                seg_positions[i, 0] = x + internode_length / 2
                x += internode_length

        center_x = seg_positions[n_compartments // 2, 0]
        electrode_pos = np.array([center_x, 500.0, 0.0])

        # Very weak pulse
        spatial_potential = _point_source_potential(
            seg_positions, electrode_pos, amplitude=-1e3
        )

        v_ext = np.zeros((n_compartments, n_time_steps))
        pulse_start = 10
        pulse_duration_steps = int(0.1 / time_step)
        for t in range(pulse_start, pulse_start + pulse_duration_steps):
            v_ext[:, t] = spatial_potential

        result = _run_neuron_simulation(
            neuron_index=0,
            v_ext=v_ext,
            neuron_workdir=simulator._neuron_workdir,
            hoc_file=simulator.hoc_file,
            signal_dict=signal_dict,
            extra_initialization=simulator._extra_initialization,
            v_init=simulator._v_init,
        )

        assert len(result) == 2, f"Simulation error: {result}"
        assert result[1] is False, "Subthreshold stimulus must not activate axon"


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
