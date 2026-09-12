"""Tests for the standalone ``ImpedanceAnalyzer``.

The analyzer computes the NxN admittance matrix Y and impedance matrix
Z = Y^-1 via the superposition approach, decoupled from
``VolumeConductor.run_full_analysis``. The admittance-matrix code lived
on ``VolumeConductor`` before P2 and was auto-triggered on
``n_active + n_floating > 2``; it now lives in its own module and has
no coupling to the stimulation pipeline.
"""

import json
import os
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

import ossdbs
from ossdbs.fem.analysis import ImpedanceAnalyzer
from ossdbs.fem.volume_conductor.floating import VolumeConductorFloating
from ossdbs.fem.volume_conductor.floating_impedance import (
    VolumeConductorFloatingImpedance,
)
from ossdbs.fem.volume_conductor.nonfloating import VolumeConductorNonFloating
from ossdbs.utils.settings import Settings


@pytest.fixture
def settings_fixture():
    json_path = os.path.join(
        os.getcwd(), "input_test_cases/input_case1/input_homogeneous.json"
    )
    with open(json_path) as file:
        settings = json.load(file)
    return Settings(settings).complete_settings()


def _build_volume_conductor(settings: dict):
    """Build a VolumeConductor from a settings dict (case-1 homogeneous)."""
    local_settings = deepcopy(settings)
    local_settings["MaterialDistribution"]["MRIPath"] = os.path.join(
        os.getcwd(), "input_files/sub-John_Doe/JD_segmask.nii.gz"
    )
    local_settings["MaterialDistribution"]["DiffusionTensorActive"] = True
    local_settings["MaterialDistribution"]["DTIPath"] = os.path.join(
        os.getcwd(), "input_files/sub-John_Doe/JD_DTI_NormMapping.nii.gz"
    )
    mri_image, _ = ossdbs.load_images(local_settings)
    brain_region = ossdbs.create_bounding_box(local_settings["BrainRegion"])
    electrodes = ossdbs.generate_electrodes(local_settings)
    brain = ossdbs.BrainGeometry(local_settings["BrainRegion"]["Shape"], brain_region)
    geometry = ossdbs.ModelGeometry(brain, electrodes)
    ossdbs.set_contact_and_encapsulation_layer_properties(local_settings, geometry)

    dielectric_model = ossdbs.prepare_dielectric_properties(local_settings)
    materials = local_settings["MaterialDistribution"]["MRIMapping"]
    conductivity = ossdbs.ConductivityCF(
        mri_image,
        brain_region,
        dielectric_model,
        materials,
        geometry.encapsulation_layers,
        complex_data=local_settings["EQSMode"],
    )

    solver = ossdbs.prepare_solver(local_settings)

    floating_mode = geometry.get_floating_mode()
    vc_classes = {
        "Floating": VolumeConductorFloating,
        "FloatingImpedance": VolumeConductorFloatingImpedance,
        "NonFloating": VolumeConductorNonFloating,
    }
    VolumeConductorClass = vc_classes.get(floating_mode, VolumeConductorNonFloating)
    return VolumeConductorClass(
        geometry,
        conductivity,
        solver,
        local_settings["FEMOrder"],
        local_settings["Mesh"],
    )


class TestImpedanceAnalyzer:
    def test_compute_populates_matrices_with_expected_shape(self, settings_fixture):
        """``compute`` fills admittance/impedance with shape (n_freq, N, N)
        and Y @ Z approximately the identity.
        """
        vc = _build_volume_conductor(settings_fixture)
        analyzer = ImpedanceAnalyzer(vc)
        frequencies = [1000.0, 10000.0]
        analyzer.compute(frequencies)

        N = len(analyzer.contact_names)
        assert N == 2, f"case-1 fixture has 2 active contacts; got {N}"

        assert analyzer.admittance_matrices.shape == (len(frequencies), N, N)
        assert analyzer.impedance_matrices.shape == (len(frequencies), N, N)

        for freq_idx in range(len(frequencies)):
            Y = analyzer.admittance_matrices[freq_idx]
            Z = analyzer.impedance_matrices[freq_idx]
            product = Y @ Z
            assert np.allclose(product, np.eye(N), atol=1e-6), (
                f"Z must invert Y at frequency index {freq_idx}"
            )

    def test_results_unavailable_before_compute(self, settings_fixture):
        vc = _build_volume_conductor(settings_fixture)
        analyzer = ImpedanceAnalyzer(vc)
        with pytest.raises(RuntimeError, match="compute"):
            _ = analyzer.admittance_matrices
        with pytest.raises(RuntimeError, match="compute"):
            _ = analyzer.impedance_matrices
        with pytest.raises(RuntimeError, match="compute"):
            _ = analyzer.frequencies

    def test_export_writes_flat_csv_schema(self, settings_fixture, tmp_path):
        """``export`` writes two CSVs with columns
        ``freq,row,col,real,imag``.
        """
        vc = _build_volume_conductor(settings_fixture)
        analyzer = ImpedanceAnalyzer(vc)
        frequencies = [1000.0]
        analyzer.compute(frequencies)
        analyzer.export(str(tmp_path))

        expected_cols = ["freq", "row", "col", "real", "imag"]
        N = len(analyzer.contact_names)

        for label in ("admittance_matrix", "impedance_matrix"):
            df = pd.read_csv(os.path.join(tmp_path, f"{label}.csv"))
            assert list(df.columns) == expected_cols, (
                f"{label}.csv column order mismatch"
            )
            assert len(df) == len(frequencies) * N * N, (
                f"{label}.csv row count mismatch"
            )

    def test_export_before_compute_raises(self, settings_fixture, tmp_path):
        vc = _build_volume_conductor(settings_fixture)
        analyzer = ImpedanceAnalyzer(vc)
        with pytest.raises(RuntimeError, match="compute"):
            analyzer.export(str(tmp_path))

    def test_analyzer_does_not_mutate_stimulation_state(self, settings_fixture):
        """Running the analyzer must not populate ``VolumeConductor``
        stimulation state (``_impedances``, ``_stimulation_variable``,
        ``_scale_factor``). These belong to ``run_full_analysis``.
        """
        vc = _build_volume_conductor(settings_fixture)
        analyzer = ImpedanceAnalyzer(vc)
        analyzer.compute([1000.0])

        assert vc._impedances is None
        assert vc._stimulation_variable is None
        assert vc._free_stimulation_variable is None
