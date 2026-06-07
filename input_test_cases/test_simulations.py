# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

"""Pytest-based simulation tests for OSS-DBS.

These tests run full simulations and compare outputs against expected results.
They are marked with @pytest.mark.simulation and excluded from regular pytest
runs via addopts in pyproject.toml.

Usage:
    # Run all simulation tests
    pytest input_test_cases/test_simulations.py -m simulation -v

    # Run only fast tests
    pytest input_test_cases/test_simulations.py -m "simulation and not slow" -v

    # Run specific category
    pytest input_test_cases/test_simulations.py -k "brain_material" -v
"""

import json
import logging
import os
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import pytest

from ossdbs.api import run_PAM
from ossdbs.main import main_run

_logger = logging.getLogger(__name__)

BASE = Path(__file__).parent

# Tolerances
IMPEDANCE_ATOL = 1e-2
VTA_DICE_ATOL = 1e-2
VTA_VOLUME_RTOL = 1e-2
FLOATING_POTENTIAL_ATOL = 1e-2

# ---------------------------------------------------------------------------
# Test case definitions
# ---------------------------------------------------------------------------

# Each entry: (test_id, input_dir, input_json, marks, checks)
# checks is a dict mapping check type -> {output_file, desired_file}
# For simple impedance-only checks, just list "impedance".

TEST_CASES = [
    # Case 1: Brain material
    {
        "id": "brain_material_homogeneous",
        "input_dir": "input_case1",
        "input_json": "input_case1/input_homogeneous.json",
        "marks": [],
        "checks": ["impedance"],
    },
    {
        "id": "brain_material_inhomogeneous",
        "input_dir": "input_case1",
        "input_json": "input_case1/input_inhomogeneous.json",
        "marks": [],
        "checks": ["impedance"],
    },
    # Case 2: Custom parameters
    {
        "id": "custom_electrode",
        "input_dir": "input_case2",
        "input_json": "input_case2/input_custom_electrode.json",
        "marks": [],
        "checks": ["impedance"],
    },
    {
        "id": "custom_material",
        "input_dir": "input_case2",
        "input_json": "input_case2/input_custom_material.json",
        "marks": [],
        "checks": ["impedance"],
    },
    # Case 3: Case grounding
    {
        "id": "case_grounding_qs",
        "input_dir": "input_case3",
        "input_json": "input_case3/input_case_grounding.json",
        "marks": [],
        "checks": ["impedance"],
    },
    {
        "id": "case_grounding_eqs",
        "input_dir": "input_case3",
        "input_json": "input_case3/input_case_grounding_EQS.json",
        "marks": [],
        "checks": ["impedance"],
    },
    # Case 4: Current controlled
    {
        "id": "current_controlled",
        "input_dir": "input_case4",
        "input_json": "input_case4/input_current_controlled.json",
        "marks": [],
        "checks": ["impedance"],
    },
    {
        "id": "multi_current",
        "input_dir": "input_case4",
        "input_json": "input_case4/input_multi_current.json",
        "marks": [],
        "checks": ["floating_potentials"],
    },
    # Case 5: Stimulation signals
    {
        "id": "stimulation_signal",
        "input_dir": "input_case5",
        "input_json": "input_case5/input_stimulation_signal.json",
        "marks": ["slow"],
        "checks": ["impedance"],
    },
    {
        "id": "stimulation_signal_counterpulse",
        "input_dir": "input_case5",
        "input_json": "input_case5/input_stimulation_signal_counterpulse.json",
        "marks": ["slow"],
        "checks": ["impedance"],
    },
    # Case 6: Floating contacts
    {
        "id": "floating_contacts",
        "input_dir": "input_case6",
        "input_json": "input_case6/input_floating.json",
        "marks": [],
        "checks": ["impedance"],
    },
    # Case 7: VTA
    {
        "id": "vta_standard",
        "input_dir": "input_case7",
        "input_json": "input_case7/input_vta.json",
        "marks": ["slow", "vta"],
        "checks": ["impedance", "vta_nifti", "vta_volume"],
    },
    {
        "id": "vta_out_of_core",
        "input_dir": "input_case7",
        "input_json": "input_case7/input_vta_out_of_core.json",
        "marks": ["slow", "vta"],
        "checks": ["impedance", "vta_nifti", "vta_volume"],
    },
    # Case 8: Pathway Activation Modeling
    {
        "id": "pam_standard",
        "input_dir": "input_case8",
        "input_json": "input_case8/input_pathway.json",
        "marks": ["slow", "requires_neuron", "pam"],
        "checks": ["impedance", "pathway_activation"],
    },
    {
        "id": "pam_out_of_core",
        "input_dir": "input_case8",
        "input_json": "input_case8/input_pathway_out_of_core.json",
        "marks": ["slow", "requires_neuron", "pam"],
        "checks": ["impedance", "pathway_activation"],
    },
    # Case 9: PAM with StimSets (produces per-contact results)
    {
        "id": "pam_stimsets",
        "input_dir": "input_case9",
        "input_json": "input_case9/input_pathway.json",
        "marks": ["slow", "requires_neuron", "pam"],
        "checks": ["stimsets_impedance", "stimsets_pathway_activation"],
    },
    # Case 10: Floating with surface impedance
    {
        "id": "floating_impedance",
        "input_dir": "input_case10",
        "input_json": "input_case10/input_floating_impedance.json",
        "marks": ["slow", "surface_impedance"],
        "checks": ["impedance"],
    },
    {
        "id": "floating_impedance_noground",
        "input_dir": "input_case10",
        "input_json": "input_case10/input_floating_impedance_noground.json",
        "marks": ["slow", "surface_impedance"],
        "checks": ["floating_potentials", "currents"],
    },
    # Case 11: Surface impedance (Lempka2009 models)
    {
        "id": "surface_impedance_homogeneous",
        "input_dir": "input_case11",
        "input_json": (
            "input_case11/input_homogeneous_surface_impedance_lempka2009.json"
        ),
        "marks": ["slow", "surface_impedance"],
        "checks": ["impedance"],
    },
    {
        "id": "surface_impedance_lempka",
        "input_dir": "input_case11",
        "input_json": "input_case11/input_surface_impedance_lempka2009.json",
        "marks": ["slow", "surface_impedance"],
        "checks": ["impedance"],
    },
]

# Map mark names to pytest mark objects
MARK_MAP = {
    "slow": pytest.mark.slow,
    "requires_neuron": pytest.mark.requires_neuron,
    "vta": pytest.mark.vta,
    "pam": pytest.mark.pam,
    "floating": pytest.mark.floating,
    "surface_impedance": pytest.mark.surface_impedance,
}


def _make_params():
    """Build pytest.param objects from TEST_CASES."""
    params = []
    for tc in TEST_CASES:
        marks = [MARK_MAP[m] for m in tc["marks"] if m in MARK_MAP]
        marks.insert(0, pytest.mark.simulation)
        params.append(pytest.param(tc, id=tc["id"], marks=marks))
    return params


# ---------------------------------------------------------------------------
# Path patching (same logic as old run_all_tests.py)
# ---------------------------------------------------------------------------


def _prepare_settings(input_dir: str, input_json: str) -> dict:
    """Load JSON and patch paths to be absolute."""
    json_path = BASE / input_json
    with open(json_path) as f:
        settings = json.load(f)

    settings["StimulationFolder"] = str(BASE)

    # MRI path: resolve ../../input_files relative to the JSON's directory
    mri_path = settings["MaterialDistribution"]["MRIPath"]
    settings["MaterialDistribution"]["MRIPath"] = str(
        (json_path.parent / mri_path).resolve()
    )

    # DTI path (if set)
    dti_path = settings.get("MaterialDistribution", {}).get("DTIPath", "")
    if dti_path:
        settings["MaterialDistribution"]["DTIPath"] = str(
            (json_path.parent / dti_path).resolve()
        )

    # OutputPath: make absolute under input_dir
    settings["OutputPath"] = str(BASE / input_dir / settings["OutputPath"])

    # Pathway FileName
    pathway_file = settings.get("PointModel", {}).get("Pathway", {}).get("FileName", "")
    if pathway_file:
        settings["PointModel"]["Pathway"]["FileName"] = str(
            (json_path.parent / pathway_file).resolve()
        )

    # StimSets file
    if settings.get("StimSets", {}).get("StimSetsFile"):
        settings["StimSets"]["StimSetsFile"] = str(
            BASE / input_dir / settings["StimSets"]["StimSetsFile"]
        )

    # PathwayFile for run_PAM (pathway_parameters.json next to data.h5)
    pathway_fn = settings.get("PointModel", {}).get("Pathway", {}).get("FileName", "")
    if pathway_fn and settings.get("PointModel", {}).get("Pathway", {}).get("Active"):
        pathway_params = str(
            Path(settings["PointModel"]["Pathway"]["FileName"]).parent
            / "pathway_parameters.json"
        )
        settings["PathwayFile"] = pathway_params

    return settings


# ---------------------------------------------------------------------------
# Comparison helpers
# ---------------------------------------------------------------------------


def _compare_csv(output_csv: str, desired_csv: str, atol: float = IMPEDANCE_ATOL):
    """Compare two CSV files using np.allclose."""
    actual = pd.read_csv(output_csv).to_numpy(dtype=float)
    desired = pd.read_csv(desired_csv).to_numpy(dtype=float)
    if not np.allclose(actual, desired, atol=atol):
        abs_err = float(np.max(np.abs(actual - desired)))
        nonzero = np.abs(desired) > np.finfo(float).eps
        rel_err = (
            float(
                np.max(np.abs((actual[nonzero] - desired[nonzero]) / desired[nonzero]))
            )
            if np.any(nonzero)
            else 0.0
        )
        pytest.fail(
            f"CSV mismatch: max abs error={abs_err:.6g}, "
            f"max rel error={rel_err * 100:.4g}%\n"
            f"  output:  {output_csv}\n"
            f"  desired: {desired_csv}"
        )


def _read_nifti_vta_points(filename: str) -> set:
    """Read NIfTI and return set of VTA voxel coordinates."""
    image = nib.load(filename)
    affine = image.affine
    data = image.get_fdata()
    points = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            for k in range(data.shape[2]):
                if data[i, j, k] > 0.1:
                    x = affine[0][3] + i * affine[0][0]
                    y = affine[1][3] + j * affine[1][1]
                    z = affine[2][3] + k * affine[2][2]
                    points.append((x, y, z))
    return set(points)


def _compare_vta_nifti(output_nii: str, desired_nii: str):
    """Compare VTA NIfTI files via Dice coefficient."""
    actual = _read_nifti_vta_points(output_nii)
    desired = _read_nifti_vta_points(desired_nii)
    intersection = len(actual.intersection(desired))
    dice = (2 * intersection) / (len(actual) + len(desired))
    if not np.isclose(dice, 1.0, atol=VTA_DICE_ATOL):
        pytest.fail(
            f"VTA Dice coefficient {dice:.4f} != 1.0 (atol={VTA_DICE_ATOL})\n"
            f"  output:  {output_nii}\n"
            f"  desired: {desired_nii}"
        )


def _compare_vta_volume(output_dir: str, desired_dir: str):
    """Compare VTA volume from Lattice.json."""
    output_json = os.path.join(output_dir, "Lattice.json")
    desired_json = os.path.join(desired_dir, "Lattice.json")
    with open(output_json) as f:
        actual_vol = json.load(f)["volume"]
    with open(desired_json) as f:
        desired_vol = json.load(f)["volume"]
    if not np.isclose(actual_vol, desired_vol, rtol=VTA_VOLUME_RTOL):
        pytest.fail(
            f"VTA volume mismatch: {actual_vol:.6f} vs {desired_vol:.6f} "
            f"(rtol={VTA_VOLUME_RTOL})"
        )


def _run_pam_step(settings: dict):
    """Run the PAM (pathway activation modeling) step after FEM."""
    pam_settings = dict(settings)
    pam_settings["Scaling"] = 1.0
    pam_settings["ScalingIndex"] = None
    pam_settings["CurrentVector"] = None
    pam_settings.setdefault("StimSets", {"Active": False, "StimSetsFile": None})
    run_PAM(pam_settings)


def _compare_pathway_activation(output_dir: str, desired_dir: str):
    """Compare percent_activated from Pathway_status*.json."""
    output_files = sorted(
        f
        for f in os.listdir(output_dir)
        if f.startswith("Pathway_status") and f.endswith(".json")
    )
    assert output_files, f"No Pathway_status*.json found in {output_dir}"
    for fname in output_files:
        output_json = os.path.join(output_dir, fname)
        desired_json = os.path.join(desired_dir, fname)
        assert os.path.exists(desired_json), f"Missing baseline: {desired_json}"
        with open(output_json) as f:
            actual = json.load(f)
        with open(desired_json) as f:
            desired = json.load(f)
        if actual["percent_activated"] != desired["percent_activated"]:
            pytest.fail(
                f"Pathway activation mismatch in {fname}: "
                f"{actual['percent_activated']}% vs {desired['percent_activated']}%"
            )


# ---------------------------------------------------------------------------
# Main test
# ---------------------------------------------------------------------------


# ruff: noqa: C901
@pytest.mark.parametrize("test_case", _make_params())
def test_simulation(test_case):
    """Run simulation and validate outputs against desired results."""
    input_dir = test_case["input_dir"]
    input_json = test_case["input_json"]
    checks = test_case["checks"]

    _logger.info("Running: %s", test_case["id"])

    # Run simulation
    settings = _prepare_settings(input_dir, input_json)
    output_dir = settings["OutputPath"]
    main_run(settings)

    # Run PAM step if any check requires it
    pam_checks = {"pathway_activation", "stimsets_pathway_activation"}
    if pam_checks.intersection(checks):
        _run_pam_step(settings)

    # Desired output directory
    desired_dir = str(BASE / "desired_output" / input_dir / Path(output_dir).name)

    # Run checks
    for check in checks:
        if check == "impedance":
            _compare_csv(
                os.path.join(output_dir, "impedance.csv"),
                os.path.join(desired_dir, "impedance.csv"),
            )
        elif check == "stimsets_impedance":
            # StimSets produces per-contact result dirs (e.g. Results_PAME1C1)
            input_dir_abs = str(BASE / input_dir)
            desired_base = str(BASE / "desired_output" / input_dir)
            stimset_dirs = sorted(
                d
                for d in os.listdir(input_dir_abs)
                if d.startswith("Results_")
                and d != Path(output_dir).name
                and os.path.isfile(os.path.join(input_dir_abs, d, "impedance.csv"))
            )
            assert stimset_dirs, f"No StimSets result dirs found in {input_dir_abs}"
            for sd in stimset_dirs:
                _compare_csv(
                    os.path.join(input_dir_abs, sd, "impedance.csv"),
                    os.path.join(desired_base, sd, "impedance.csv"),
                )
        elif check == "vta_nifti":
            # Try .nii.gz first, then .nii
            for ext in [".nii.gz", ".nii"]:
                output_nii = os.path.join(output_dir, f"VTA_solution_Lattice{ext}")
                desired_nii = os.path.join(desired_dir, f"VTA_solution_Lattice{ext}")
                if os.path.exists(output_nii) and os.path.exists(desired_nii):
                    _compare_vta_nifti(output_nii, desired_nii)
                    break
            else:
                pytest.fail(f"VTA NIfTI not found in {output_dir}")
        elif check == "vta_volume":
            _compare_vta_volume(output_dir, desired_dir)
        elif check == "floating_potentials":
            _compare_csv(
                os.path.join(output_dir, "floating_potentials.csv"),
                os.path.join(desired_dir, "floating_potentials.csv"),
                atol=FLOATING_POTENTIAL_ATOL,
            )
        elif check == "impedance_matrix":
            _compare_csv(
                os.path.join(output_dir, "impedance_matrix.csv"),
                os.path.join(desired_dir, "impedance_matrix.csv"),
            )
        elif check == "currents":
            _compare_csv(
                os.path.join(output_dir, "currents.csv"),
                os.path.join(desired_dir, "currents.csv"),
            )
        elif check == "pathway_activation":
            _compare_pathway_activation(output_dir, desired_dir)
        elif check == "stimsets_pathway_activation":
            # StimSets writes Pathway_status files to the main OutputPath
            _compare_pathway_activation(output_dir, desired_dir)
