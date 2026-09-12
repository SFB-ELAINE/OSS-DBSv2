# Copyright 2026 OSS-DBS contributors
# SPDX-License-Identifier: GPL-3.0-or-later
"""EasyVVUQ runner for a small OSS-DBS impedance uncertainty example.

The runner is intentionally independent of EasyVVUQ. EasyVVUQ only creates
``input.json`` files with three uncertain parameters; this script converts those
parameters into a valid OSS-DBS input dictionary, runs OSS-DBS, and writes the
``output.csv`` file decoded by EasyVVUQ.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from ossdbs.dielectric_model.colecole4 import BloodColeCole4Default
from ossdbs.main import main_run

REPO_ROOT = Path(__file__).resolve().parents[3]
BASE_INPUT = REPO_ROOT / "input_files" / "inputTest.json"
DEFAULT_FREQUENCIES = [
    3000.0,
    6207.41424334437,
    12843.997196158187,
    700716.4407270363,
    3000000.0,
]


def _as_float(
    parameters: dict[str, Any], name: str, default: float | None = None
) -> float:
    """Extract a parameter value as float."""
    value = parameters.get(name, default)
    if value is None:
        raise KeyError(f"Missing required EasyVVUQ parameter: {name}")
    return float(value)


def _cole_cole_parameters(parameters: dict[str, Any]) -> dict[str, Any]:
    """Return Cole-Cole parameters for the encapsulation layer.

    EasyVVUQ varies only three entries. The remaining Blood defaults are kept
    fixed so the example is easy to inspect.
    """
    alpha = np.array(BloodColeCole4Default.alpha, dtype=float)
    eps_delta = np.array(BloodColeCole4Default.eps_delta, dtype=float)
    tau = np.array(BloodColeCole4Default.tau, dtype=float)

    # The tutorial uses one-based names following the usual Cole-Cole notation.
    alpha[2] = _as_float(parameters, "alpha_3")
    eps_delta[2] = _as_float(
        parameters,
        "eps_delta_3",
        parameters.get("delta_epsilon_3"),
    )

    return {
        "eps_inf": float(BloodColeCole4Default.eps_inf),
        "sigma": _as_float(parameters, "sigma"),
        "alpha": alpha.tolist(),
        "eps_delta": eps_delta.tolist(),
        "tau": tau.tolist(),
    }


def build_ossdbs_input(parameters: dict[str, Any], run_dir: Path) -> dict[str, Any]:
    """Build an OSS-DBS input dictionary from EasyVVUQ parameters."""
    settings = json.loads(BASE_INPUT.read_text(encoding="utf-8"))
    input_files = REPO_ROOT / "input_files"
    output_dir = run_dir / "ossdbs_output"
    output_dir.mkdir(parents=True, exist_ok=True)

    settings["StimulationFolder"] = str(run_dir)
    settings["OutputPath"] = str(output_dir)
    settings["ComputeImpedance"] = True
    settings["ExportVTK"] = False
    settings["ExportElectrode"] = False
    settings["EQSMode"] = True
    settings["DielectricAccuracy"] = 0.0

    settings["StimulationSignal"] = {
        "CurrentControlled": True,
        "Type": "Multisine",
        "ListOfFrequencies": DEFAULT_FREQUENCIES,
    }

    material_distribution = settings["MaterialDistribution"]
    material_distribution["MRIPath"] = str(
        input_files / material_distribution["MRIPath"]
    )
    if material_distribution.get("DiffusionTensorActive") and material_distribution.get(
        "DTIPath"
    ):
        material_distribution["DTIPath"] = str(
            input_files / material_distribution["DTIPath"]
        )

    encapsulation = settings["Electrodes"][0]["EncapsulationLayer"]
    encapsulation["Thickness[mm]"] = _as_float(
        parameters, "encapsulation_thickness", 0.2
    )
    encapsulation["Material"] = "Blood"
    encapsulation["DielectricModel"] = "ColeCole4"
    encapsulation["DielectricParameters"] = _cole_cole_parameters(parameters)
    encapsulation.setdefault("MaxMeshSize", 0.1)

    return settings


def write_easyvvuq_output(ossdbs_output_dir: Path, target_csv: Path) -> None:
    """Convert OSS-DBS impedance output to the CSV format expected by EasyVVUQ."""
    impedance_file = ossdbs_output_dir / "impedance.csv"
    if not impedance_file.exists():
        raise FileNotFoundError(f"OSS-DBS did not create {impedance_file}")

    with impedance_file.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    frequency_key = "frequency" if "frequency" in rows[0] else "freq"
    with target_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["freq", "real", "imag"])
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "freq": row[frequency_key],
                    "real": row["real"],
                    "imag": row["imag"],
                }
            )


def main() -> None:
    """Run one EasyVVUQ-generated OSS-DBS simulation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_json", type=Path, help="EasyVVUQ encoded input JSON")
    args = parser.parse_args()

    input_json = args.input_json.resolve()
    run_dir = input_json.parent
    parameters = json.loads(input_json.read_text(encoding="utf-8"))

    settings = build_ossdbs_input(parameters, run_dir)
    (run_dir / "ossdbs_input.json").write_text(
        json.dumps(settings, indent=2), encoding="utf-8"
    )

    main_run(settings)
    write_easyvvuq_output(Path(settings["OutputPath"]), run_dir / "output.csv")


if __name__ == "__main__":
    main()
