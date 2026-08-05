"""Validate the StimSets VTA scale against a full direct simulation.

For a single monopolar current (contact `CONTACT_1B` at `CURRENT_MA`), this:
  1. runs a full direct simulation (non-StimSets) with that current, the
     Lattice point model, and the activation threshold set, so OSS-DBS
     also computes the NGSolve mesh-integrated VTA volume (VCM report);
  2. counts the VTA on the lattice from the direct simulation's
     `E_field_Lattice.csv` (complex modulus > threshold);
  3. counts the VTA from the existing StimSets unit solution
     (`Results_VTA_defaultE1C<contact>/E_field_Lattice.csv`) scaled by the
     same current, i.e. the StimSets superposition for this monopolar case.

If superposition + scale are correct, (2) and (3) match, and both agree
with the mesh-integrated volume (1) up to lattice-discretization error.
"""

import json
import logging
import os

import pandas as pd

import ossdbs
from ossdbs.main import main_run

CONFIG_FILE = "oss-dbs_parameters_vta.json"
CONTACT_1B = 1  # 1-based contact index to stimulate
CURRENT_MA = -2.0  # stimulation current in mA (monopolar)
THRESHOLD_VPM = 200.0
OUT_DIR = "Results_VTA_validation_direct"
UNIT_DIR = f"Results_VTA_defaultE1C{CONTACT_1B}"


def magnitude_vta_volume(csv_path, scale, spacing):
    """VTA volume (mm^3) from a lattice CSV's complex-modulus magnitude.

    The `magnitude` column is in V/mm; `scale` multiplies the unit field
    (dimensionless for a direct run = 1.0). Converted to V/m via x1000.
    """
    df = pd.read_csv(csv_path)
    magn_vpm = df["magnitude"].to_numpy() * abs(scale) * 1000.0
    return float((magn_vpm > THRESHOLD_VPM).sum()) * spacing**3, len(df)


def main():
    """Run the direct simulation and compare the three VTA estimates."""
    ossdbs.set_logger(level=logging.INFO)
    with open(CONFIG_FILE) as fp:
        cfg = json.load(fp)

    spacing = float(cfg["PointModel"]["Lattice"]["PointDistance[mm]"])
    current_a = CURRENT_MA * 1e-3

    # --- direct full simulation: one active contact, others floating ---
    cfg["StimSets"]["Active"] = False
    cfg["CalcAxonActivation"] = False
    cfg["ActivationThresholdVTA[V-per-m]"] = THRESHOLD_VPM  # -> mesh VTA volume
    cfg["OutputPath"] = OUT_DIR
    for contact in cfg["Electrodes"][0]["Contacts"]:
        if contact["Contact_ID"] == CONTACT_1B:
            contact["Active"] = True
            contact["Floating"] = False
            contact["Current[A]"] = current_a
            # non-zero reference voltage so only the ground (BrainSurface,
            # V=0) counts as the single grounded active contact in
            # current-controlled mode (see VCM._check_signal).
            contact["Voltage[V]"] = 1.0
        else:
            contact["Active"] = False
            contact["Floating"] = True
            contact["Current[A]"] = 0.0
    # ground returns the injected current
    cfg["Surfaces"][0]["Current[A]"] = -current_a

    main_run(cfg)

    # --- gather the three estimates ---
    spacing_mm = spacing
    direct_lattice_vol, n_direct = magnitude_vta_volume(
        os.path.join(OUT_DIR, "E_field_Lattice.csv"), 1.0, spacing_mm
    )
    unit_lattice_vol, n_unit = magnitude_vta_volume(
        os.path.join(UNIT_DIR, "E_field_Lattice.csv"), current_a, spacing_mm
    )
    with open(os.path.join(OUT_DIR, "VCM_report.json")) as fp:
        report = json.load(fp)
    mesh_vol = report.get("VTA_volume_mm3")

    print("\n" + "=" * 64)
    print(f"VTA scale validation  (contact {CONTACT_1B}, {CURRENT_MA} mA)")
    print("=" * 64)
    print(f"  NGSolve mesh-integrated VTA (VCM report) : {mesh_vol} mm^3")
    print(
        f"  Direct-sim lattice VTA (|E|>thr)         : {direct_lattice_vol:.2f} mm^3"
        f"  ({n_direct} pts)"
    )
    print(
        f"  StimSets unit x {current_a:+.4f} A lattice VTA : "
        f"{unit_lattice_vol:.2f} mm^3  ({n_unit} pts)"
    )
    if mesh_vol:
        print(
            f"  lattice/mesh ratio (direct)              : "
            f"{direct_lattice_vol / mesh_vol:.3f}"
        )
    print(
        f"  superposition/direct ratio (lattice)     : "
        f"{unit_lattice_vol / direct_lattice_vol:.3f}"
        if direct_lattice_vol
        else "  (direct VTA is zero)"
    )
    box_vol = n_direct * spacing_mm**3
    print(f"  lattice box volume (reference)           : {box_vol:.0f} mm^3")


if __name__ == "__main__":
    main()
