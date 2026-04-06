"""Check impedance decomposition: Z_total ≈ Z_tissue + Z_interface.

Compares the impedance from simulations with and without surface impedance.
The interface impedance is computed from the impedancefitter model and the
contact area. A small discrepancy (~1%) is expected because the Robin BC
modifies the current distribution in the tissue volume.

Usage
-----
First run the required simulations::

    ossdbs input_no_interface.json
    ossdbs input_interface_1kOhm.json
    ossdbs input_interface_1kOhm_EQS.json
    ossdbs input_interface_lempka2009.json
    ossdbs input_no_interface_lempka2009.json

Then run this script::

    python check_impedance.py
"""

import csv
import sys
from pathlib import Path

import impedancefitter as ifit
import numpy as np


def read_impedance(result_dir: str) -> complex:
    """Read impedance from a result directory."""
    path = Path(result_dir) / "impedance.csv"
    with open(path) as f:
        reader = csv.DictReader(f)
        row = next(reader)
        return complex(float(row["real"]), float(row["imag"]))


def check_pair(
    name: str,
    z_total: complex,
    z_tissue: complex,
    z_interface: complex,
    tol: float = 0.05,
) -> bool:
    """Check Z_total ≈ Z_tissue + Z_interface and print results."""
    z_reconstructed = z_tissue + z_interface
    error = abs(z_total - z_reconstructed) / abs(z_total)
    passed = error < tol

    print(f"\n{'=' * 60}")
    print(f"  {name}")
    print(f"{'=' * 60}")
    print(f"  Z_total (with interface):    {z_total:.1f} Ohm")
    print(f"  Z_tissue (without interface): {z_tissue:.1f} Ohm")
    print(f"  Z_interface (from model):     {z_interface:.1f} Ohm")
    print(f"  Z_reconstructed:              {z_reconstructed:.1f} Ohm")
    print(f"  Relative error:               {error:.2%}")
    print(f"  {'PASS' if passed else 'FAIL'} (tolerance: {tol:.0%})")
    return passed


def main():
    """Check all simulation scripts."""
    all_passed = True

    # --- 1 kOhm resistive interface (non-EQS, 100 Hz) ---
    try:
        z_total = read_impedance("Results_1kOhm")
        z_tissue = read_impedance("Results_no_interface")
        # R model: Z = R = 1000 Ohm (frequency-independent)
        ecm = ifit.get_equivalent_circuit_model("R")
        z_interface = complex(ecm.eval(omega=2 * np.pi * 100, R=1e3))
        passed = check_pair(
            "1 kOhm resistive (non-EQS)", z_total, z_tissue, z_interface
        )
        all_passed &= passed
    except FileNotFoundError as e:
        print(f"\nSkipping 1 kOhm check: {e}")

    # --- CPE interface (EQS, 1 kHz) ---
    try:
        z_total = read_impedance("Results_1kOhm_EQS")
        z_tissue = read_impedance("Results_no_interface_EQS")
        ecm = ifit.get_equivalent_circuit_model("CPE")
        z_interface = complex(ecm.eval(omega=2 * np.pi * 1e3, k=1e5, alpha=0.5))
        passed = check_pair("1 kOhm resistive (EQS)", z_total, z_tissue, z_interface)
        all_passed &= passed
    except FileNotFoundError as e:
        print(f"\nSkipping CPE check: {e}")

    # --- Lempka 2009 CPE_dl interface (EQS, 100 Hz) ---
    try:
        z_total = read_impedance("Results_LempkaImpedance")
        z_tissue = read_impedance("Results_NoInterface_Lempka")
        # CPE_dl model at 100 Hz
        ecm = ifit.get_equivalent_circuit_model("CPE_dl")
        z_interface = complex(ecm.eval(omega=2 * np.pi * 100, dl_k=1.5e6, dl_alpha=0.8))
        passed = check_pair(
            "Lempka CPE_dl (EQS, 100 Hz)", z_total, z_tissue, z_interface
        )
        all_passed &= passed
    except FileNotFoundError as e:
        print(f"\nSkipping Lempka check: {e}")

    # --- Floating contacts (current-controlled, R = 1 kOhm) ---
    try:
        fp_path = Path("Results_1kOhm_floating") / "floating_potentials.csv"
        with open(fp_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)
            v1 = float(row["E1C1_real"])
            v2 = float(row["E1C8_real"])
        current = 0.002  # A
        z_total_float = (v1 - v2) / current  # total impedance from floating potentials
        z_tissue_float = read_impedance("Results_no_interface")
        z_interface_float = 2 * 1000.0  # 2 contacts x R=1kOhm each
        passed = check_pair(
            "Floating R=1kOhm (current-controlled)",
            z_total_float,
            z_tissue_float.real,
            z_interface_float,
        )
        all_passed &= passed
    except FileNotFoundError as e:
        print(f"\nSkipping floating check: {e}")

    print(f"\n{'=' * 60}")
    if all_passed:
        print("  All checks passed.")
    else:
        print("  Some checks FAILED.")
    print(f"{'=' * 60}")

    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
