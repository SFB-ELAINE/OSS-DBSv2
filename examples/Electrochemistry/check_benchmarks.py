"""Check manual NGSolve scripts and OSS-DBS runs against reference benchmarks.

The reference CSVs in ``benchmark/`` were produced by the manual NGSolve
scripts (``Simulation_hp_EQS.py`` and ``Simulation_hp_EQS_pure_float.py``)
and are used to detect regressions in the surface-impedance formulation,
both for the standalone NGSolve scripts *and* for the corresponding
OSS-DBS JSON runs.

Usage
-----
First (re)run the manual simulations to generate fresh CSVs in the
current directory::

    python Simulation_hp_EQS.py
    python Simulation_hp_EQS_pure_float.py

Then run the OSS-DBS JSON inputs to populate the ``Results_*`` folders::

    ossdbs input_interface_1kOhm_EQS.json
    ossdbs input_interface_1kOhm_floating.json

Then run this script::

    python check_benchmarks.py

Numerical comparisons use a relative tolerance because mesh randomness in
NGSolve produces slightly different DOFs and therefore slightly different
physical quantities across runs. OSS-DBS and the manual scripts also use
opposite sign conventions for contact currents, so the comparison is
sign-invariant (the smaller of ``|c - b|`` and ``|c + b|`` is used).
"""

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
BENCHMARK_DIR = HERE / "benchmark"

# Tolerances shared across both comparison paths
TOL_MANUAL = 5e-2
TOL_OSSDBS = 5e-2

# Comparison spec for the manual NGSolve script outputs that live next to
# this file. Each column listed here is looked up directly in both CSVs.
MANUAL_BENCHMARKS = [
    {
        "name": "Simulation_hp_EQS.py",
        "file": "results_Z_EQS.csv",
        "columns": ["I", "I1", "I2", "Field"],
    },
    {
        "name": "Simulation_hp_EQS_pure_float.py",
        "file": "results_Z_pure_float.csv",
        "columns": ["I1", "I2", "V1", "V2", "Field"],
    },
]


def _to_complex(value):
    """Parse a CSV cell that may be a real number or a complex string."""
    if isinstance(value, complex):
        return value
    if isinstance(value, int | float):
        return complex(value)
    return complex(str(value).replace(" ", ""))


def _rel_err(bench: complex, current: complex, sign_invariant: bool) -> float:
    """Relative error between two (complex) values.

    When ``sign_invariant`` is true the smaller of ``|c - b|`` and
    ``|c + b|`` is used; this absorbs opposite sign conventions between
    OSS-DBS and the manual NGSolve scripts.
    """
    denom = abs(bench) if abs(bench) > 0 else 1.0
    diff = abs(current - bench)
    if sign_invariant:
        diff = min(diff, abs(current + bench))
    return diff / denom


def _print_row(label: str, bench: complex, current: complex, tol: float) -> bool:
    rel_err = _rel_err(bench, current, sign_invariant=True)
    passed = rel_err < tol
    marker = "PASS" if passed else "FAIL"
    print(
        f"  [{label:4s}] benchmark={bench!s:>32s}  "
        f"current={current!s:>32s}  rel_err={rel_err:7.2%}  "
        f"(tol={tol:.0%})  {marker}"
    )
    return passed


def _load_benchmark(file_name: str) -> pd.Series:
    path = BENCHMARK_DIR / file_name
    if not path.is_file():
        raise FileNotFoundError(f"Missing benchmark file: {path}")
    return pd.read_csv(path).iloc[0]


# ---------------------------------------------------------------------------
# Manual NGSolve scripts: compare produced CSV against benchmark CSV.
# ---------------------------------------------------------------------------


def check_manual(name: str, file_name: str, columns: list[str]) -> bool:
    """Check manual simulation files."""
    print(f"\n{'=' * 72}")
    print(f"  Manual script: {name}  ({file_name})")
    print(f"{'=' * 72}")

    try:
        bench = _load_benchmark(file_name)
    except FileNotFoundError as e:
        print(f"  {e}")
        return False

    current_path = HERE / file_name
    if not current_path.is_file():
        print(f"  MISSING current result: {current_path}")
        print(f"  Re-run the simulation: python {Path(name).name}")
        return False

    current = pd.read_csv(current_path).iloc[0]

    all_passed = True
    for col in columns:
        if col not in bench.index or col not in current.index:
            print(f"  [{col}] MISSING column - FAIL")
            all_passed = False
            continue
        b = _to_complex(bench[col])
        c = _to_complex(current[col])
        all_passed &= _print_row(col, b, c, TOL_MANUAL)
    return all_passed


# ---------------------------------------------------------------------------
# OSS-DBS runs: assemble the same physical quantities from Results_*/ and
# compare them against the manual-script benchmark CSV.
# ---------------------------------------------------------------------------


def _read_currents_mA(result_dir: Path, contact: str) -> complex:
    """Return per-contact current in mA from ``currents.csv``."""
    df = pd.read_csv(result_dir / "currents.csv").iloc[0]
    return 1e3 * complex(df[f"{contact}_real"], df[f"{contact}_imag"])


def _read_impedance(result_dir: Path) -> complex:
    df = pd.read_csv(result_dir / "impedance.csv").iloc[0]
    return complex(df["real"], df["imag"])


def _read_floating_potential(result_dir: Path, contact: str) -> complex:
    df = pd.read_csv(result_dir / "floating_potentials.csv").iloc[0]
    return complex(df[f"{contact}_real"], df[f"{contact}_imag"])


def check_ossdbs_eqs() -> bool:
    """Compare ``Results_1kOhm_EQS/`` against ``benchmark/results_Z_EQS.csv``.

    The manual script applies ``V=1 V`` between ``Contact_1`` and
    ``Contact_8`` with a CPE Robin BC. OSS-DBS reports the per-contact
    currents in Amperes (``currents.csv``) and the total impedance in Ohms
    (``impedance.csv``). We convert both to mA and compute the equivalent
    total current ``I = V / Z`` to compare against the benchmark columns
    ``I``, ``I1``, ``I2``.
    """
    name = "input_interface_1kOhm_EQS.json"
    result_dir = HERE / "Results_1kOhm_EQS"

    print(f"\n{'=' * 72}")
    print(f"  OSS-DBS: {name}  ({result_dir.name}/ vs results_Z_EQS.csv)")
    print(f"{'=' * 72}")

    try:
        bench = _load_benchmark("results_Z_EQS.csv")
    except FileNotFoundError as e:
        print(f"  {e}")
        return False

    if not result_dir.is_dir():
        print(f"  MISSING OSS-DBS result dir: {result_dir}")
        print(f"  Run: ossdbs {name}")
        return False

    try:
        Z = _read_impedance(result_dir)
        I1 = _read_currents_mA(result_dir, "E1C1")
        I2 = _read_currents_mA(result_dir, "E1C8")
    except (FileNotFoundError, KeyError) as e:
        print(f"  Failed to read OSS-DBS results: {e}")
        return False

    # Applied voltage is 1 V; impedance is in Ohm, so I_total [mA] = 1e3 / Z.
    I_total = 1e3 / Z if abs(Z) > 0 else complex("nan")

    current_values = {"I": I_total, "I1": I1, "I2": I2}
    all_passed = True
    for col, c in current_values.items():
        b = _to_complex(bench[col])
        all_passed &= _print_row(col, b, c, TOL_OSSDBS)
    return all_passed


def check_ossdbs_floating() -> bool:
    """Compare ``Results_1kOhm_floating/`` against ``results_Z_pure_float.csv``.

    The manual script drives 2 mA between ``Contact_1`` and ``Contact_8``
    with a 1 kOhm Robin BC on each contact and records the resulting
    floating potentials and contact currents. OSS-DBS reports the
    floating potentials in Volts (``floating_potentials.csv``) and the
    per-contact currents in Amperes (``currents.csv``); we compare the
    ``V1``, ``V2``, ``I1``, ``I2`` columns of the benchmark.
    """
    name = "input_interface_1kOhm_floating.json"
    result_dir = HERE / "Results_1kOhm_floating"

    print(f"\n{'=' * 72}")
    print(f"  OSS-DBS: {name}  ({result_dir.name}/ vs results_Z_pure_float.csv)")
    print(f"{'=' * 72}")

    try:
        bench = _load_benchmark("results_Z_pure_float.csv")
    except FileNotFoundError as e:
        print(f"  {e}")
        return False

    if not result_dir.is_dir():
        print(f"  MISSING OSS-DBS result dir: {result_dir}")
        print(f"  Run: ossdbs {name}")
        return False

    try:
        V1 = _read_floating_potential(result_dir, "E1C1")
        V2 = _read_floating_potential(result_dir, "E1C8")
        I1 = _read_currents_mA(result_dir, "E1C1")
        I2 = _read_currents_mA(result_dir, "E1C8")
    except (FileNotFoundError, KeyError) as e:
        print(f"  Failed to read OSS-DBS results: {e}")
        return False

    current_values = {"V1": V1, "V2": V2, "I1": I1, "I2": I2}
    all_passed = True
    for col, c in current_values.items():
        b = _to_complex(bench[col])
        all_passed &= _print_row(col, b, c, TOL_OSSDBS)
    return all_passed


def main():
    """Run every benchmark comparison and exit non-zero on any failure."""
    all_passed = True

    # Manual NGSolve scripts
    for spec in MANUAL_BENCHMARKS:
        all_passed &= check_manual(spec["name"], spec["file"], spec["columns"])

    # OSS-DBS runs against the same benchmarks
    all_passed &= check_ossdbs_eqs()
    all_passed &= check_ossdbs_floating()

    print(f"\n{'=' * 72}")
    if all_passed:
        print("  All benchmark checks passed.")
    else:
        print("  Some benchmark checks FAILED.")
    print(f"{'=' * 72}")

    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
