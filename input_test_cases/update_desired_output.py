"""Update desired_output/ from current simulation results.

Run from input_test_cases/ after running simulations to refresh baselines:
    cd input_test_cases && python update_desired_output.py
"""

import os
import shutil

# Files to copy from Results_* directories when present
COPY_FILES = [
    "impedance.csv",
    "VTA_solution_Lattice.nii.gz",
    "Lattice.json",
    "floating_potentials.csv",
    "impedance_matrix.csv",
    "currents.csv",
]

for test_dir in sorted(os.listdir(".")):
    if not os.path.isdir(test_dir) or "input_case" not in test_dir:
        continue
    for subdir in sorted(os.listdir(test_dir)):
        source_dir = os.path.join(test_dir, subdir)
        if not os.path.isdir(source_dir) or "Result" not in subdir:
            continue
        result_dir = os.path.join("desired_output", test_dir, subdir)
        if not os.path.isdir(result_dir):
            os.makedirs(result_dir, exist_ok=True)
        copied = []
        for filename in COPY_FILES:
            src = os.path.join(source_dir, filename)
            if os.path.exists(src):
                shutil.copy(src, os.path.join(result_dir, filename))
                copied.append(filename)
        # Also copy Pathway_status*.json files (PAM results)
        for filename in os.listdir(source_dir):
            if filename.startswith("Pathway_status") and filename.endswith(".json"):
                shutil.copy(
                    os.path.join(source_dir, filename),
                    os.path.join(result_dir, filename),
                )
                copied.append(filename)
        if copied:
            print(f"{result_dir}: {', '.join(copied)}")
