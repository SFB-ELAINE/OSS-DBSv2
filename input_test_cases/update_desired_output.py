import os
import shutil

for test_dir in os.listdir("."):
    if not os.path.isdir(test_dir):
        continue
    if "input_case" in test_dir:
        subdirs = os.listdir(test_dir)
        for subdir in subdirs:
            source_dir = os.path.join(test_dir, subdir)
            if not os.path.isdir(source_dir):
                continue
            if "Result" in subdir:
                result_dir = os.path.join("desired_output", test_dir, subdir)
                if not os.path.isdir(result_dir):
                    continue
                print(result_dir)
                shutil.copy(
                    os.path.join(source_dir, "impedance.csv"),
                    os.path.join(result_dir, "impedance.csv"),
                )
                try:
                    shutil.copy(
                        os.path.join(source_dir, "VTA_solution_Lattice.nii"),
                        os.path.join(result_dir, "VTA_solution_Lattice.nii"),
                    )
                except FileNotFoundError:
                    # if file is not there, we don't need to copy it
                    pass
