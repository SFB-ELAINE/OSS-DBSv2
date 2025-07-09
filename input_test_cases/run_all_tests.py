import json
import logging
import os

import nibabel as nib
import numpy as np
import pandas as pd

from ossdbs.main import main_run

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)

INPUT_DATA = [
    {
        "input_dir": "input_case1",
        "input_json": "input_case1/input_homogeneous.json",
        "output_csv": "input_case1/Results_homogeneous/impedance.csv",
        "desired_csv": "desired_output/input_case1/Results_homogeneous/impedance.csv",
    },
    {
        "input_dir": "input_case1",
        "input_json": "input_case1/input_inhomogeneous.json",
        "output_csv": "input_case1/Results_inhomogeneous/impedance.csv",
        "desired_csv": "desired_output/input_case1/Results_inhomogeneous/impedance.csv",
    },
    {
        "input_dir": "input_case3",
        "input_json": "input_case3/input_case_grounding.json",
        "output_csv": "input_case3/Results_case_grounding/impedance.csv",
        "desired_csv": "desired_output/input_case3/Results_case_grounding/\
impedance.csv",
    },
    {
        "input_dir": "input_case3",
        "input_json": "input_case3/input_case_grounding_EQS.json",
        "output_csv": "input_case3/Results_case_grounding_EQS/impedance.csv",
        "desired_csv": "desired_output/input_case3/Results_case_grounding_EQS/\
impedance.csv",
    },
    {
        "input_dir": "input_case4",
        "input_json": "input_case4/input_current_controlled.json",
        "output_csv": "input_case4/Results_current_controlled/impedance.csv",
        "desired_csv": "desired_output/input_case4/Results_current_controlled/\
impedance.csv",
    },
    {
        "input_dir": "input_case5",
        "input_json": "input_case5/input_stimulation_signal.json",
        "output_csv": "input_case5/Results_signal/impedance.csv",
        "desired_csv": "desired_output/input_case5/Results_signal/impedance.csv",
    },
    {
        "input_dir": "input_case6",
        "input_json": "input_case6/input_floating.json",
        "output_csv": "input_case6/Results_floating/impedance.csv",
        "desired_csv": "desired_output/input_case6/Results_floating/impedance.csv",
    },
    {
        "input_dir": "input_case7",
        "input_json": "input_case7/input_vta.json",
        "output_csv": "input_case7/Results_VTA/impedance.csv",
        "desired_csv": "desired_output/input_case7/Results_VTA/impedance.csv",
        "output_nii": "input_case7/Results_VTA/VTA_solution_Lattice.nii",
        "desired_nii": "desired_output/input_case7/Results_VTA/\
VTA_solution_Lattice.nii",
    },
    {
        "input_dir": "input_case7",
        "input_json": "input_case7/input_vta_out_of_core.json",
        "output_csv": "input_case7/Results_VTA_OOC/impedance.csv",
        "desired_csv": "desired_output/input_case7/Results_VTA_OOC/impedance.csv",
        "output_nii": "input_case7/Results_VTA_OOC/VTA_solution_Lattice.nii",
        "desired_nii": "desired_output/input_case7/Results_VTA_OOC/\
VTA_solution_Lattice.nii",
    },
    {
        "input_dir": "input_case8",
        "input_json": "input_case8/input_pathway_out_of_core.json",
        "output_csv": "input_case8/Results_PAM_OOC/impedance.csv",
        "desired_csv": "desired_output/input_case8/Results_PAM_OOC/impedance.csv",
    },
    {
        "input_dir": "input_case8",
        "input_json": "input_case8/input_pathway.json",
        "output_csv": "input_case8/Results_PAM/impedance.csv",
        "desired_csv": "desired_output/input_case8/Results_PAM/impedance.csv",
    },
    {
        "input_dir": "input_case2",
        "input_json": "input_case2/input_custom_electrode.json",
        "output_csv": "input_case2/Results_electrode/impedance.csv",
        "desired_csv": "desired_output/input_case2/Results_electrode/impedance.csv",
    },
    {
        "input_dir": "input_case2",
        "input_json": "input_case2/input_custom_material.json",
        "output_csv": "input_case2/Results_material/impedance.csv",
        "desired_csv": "desired_output/input_case2/Results_material/impedance.csv",
    },
]

# Tolerances for test comparisons
IMPEDANCE_ATOL = 1e-2
VTA_ATOL = 1e-2


def run_test() -> pd.DataFrame:
    """
    Runs all test cases defined in INPUT_DATA:
    - Runs simulations for each input JSON.
    - Compares output impedance and VTA results to desired outputs.
    - Collects and logs results for each test.
    - Writes a summary CSV file with the results.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the results of all tests.
    """
    impedance_test_passed = True
    vta_test_passed = True
    test_results = []

    # run all simulations
    for data in INPUT_DATA:
        logging.info("Running %s", data["input_json"])
        logging.info("####################")
        _run_simulation(data["input_dir"], data["input_json"])

    # test everything
    for data in INPUT_DATA:
        tests_appended = False
        try:
            # check impedances
            if "desired_csv" in data:
                logging.info("Testing impedance of %s", data["input_json"])
                impedance_test_passed = _test_impedance(
                    data["output_csv"], data["desired_csv"]
                )
                test_results.append(
                    {
                        "input_file": data["input_json"],
                        "impedance_test": "Passed"
                        if impedance_test_passed
                        else "Failed",
                        "VTA_test": "-",
                    }
                )
                tests_appended = True

            if "desired_nii" in data:
                logging.info("Testing VTA of %s", data["input_json"])
                vta_test_passed = _test_VTA(data["output_nii"], data["desired_nii"])

                if tests_appended:
                    test_results[-1]["VTA_test"] = (
                        "Passed" if vta_test_passed else "Failed"
                    )
                else:
                    test_results.append(
                        {
                            "input_file": data["input_json"],
                            "impedance_test": "-",
                            "VTA_test": "Passed" if vta_test_passed else "Failed",
                        }
                    )
        except Exception as e:
            test_results.append(
                {
                    "input_file": data["input_json"],
                    "impedance_test": "Error",
                    "VTA_test": "Error",
                    "error": str(e),
                }
            )

    df = pd.DataFrame(test_results)
    output_csv_file = "test_input_cases_results.csv"
    df.to_csv(output_csv_file, index=False)
    return df


def _run_simulation(input_dir: str, input_json: str):
    """
    Runs a simulation with a given JSON file and updates input paths as needed.

    Parameters
    ----------
    input_dir : str
        Directory where the input files are located.
    input_json : str
        Path to the JSON file with input parameters.
    """
    with open(input_json) as json_file:
        input_settings = json.load(json_file)

    # add the stimulation folder
    input_settings["StimulationFolder"] = os.path.dirname(input_dir)

    # update folder with MRI
    input_settings["MaterialDistribution"]["MRIPath"] = input_settings[
        "MaterialDistribution"
    ]["MRIPath"].replace("../..", "..")

    # update OutputPath
    input_settings["OutputPath"] = os.path.join(input_dir, input_settings["OutputPath"])

    # update Pathway FileName
    if input_settings["PointModel"]["Pathway"]["FileName"]:
        input_settings["PointModel"]["Pathway"]["FileName"] = input_settings[
            "PointModel"
        ]["Pathway"]["FileName"].replace("../..", "..")

    main_run(input_settings)


def _test_impedance(output_csv: str, desired_csv: str) -> bool:
    """
    Tests impedance by comparing the actual output with the desired output.

    Parameters
    ----------
    output_csv : str
        Path to CSV file with actual output values.
    desired_csv : str
        Path to CSV file with desired output values.

    Returns
    -------
    bool
        True if actual and desired values are close within tolerance, else False.
    """
    with open(output_csv) as csv_file:
        actual = pd.read_csv(csv_file).to_numpy()

    with open(desired_csv) as csv_file:
        desired = pd.read_csv(csv_file).to_numpy()
    return np.allclose(actual, desired, atol=IMPEDANCE_ATOL)


def _test_VTA(output_nii: str, desired_nii: str) -> bool:
    """
    Tests VTA by comparing the actual output with the desired output
    using the Dice coefficient.

    Parameters
    ----------
    output_nii : str
        Path to NIfTI file with the output data.
    desired_nii : str
        Path to NIfTI file with the desired data.

    Returns
    -------
    bool
        True if Dice coefficient is close to 1, else False.
    """
    actual = _readNIfTI(output_nii)
    desired = _readNIfTI(desired_nii)
    return np.isclose(_compute_dice_coefficient(actual, desired), 1, atol=VTA_ATOL)


def _readNIfTI(filename: str) -> set:
    """
    Reads a NIfTI file and returns a set of coordinates for voxels above threshold.

    Parameters
    ----------
    filename : str
        Path to the NIfTI file.

    Returns
    -------
    set
        Set of (X, Y, Z) coordinates for voxels with value > 0.1.
    """
    image = nib.load(filename)
    affine = image.affine
    data_shape = image.get_fdata().shape
    data = image.get_fdata()
    vta_points = []

    end_voxel = [data_shape[0], data_shape[1], data_shape[2]]
    counter = 0
    for i in range(0, end_voxel[0]):
        for j in range(0, end_voxel[1]):
            for k in range(0, end_voxel[2]):
                if data[i, j, k] > 0.1:
                    counter += 1
                    X_coord = affine[0][3] + i * affine[0][0]
                    Y_coord = affine[1][3] + j * affine[1][1]
                    Z_coord = affine[2][3] + k * affine[2][2]
                    vta_points.append((X_coord, Y_coord, Z_coord))

    return set(vta_points)


def _compute_dice_coefficient(set1: set, set2: set) -> float:
    """
    Computes the Dice coefficient between two sets.

    Parameters
    ----------
    set1 : set
        First set of points.
    set2 : set
        Second set of points.

    Returns
    -------
    float
        Dice coefficient between set1 and set2.
    """
    intersection = set1.intersection(set2)
    intersection_size = len(intersection)
    size1 = len(set1)
    size2 = len(set2)
    return (2 * intersection_size) / (size1 + size2)


def check_tests(df: pd.DataFrame) -> bool:
    """
    Check if all tests were successful.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing test results.

    Returns
    -------
    bool
        True if all tests passed, False otherwise.
    """
    # Print a concise summary table at the end
    summary_cols = ["input_file", "impedance_test", "VTA_test"]
    summary = df[summary_cols] if all(col in df.columns for col in summary_cols) else df
    logging.info("\n" + "=" * 60)
    logging.info("Test Summary:")
    logging.info("\n" + summary.to_string(index=False) + "\n" + "=" * 60)

    # Log failed tests only once
    failed = df[(df["impedance_test"] == "Failed") | (df["VTA_test"] == "Failed")]
    if not failed.empty:
        for _, row in failed.iterrows():
            logging.error(
                f"Test {row.get('input_file', row.get('input_file', ''))} failed | "
                f"Impedance: {row.get('impedance_test')} | VTA: {row.get('VTA_test')}"
            )
        return False
    return True


if __name__ == "__main__":
    tests_df = run_test()
    all_tests_passed = check_tests(tests_df)
    if not all_tests_passed:
        open("fail_oss.txt", "w").close()
