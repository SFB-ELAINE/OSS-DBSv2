import os
import subprocess

import nibabel as nib
import numpy as np
import pandas as pd


def test_impedance(input_dir: str, input_json: str, output_dir: str, desired_dir: str):
    """
    Tests impedance by comparing the actual output with the desired output.

    Parameters
    ----------
    input_dir : str
        Directory where the input files are located.
    input_json : str
        JSON file with input parameters.
    output_dir : str
        Directory where the output files are generated.
    desired_dir : str
        Directory where the desired output files are located.
    """
    os.chdir(input_dir)
    subprocess.run(["ossdbs", input_json], check=True)
    os.chdir(output_dir)
    actual = pd.read_csv("impedance.csv").to_numpy()
    os.chdir(desired_dir)
    desired = pd.read_csv("impedance.csv").to_numpy()
    os.chdir("../../..")
    np.testing.assert_allclose(actual, desired, atol=1e-5)


def test_VTA(output_dir: str, output_nii: str, desired_dir: str, desired_nii: str):
    """
    Tests VTA by comparing the actual output with the desired output
    using the Dice coefficient.

    Parameters
    ----------
    output_dir : str
        Directory where the output files are located.
    output_nii : str
        NIfTI file with the output data.
    desired_dir : str
        Directory where the desired output files are located.
    desired_nii : str
        NIfTI file with the desired data.
    """
    os.chdir(output_dir)
    actual = _readNIfTI(output_nii)
    os.chdir(desired_dir)
    desired = _readNIfTI(desired_nii)
    os.chdir("../../..")
    assert _compute_dice_coefficient(actual, desired) == 1


def _readNIfTI(filename: str) -> set:
    """
    Reads a NIfTI file.

    Parameters
    ----------
    filename : str
        Path to the NIfTI file.
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


def _compute_dice_coefficient(set1, set2) -> float:
    """
    Computes the Dice coefficient.

    Parameters
    ----------
    set1 : set
        First set of points.
    set2 : set
        Second set of points.
    """
    intersection = set1.intersection(set2)
    intersection_size = len(intersection)
    size1 = len(set1)
    size2 = len(set2)
    return (2 * intersection_size) / (size1 + size2)


# input_case1
test_impedance(
    input_dir="input_case1",
    input_json="input_homogeneous.json",
    output_dir="Results_homogeneous",
    desired_dir="../../desired_output/input_case1/Results_homogeneous",
)

test_impedance(
    input_dir="input_case1",
    input_json="input_inhomogeneous.json",
    output_dir="Results_inhomogeneous",
    desired_dir="../../desired_output/input_case1/Results_inhomogeneous",
)

# input_case2
test_impedance(
    input_dir="input_case2",
    input_json="input_custom_electrode.json",
    output_dir="Results_electrode",
    desired_dir="../../desired_output/input_case2/Results_electrode",
)

test_impedance(
    input_dir="input_case2",
    input_json="input_custom_material.json",
    output_dir="Results_material",
    desired_dir="../../desired_output/input_case2/Results_material",
)

# input_case3
test_impedance(
    input_dir="input_case3",
    input_json="input_case_grounding.json",
    output_dir="Results_case_grounding",
    desired_dir="../../desired_output/input_case3/Results_case_grounding",
)

test_impedance(
    input_dir="input_case3",
    input_json="input_case_grounding_EQS.json",
    output_dir="Results_case_grounding_EQS",
    desired_dir="../../desired_output/input_case3/Results_case_grounding_EQS",
)

# input_case4
test_impedance(
    input_dir="input_case4",
    input_json="input_current_controlled.json",
    output_dir="Results_current_controlled",
    desired_dir="../../desired_output/input_case4/Results_current_controlled",
)

# input_case5
test_impedance(
    input_dir="input_case5",
    input_json="input_stimulation_signal.json",
    output_dir="Results_signal",
    desired_dir="../../desired_output/input_case5/Results_signal",
)

# input_case6
test_impedance(
    input_dir="input_case6",
    input_json="input_floating.json",
    output_dir="Results_floating",
    desired_dir="../../desired_output/input_case6/Results_floating",
)

# input_case7
test_impedance(
    input_dir="input_case7",
    input_json="input_vta.json",
    output_dir="Results_VTA",
    desired_dir="../../desired_output/input_case7/Results_VTA",
)

test_impedance(
    input_dir="input_case7",
    input_json="input_vta_out_of_core.json",
    output_dir="Results_VTA_OOC",
    desired_dir="../../desired_output/input_case7/Results_VTA_OOC",
)

test_VTA(
    output_dir="input_case7/Results_VTA",
    output_nii="VTA_solution_Lattice.nii",
    desired_dir="../../desired_output/input_case7/Results_VTA",
    desired_nii="VTA_solution_Lattice.nii",
)

test_VTA(
    output_dir="input_case7/Results_VTA_OOC",
    output_nii="VTA_solution_Lattice.nii",
    desired_dir="../../desired_output/input_case7/Results_VTA_OOC",
    desired_nii="VTA_solution_Lattice.nii",
)

# input_case8
test_impedance(
    input_dir="input_case8",
    input_json="input_pathway_out_of_core.json",
    output_dir="Results_PAM_OOC",
    desired_dir="../../desired_output/input_case8/Results_PAM_OOC",
)

test_impedance(
    input_dir="input_case8",
    input_json="input_pathway.json",
    output_dir="Results_PAM",
    desired_dir="../../desired_output/input_case8/Results_PAM",
)
