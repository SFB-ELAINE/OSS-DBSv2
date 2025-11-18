"""
Example masking a DTI image based on the material distribution.

Outcome:
    Returns anisotropy tensors with and without masking for
    randomly selected points.
    Note that if MRI or DTI are well aligned or not, either
    almost none or almost all points will be masked.
"""

import random

import numpy as np

from ossdbs.api import generate_mesh, prepare_dielectric_properties
from ossdbs.fem.volume_conductor.conductivity import ConductivityCF
from ossdbs.utils.nifti1image import DiffusionTensorImage, MagneticResonanceImage

# Load settings for the simulation
settings = {
    "BrainRegion": {
        "Center": {"x[mm]": -11.23, "y[mm]": -16.46, "z[mm]": -10.01},
        "Dimension": {"x[mm]": 50.0, "y[mm]": 50.0, "z[mm]": 50.0},
        "Shape": "Box",
    },
    "Electrodes": [
        {
            "Name": "BostonScientificVerciseDirected",
            "CustomParameters": None,
            "Rotation[Degrees]": 0.0,
            "Direction": {"x[mm]": -0.40, "y[mm]": 0.67, "z[mm]": 0.62},
            "TipPosition": {"x[mm]": -11.23, "y[mm]": -16.46, "z[mm]": -10.01},
        }
    ],
    "MaterialDistribution": {
        "MRIPath": "../../input_files/Butenko_segmask.nii.gz",
        "MRIMapping": {
            "Unknown": 0,
            "CSF": 1,
            "White matter": 2,
            "Gray matter": 3,
            "Blood": 4,
        },
        "DiffusionTensorActive": True,
        "DTIPath": "../../input_files/IITMeanTensor_NormMapping.nii.gz",
    },
    "DielectricModel": {"Type": "ColeCole4", "CustomParameters": None},
    "Mesh": {
        "LoadMesh": False,
        "LoadPath": "",
        "MeshingHypothesis": {
            "Type": "Default",
            "MaxMeshSize": 1000000.0,
            "MeshSizeFilename": "",
        },
        "SaveMesh": False,
    },
    "ExportElectrode": False,
}

# Load MRI image and create brain geometry
print("Loading MRI image...")
mri_image = MagneticResonanceImage(settings["MaterialDistribution"]["MRIPath"])

# Loading DTI image
dti_image = DiffusionTensorImage(settings["MaterialDistribution"]["DTIPath"])

# Prepare dielectric properties
print("Preparing dielectric properties...")
dielectric_properties = prepare_dielectric_properties(settings)

# Generate mesh for the simulation
print("Generating mesh...")
mesh = generate_mesh(settings)
ngmesh = mesh.ngsolvemesh

# Create conductivity distributions (with and without masking)
print("Creating conductivity distributions...")
conductivity = ConductivityCF(
    mri_image=mri_image,
    brain_bounding_box=mri_image.bounding_box,
    dielectric_properties=dielectric_properties,
    materials=settings["MaterialDistribution"]["MRIMapping"],
    dti_image=dti_image,
    wm_masking=False,
)
conductivity_distribution = conductivity(mesh=mesh, frequency=10000.0)

conductivity_masked = ConductivityCF(
    mri_image=mri_image,
    brain_bounding_box=mri_image.bounding_box,
    dielectric_properties=dielectric_properties,
    materials=settings["MaterialDistribution"]["MRIMapping"],
    dti_image=dti_image,
    wm_masking=True,
)
conductivity_distribution_masked = conductivity_masked(mesh=mesh, frequency=10000.0)

# Generate random test points within the brain region
print("Generating random test points...")
test_indices = [
    (
        round(
            random.uniform(
                settings["BrainRegion"]["Center"]["x[mm]"]
                - settings["BrainRegion"]["Dimension"]["x[mm]"] / 2,
                settings["BrainRegion"]["Center"]["x[mm]"]
                + settings["BrainRegion"]["Dimension"]["x[mm]"] / 2,
            ),
            3,
        ),
        round(
            random.uniform(
                settings["BrainRegion"]["Center"]["y[mm]"]
                - settings["BrainRegion"]["Dimension"]["y[mm]"] / 2,
                settings["BrainRegion"]["Center"]["y[mm]"]
                + settings["BrainRegion"]["Dimension"]["y[mm]"] / 2,
            ),
            3,
        ),
        round(
            random.uniform(
                settings["BrainRegion"]["Center"]["z[mm]"]
                - settings["BrainRegion"]["Dimension"]["z[mm]"] / 2,
                settings["BrainRegion"]["Center"]["z[mm]"]
                + settings["BrainRegion"]["Dimension"]["z[mm]"] / 2,
            ),
            3,
        ),
    )
    for i in range(20)
]

# Compare conductivity tensors at test points
print("Comparing conductivity tensors...")
same_count = 0
diff_count = 0

for point in test_indices:
    try:
        sigma = conductivity_distribution(ngmesh(point))
        sigma_masked = conductivity_distribution_masked(ngmesh(point))
        if np.allclose(sigma, sigma_masked, rtol=1e-6, atol=1e-12):
            same_count += 1
        else:
            diff_count += 1
            print(
                f"Point {point} differs (off-diagonal scaling in CSF and GM "
                f"for sigma but not sigma_masked):"
            )
            print(f"  sigma:\n {np.round(sigma, 4)}")
            print(f"  sigma_masked:\n {np.round(sigma_masked, 4)}")
    except Exception as e:
        print(f"Cannot process point {point}: {e}")

# Print summary of results
print("\nSummary:")
print(f"  Points where sigma == sigma_masked: {same_count}")
print(f"  Points where sigma != sigma_masked: {diff_count}")
