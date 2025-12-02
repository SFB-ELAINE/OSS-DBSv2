"""
Example masking a DTI image based on the material distribution.
This option can be selected to make sure that anisotropy tensors
are only considered in white matter regions.

Outcome:
    Exports two conductivity fields as VTK (.vtu) files:
        conductivity.vtu          # unmasked
        conductivity_masked.vtu   # masked

Visualization in ParaView:
    Load both files, select them together, and apply:
        Filters -> Alphabetical -> Append Attributes
    Then apply a Python Calculator to the Append Attributes output with:
        conductivity_masked_real - conductivity_real

    Negative values in the result indicate voxels where the unmasked conductivity
    had unphysically high/scaled conductivity compared to the masked conductivty.
"""

import os

import ngsolve
import numpy as np

from ossdbs.api import (
    create_bounding_box,
    generate_mesh,
    load_images,
    prepare_dielectric_properties,
)
from ossdbs.fem.volume_conductor.conductivity import ConductivityCF
from ossdbs.utils.vtk_export import FieldSolution

# Load settings for the simulation
settings = {
    "BrainRegion": {
        "Center": {"x[mm]": -13.99, "y[mm]": -7.73, "z[mm]": -7.91},
        "Dimension": {"x[mm]": 50.0, "y[mm]": 50.0, "z[mm]": 50.0},
        "Shape": "Sphere",
    },
    "Electrodes": [
        {
            "Name": "BostonScientificVerciseDirected",
            "CustomParameters": None,
            "Rotation[Degrees]": 0.0,
            "Direction": {"x[mm]": -0.45, "y[mm]": 0.65, "z[mm]": 0.61},
            "TipPosition": {"x[mm]": -13.99, "y[mm]": -7.73, "z[mm]": -7.91},
        }
    ],
    "MaterialDistribution": {
        "MRIPath": "../../input_files/sub-John_Doe/JD_segmask.nii.gz",
        "MRIMapping": {
            "Unknown": 0,
            "CSF": 3,
            "White matter": 2,
            "Gray matter": 1,
            "Blood": 4,
        },
        "DiffusionTensorActive": True,
        "DTIPath": "../../input_files/sub-John_Doe/JD_DTI_NormMapping.nii.gz",
    },
    "DielectricModel": {"Type": "ColeCole4", "CustomParameters": None},
    "Mesh": {
        "LoadMesh": False,
        "LoadPath": "",
        "MeshingHypothesis": {
            "Type": "Fine",
            "MaxMeshSize": 1000.0,
            "MeshSizeFilename": "",
        },
        "SaveMesh": False,
    },
    "OutputPath": "./",
}

# Load images
print("Load images...")
mri_image, dti_image = load_images(settings)

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
    brain_bounding_box=create_bounding_box(settings["BrainRegion"]),
    dielectric_properties=dielectric_properties,
    materials=settings["MaterialDistribution"]["MRIMapping"],
    dti_image=dti_image,
    wm_masking=False,
)
conductivity_distribution = conductivity(mesh=mesh, frequency=10000.0)

# Naming convention by ParaView!
cf_list = (
    conductivity_distribution[0],  # xx
    conductivity_distribution[4],  # yy
    conductivity_distribution[8],  # zz
    conductivity_distribution[1],  # xy
    conductivity_distribution[5],  # yz
    conductivity_distribution[2],  # xz
)
conductivity_export = ngsolve.CoefficientFunction(cf_list, dims=(6,))
FieldSolution(conductivity_export, "conductivity", ngmesh, False).save(
    os.path.join(settings["OutputPath"], "conductivity")
)

conductivity_masked = ConductivityCF(
    mri_image=mri_image,
    brain_bounding_box=create_bounding_box(settings["BrainRegion"]),
    dielectric_properties=dielectric_properties,
    materials=settings["MaterialDistribution"]["MRIMapping"],
    dti_image=dti_image,
    wm_masking=True,
)
conductivity_distribution_masked = conductivity_masked(mesh=mesh, frequency=10000.0)

# Naming convention by ParaView!
cf_list = (
    conductivity_distribution_masked[0],  # xx
    conductivity_distribution_masked[4],  # yy
    conductivity_distribution_masked[8],  # zz
    conductivity_distribution_masked[1],  # xy
    conductivity_distribution_masked[5],  # yz
    conductivity_distribution_masked[2],  # xz
)
conductivity_export_masked = ngsolve.CoefficientFunction(cf_list, dims=(6,))
FieldSolution(conductivity_export_masked, "conductivity_masked", ngmesh, False).save(
    os.path.join(settings["OutputPath"], "conductivity_masked")
)

FieldSolution(conductivity.material_distribution(mesh), "material", ngmesh, False).save(
    os.path.join(settings["OutputPath"], "material")
)

print("Stored conductivity distributions in OutputPath.")

test_point = (-14.14, -1.63, -16.99)

sigma = conductivity_distribution(ngmesh(test_point))
sigma_masked = conductivity_distribution_masked(ngmesh(test_point))

print(
    "Conductivity at vertices of the element including point",
    test_point,
    "without masking:",
)
print("WM:\n", np.round(np.array(sigma[0], dtype=float).reshape((3, 3)), 8))
print("CSF:\n", np.round(np.array(sigma[1], dtype=float).reshape((3, 3)), 8))
print("GM:\n", np.round(np.array(sigma[2], dtype=float).reshape((3, 3)), 8))

print(
    "Conductivity at vertices of the element including point",
    test_point,
    "with masking:",
)
print("WM:\n", np.round(np.array(sigma_masked[0], dtype=float).reshape((3, 3)), 8))
print("CSF:\n", np.round(np.array(sigma_masked[1], dtype=float).reshape((3, 3)), 8))
print("GM:\n", np.round(np.array(sigma_masked[2], dtype=float).reshape((3, 3)), 8))

print(
    "No masking was applied to WM, CSF tensor remains isotropic (masking not needed), "
    "and masking was applied to GM."
)
