"""
The interface with the BrainGeometry
and the respective entries from the input
dictionary are shown.
The different choices of modelling
the brain are shown.
The different geometries are saved
to STEP files.
"""

import netgen.occ as occ
from ngsolve import Draw, Mesh, TaskManager

import ossdbs

settings = {
    "MaterialDistribution": {"MRIPath": "segmask.nii.gz"},
    "BrainRegion": {
        "Center": {"x[mm]": 5, "y[mm]": 14, "z[mm]": -4.5},
        "Dimension": {"x[mm]": 50.0, "y[mm]": 50.0, "z[mm]": 50.0},
    },
}

nifti_image = ossdbs.MagneticResonanceImage(settings["MaterialDistribution"]["MRIPath"])
box = ossdbs.BrainGeometry("Box", nifti_image.bounding_box)
box.geometry.WriteStep("box.step")
region_parameters = settings["BrainRegion"]
region_of_interest = ossdbs.create_bounding_box(region_parameters)
box = ossdbs.BrainGeometry("Box", region_of_interest)
box.geometry.WriteStep("box_ROI.step")

geometry = ossdbs.BrainGeometry("Ellipsoid", region_of_interest)
geometry.geometry.WriteStep("geo_ROI_ellipsoid.step")
geometry = ossdbs.BrainGeometry("Ellipsoid", nifti_image.bounding_box)
geometry.geometry.WriteStep("geo_full_ellipsoid.step")

geometry = ossdbs.BrainGeometry("Sphere", region_of_interest)
geometry.geometry.WriteStep("geo_ROI_sphere.step")
geometry = ossdbs.BrainGeometry("Sphere", nifti_image.bounding_box)
geometry.geometry.WriteStep("geo_full_sphere.step")

occgeo = occ.OCCGeometry(geometry.geometry)
with TaskManager():
    mesh = Mesh(occgeo.GenerateMesh())
Draw(mesh)
print(mesh.GetBoundaries())
print(mesh.GetMaterials())
