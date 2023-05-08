"""
The interface with the BrainGeometry
and the respective entries from the input
dictionary are shown.
The different choices of modelling
the brain are shown.
The different geometries are saved
to STEP files.
"""
import ossdbs
from ossdbs.factories import BrainRegionFactory
from ngsolve import Draw, Mesh, TaskManager
import netgen.occ as occ

settings = \
    {"MaterialDistribution":
        {"MRIPath": "homogeneous.nii"
         },
     "BrainRegion":
        {"Center": {"x[mm]": 5, "y[mm]": 14, "z[mm]": -4.5},
         "Dimension": {"x[mm]": 50.0, "y[mm]": 50.0, "z[mm]": 50.0}
         }
     }

nifti_image = ossdbs.Nifti1Image(settings['MaterialDistribution']['MRIPath'])
box = ossdbs.BrainGeometry(nifti_image.bounding_box(), "Box")
box.geometry.WriteStep("box.step")
region_parameters = settings['BrainRegion']
region_of_interest = BrainRegionFactory.create(region_parameters)
box = ossdbs.BrainGeometry(region_of_interest, "Box")
box.geometry.WriteStep("box_ROI.step")

geometry = ossdbs.BrainGeometry(region_of_interest, "Ellipsoid")
geometry.geometry.WriteStep("geo_ROI_ellipsoid.step")
geometry = ossdbs.BrainGeometry(nifti_image.bounding_box(), "Ellipsoid")
geometry.geometry.WriteStep("geo_full_ellipsoid.step")

geometry = ossdbs.BrainGeometry(region_of_interest, "Sphere")
geometry.geometry.WriteStep("geo_ROI_sphere.step")
geometry = ossdbs.BrainGeometry(nifti_image.bounding_box(), "Sphere")
geometry.geometry.WriteStep("geo_full_sphere.step")

occgeo = occ.OCCGeometry(geometry.geometry)
with TaskManager():
    mesh = Mesh(occgeo.GenerateMesh())
Draw(mesh)
print(mesh.GetBoundaries())
print(mesh.GetMaterials())
