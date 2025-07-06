"""
Example of a full model geometry
with one Vercise electrode without
encapsulation layer.
"""

from pprint import pprint

from ngsolve import Draw, Mesh, TaskManager

import ossdbs

settings = {
    "Electrodes": [
        {
            "Name": "BostonScientificVercise",
            "Rotation[Degrees]": 0,
            "Direction": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 1},
            "TipPosition": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 0},
            "EncapsulationLayer": {
                "Thickness[mm]": 0.5,
            },
        },
        {
            "Name": "BostonScientificVercise",
            "Rotation[Degrees]": 0,
            "Direction": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 1},
            "TipPosition": {"x[mm]": 10, "y[mm]": 0, "z[mm]": 0},
            "EncapsulationLayer": {
                "Thickness[mm]": 1.0,
            },
        },
    ],
    "MaterialDistribution": {"MRIPath": "segmask.nii.gz"},
    "BrainRegion": {
        "Center": {"x[mm]": 5, "y[mm]": 14, "z[mm]": -4.5},
        "Dimension": {"x[mm]": 50.0, "y[mm]": 50.0, "z[mm]": 50.0},
        "Shape": "Ellipsoid",
    },
}

model_geometry = ossdbs.generate_model_geometry(settings)

occgeo = model_geometry.geometry
with TaskManager():
    mesh = Mesh(occgeo.GenerateMesh())
Draw(mesh)
bnd_dict = {}
for idx, contact in enumerate(model_geometry.contacts):
    bnd_dict[contact.name] = idx
bnd_dict["Body"] = -1
bnd_dict["EncapsulationLayerSurface_1"] = -2
bnd_dict["EncapsulationLayerSurface_2"] = -3
Draw(mesh.BoundaryCF(bnd_dict), mesh, "BND")
print(mesh.GetBoundaries())
print(mesh.GetMaterials())
# List of contacts used later to impose BCs
pprint(model_geometry.contacts)
