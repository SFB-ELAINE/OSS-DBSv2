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
                "Thickness[mm]": 0.0,  # indicates that no encapsulation is modelled
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

brain = ossdbs.generate_brain_model(settings)
electrodes = ossdbs.generate_electrodes(settings)
model_geometry = ossdbs.ModelGeometry(brain, electrodes)

occgeo = model_geometry.geometry
with TaskManager():
    mesh = Mesh(occgeo.GenerateMesh())
Draw(mesh)
print(mesh.GetBoundaries())
print(mesh.GetMaterials())
# List of contacts used later to impose BCs
pprint(model_geometry.contacts)
