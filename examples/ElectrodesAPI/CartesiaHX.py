"""
A single vercise electrode
with an encapsulation layer is generated.
It is meshed and the contacts
are colored to highlight their location.
"""

import netgen.occ as occ
from ngsolve import Draw, Mesh, TaskManager

import ossdbs

settings = {
    "Electrodes": [
        {
            "Name": "BostonScientificCartesiaHX",
            "Rotation[Degrees]": 0,
            "Direction": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 1},
            "TipPosition": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 0},
            "EncapsulationLayer": {
                "Thickness[mm]": 0.1,
            },
        },
    ],
    "BrainRegion": {
        "Center": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 0},
        "Dimension": {"x[mm]": 50.0, "y[mm]": 50.0, "z[mm]": 50.0},
    },
    "OutputPath": "./",
    "ExportElectrode": True,
}

electrodes = ossdbs.generate_electrodes(settings)
vercise_settings = settings["Electrodes"][0]
vercise = electrodes[0]
vercise.geometry.WriteStep("CartesiaHX.stp")
occgeo = occ.OCCGeometry(vercise.geometry)
with TaskManager():
    mesh = Mesh(occgeo.GenerateMesh()).Curve(2)
Draw(mesh)
print(mesh.GetBoundaries())
print(mesh.GetMaterials())


bnd_dict = {}
for idx, contact in enumerate(vercise.boundaries):
    bnd_dict[contact] = idx

Draw(mesh.BoundaryCF(bnd_dict), mesh, "bnd")
