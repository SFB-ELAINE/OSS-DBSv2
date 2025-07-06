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
            "Name": "BostonScientificVercise",
            "Rotation[Degrees]": 0,
            "Direction": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 1},
            "TipPosition": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 0},
            "EncapsulationLayer": {
                "Thickness[mm]": 0.1,
            },
        },
    ],
    "ExportElectrode": False,
}

electrodes = ossdbs.generate_electrodes(settings)
vercise_settings = settings["Electrodes"][0]
vercise = electrodes[0]
encap = vercise.encapsulation_geometry(
    vercise_settings["EncapsulationLayer"]["Thickness[mm]"]
)
occgeo = occ.OCCGeometry(occ.Glue([vercise.geometry, encap]))
with TaskManager():
    mesh = Mesh(occgeo.GenerateMesh()).Curve(2)
Draw(mesh)
print(mesh.GetBoundaries())
print(mesh.GetMaterials())

# Clip in the center and show the scalar function to see the layer
Draw(mesh.MaterialCF({"EncapsulationLayer": 1.0}), mesh, "Material")
