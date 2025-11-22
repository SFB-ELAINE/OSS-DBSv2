"""
Two vercise electrodes
with an encapsulation layer
is generated.
The encapsulation layers have
different thicknesses.
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
    "ExportElectrode": False,
}

electrodes = ossdbs.generate_electrodes(settings)
geo = None
for idx, vercise_settings in enumerate(settings["Electrodes"]):
    vercise = electrodes[idx]
    encap = vercise.encapsulation_geometry(
        vercise_settings["EncapsulationLayer"]["Thickness[mm]"]
    )
    if geo is None:
        geo = occ.Glue([vercise.geometry, encap])
    else:
        geo = occ.Glue([geo, vercise.geometry, encap])
occgeo = occ.OCCGeometry(geo)
with TaskManager():
    mesh = Mesh(occgeo.GenerateMesh())
Draw(mesh)
print(mesh.GetBoundaries())
print(mesh.GetMaterials())
# Clip in the center and show the scalar function to see the layer
Draw(mesh.MaterialCF({"EncapsulationLayer": 1.0}, default=0), mesh, "Material")
