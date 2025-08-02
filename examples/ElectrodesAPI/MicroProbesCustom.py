"""
Generates MicroProbes electrode for rodents.
Custom parameters are used.
"""

import netgen.occ as occ
from ngsolve import Draw, Mesh, TaskManager

import ossdbs

parameters = {
    "exposed_wire": 0.05,
    "contact_radius": 0.1125,
    "lead_radius": 0.1125,
    "total_length": 13.3,
    "wire_radius": 0.09,
}

settings = {
    "Electrodes": [
        {
            "Name": "MicroProbesRodentElectrodeCustom",
            "CustomParameters": parameters,
            "Rotation[Degrees]": 0,
            "Direction": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 1},
            "TipPosition": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 0},
        }
    ],
    "ExportElectrode": False,
}

electrodes = ossdbs.generate_electrodes(settings)
electrode = electrodes[0]
occgeo = occ.OCCGeometry(electrode.geometry)
with TaskManager():
    mesh = Mesh(occgeo.GenerateMesh()).Curve(2)
Draw(mesh)
print(mesh.GetBoundaries())

bnd_dict = {}
for idx, contact in enumerate(electrode.boundaries):
    bnd_dict[contact] = idx

Draw(mesh.BoundaryCF(bnd_dict), mesh, "bnd")
