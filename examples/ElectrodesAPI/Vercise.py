"""
A single vercise electrode
is generated.
It is meshed and the contacts
are colored to highlight their location.
"""

import ossdbs
from ngsolve import Draw, Mesh
import netgen.occ as occ

settings = {"Electrodes":
            [{"Name": "BostonScientificVercise",
              "Rotation[Degrees]": 0,
              "Direction": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 1},
              "TipPosition": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 0},
              },
             ]
            }

electrodes = ossdbs.generate_electrodes(settings)
vercise = electrodes[0]
occgeo = occ.OCCGeometry(vercise.geometry)
mesh = Mesh(occgeo.GenerateMesh())
Draw(mesh)
print(mesh.GetBoundaries())
print(mesh.GetMaterials())

bnd_dict = {}
for idx, contact in enumerate(vercise.boundaries):
    bnd_dict[contact] = idx

Draw(mesh.BoundaryCF(bnd_dict), mesh, "bnd")
