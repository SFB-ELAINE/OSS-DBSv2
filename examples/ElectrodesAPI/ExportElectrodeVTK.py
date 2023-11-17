"""
A single electrode is generated and exported as VTK file.
It is meshed and the contacts are colored to highlight their location.
"""

import ossdbs
from ngsolve import Mesh, VTKOutput, BND
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
electrode = electrodes[0]
occgeo = occ.OCCGeometry(electrode.geometry)
mesh = Mesh(occgeo.GenerateMesh())

bnd_dict = {}
for idx, contact in enumerate(electrode.boundaries):
    bnd_dict[contact] = idx

boundary_cf = mesh.BoundaryCF(bnd_dict, default=-1)

VTKOutput(
  ma=mesh,
  coefs=[boundary_cf],
  names=["boundaries"],
  filename="electrode",
  subdivision=0,
).Do(vb=BND)
