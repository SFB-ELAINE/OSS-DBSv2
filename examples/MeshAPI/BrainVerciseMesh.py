"""
Example of a full model geometry
with one Vercise electrode without
encapsulation layer.
The mesh is generated and then
refined at one surface.
"""
import ossdbs
from ngsolve import Draw, Redraw

settings = \
    {"Electrodes":
        [{"Name": "BostonScientificVercise",
          "Rotation[Degrees]": 0,
          "Direction": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 1},
          "TipPosition": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 0},
          "EncapsulationLayer":
              {"Thickness[mm]": 0.0,  # indicates that no encapsulation is modelled
               },

          },
         ],
     "MaterialDistribution":
        {"MRIPath": "../BrainGeometryAPI/segmask.nii.gz"
         },
     "BrainRegion":
        {"Center": {"x[mm]": 5, "y[mm]": 14, "z[mm]": -4.5},
         "Dimension": {"x[mm]": 50.0, "y[mm]": 50.0, "z[mm]": 50.0},
         "Shape": "Ellipsoid"
         },
     "Mesh":
        {
         "LoadMesh": False,
         "SaveMesh": False
         }
     }

mesh = ossdbs.generate_mesh(settings)
print(mesh.boundaries)
print(mesh.materials)
print(mesh.ngsolvemesh.ne)
Draw(mesh.ngsolvemesh)
# Clip the drawing to observe the refinement
input("Hit enter to refine mesh")
mesh.refine_at_boundaries(["E1C1"])
print(mesh.ngsolvemesh.ne)
Redraw()

input("Hit enter to refine mesh")
mesh.refine_at_boundaries(["E1C1"])
print(mesh.ngsolvemesh.ne)
Redraw()

input("Hit to curve mesh")
mesh.curve(2)
Redraw()
