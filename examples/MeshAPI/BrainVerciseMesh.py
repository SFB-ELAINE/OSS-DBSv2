"""
Example of a full model geometry
with one Vercise electrode without
encapsulation layer.
The mesh is generated and then
refined at one surface.
"""

from ngsolve import Draw, Redraw

import ossdbs

settings = {
    "BrainRegion": {
        "Center": {"x[mm]": -9.48, "y[mm]": 11.61, "z[mm]": 4.68},
        "Dimension": {"x[mm]": 40.0, "y[mm]": 40.0, "z[mm]": 40.0},
        "Shape": "Ellipsoid",
    },
    "Electrodes": [
        {
            "Name": "BostonScientificVercise",
            "Rotation[Degrees]": 0,
            "Direction": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 1},
            "TipPosition": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 0},
            "EncapsulationLayer": {"Thickness[mm]": 0.0},
        },
    ],
    "MaterialDistribution": {"MRIPath": "../../input_files/Butenko_segmask.nii.gz"},
    "Mesh": {
        "LoadMesh": False,
        "SaveMesh": False,
    },
    "ExportElectrode": False,
}

mesh = ossdbs.generate_mesh(settings)
print(mesh.boundaries)
print(mesh.materials)
print("Number of elements:", mesh.ngsolvemesh.ne)
Draw(mesh.ngsolvemesh)
# Clip the drawing to observe the refinement
input("Hit enter to refine mesh")
mesh.refine_at_boundaries(["E1C1"])
print("Number of elements:", mesh.ngsolvemesh.ne)
Redraw()

input("Hit enter to refine mesh")
mesh.refine_at_boundaries(["E1C1"])
print("Number of elements:", mesh.ngsolvemesh.ne)
Redraw()

input("Hit to curve mesh")
mesh.curve(2)
Redraw()
