"""
Example of a full model geometry
with one Vercise electrode with
encapsulation layer.
The mesh is generated and
local mesh refinements are taken
into account during the meshing
process.
"""

from ngsolve import Draw

import ossdbs

ossdbs.set_logger()

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
            "TipPosition": {"x[mm]": -9.48, "y[mm]": 11.61, "z[mm]": 4.68},
            "EncapsulationLayer": {"Thickness[mm]": 1.0, "MaxMeshSize": 0.2},
            "Contacts": [
                {
                    "Contact_ID": 1,
                    "Active": True,
                    "Voltage[V]": 1.0,
                    "MaxMeshSizeEdge": 0.05,
                },
                {
                    "Contact_ID": 3,
                    "Active": True,
                    "Voltage[V]": 0.0,
                    "MaxMeshSizeEdge": 0.05,
                },
            ],
        },
    ],
    "MaterialDistribution": {"MRIPath": "../../input_files/Butenko_segmask.nii.gz"},
    "Mesh": {"LoadMesh": False, "SaveMesh": False},
    "ExportElectrode": False,
}

mesh = ossdbs.generate_mesh(settings)
print(mesh.boundaries)
print(mesh.materials)
print(mesh.ngsolvemesh.ne)

Draw(mesh.ngsolvemesh)

bnd_dict = {}
for idx in range(1, 9):
    bnd_dict[f"E1C{idx}"] = idx

Draw(mesh.boundary_coefficients(bnd_dict), mesh.ngsolvemesh, "bnd")
