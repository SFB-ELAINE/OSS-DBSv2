"""
Example of a full model geometry
with one Vercise electrode without
encapsulation layer.
The mesh is generated and
local mesh refinements are taken
into account during the meshing
process.
"""
from ngsolve import Draw

import ossdbs

settings = {
    "Electrodes": [
        {
            "Name": "BostonScientificVerciseDirected",
            "Rotation[Degrees]": 0,
            "Direction": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 1},
            "TipPosition": {"x[mm]": 0, "y[mm]": 15, "z[mm]": -3},
            "EncapsulationLayer": {
                "Thickness[mm]": 0.0,  # indicates that no encapsulation is modelled
            },
            "Contacts": [
                {
                    "Contact_ID": 1,
                    "Active": True,
                    "Current[A]": 0.0,
                    "Voltage[V]": 1.0,
                    "Floating": False,
                    "SurfaceImpedance[Ohmm]": {"real": 0.0, "imag": 0.0},
                    "MaxMeshSizeEdge": 0.05,
                },
                {
                    "Contact_ID": 2,
                    "MaxMeshSizeEdge": 0.05,
                },
                {
                    "Contact_ID": 3,
                    "Active": True,
                    "Current[A]": 0.0,
                    "Voltage[V]": 0.0,
                    "Floating": False,
                    "SurfaceImpedance[Ohmm]": {"real": 0.0, "imag": 0.0},
                    "MaxMeshSizeEdge": 0.05,
                },
                {
                    "Contact_ID": 4,
                    "MaxMeshSizeEdge": 0.05,
                },
                {
                    "Contact_ID": 5,
                    "Active": False,
                    "Current[A]": 0.0,
                    "Voltage[V]": 0.0,
                    "Floating": False,
                    "SurfaceImpedance[Ohmm]": {"real": 0.0, "imag": 0.0},
                    "MaxMeshSize": 0.25,
                },
                {
                    "Contact_ID": 6,
                    "MaxMeshSizeEdge": 0.05,
                },
                {
                    "Contact_ID": 7,
                    "MaxMeshSizeEdge": 0.05,
                },
                {
                    "Contact_ID": 8,
                    "MaxMeshSizeEdge": 0.05,
                },
            ],
        }
    ],
    "MaterialDistribution": {"MRIPath": "../BrainGeometryAPI/segmask.nii.gz"},
    "BrainRegion": {
        "Center": {"x[mm]": 5, "y[mm]": 14, "z[mm]": -4.5},
        "Dimension": {"x[mm]": 50.0, "y[mm]": 50.0, "z[mm]": 50.0},
        "Shape": "Ellipsoid",
    },
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
