"""
Example of a full model geometry
with one Vercise electrode without
encapsulation layer.
The mesh is generated and then
refined at one surface.
"""

from ngsolve import Draw

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
            "TipPosition": {"x[mm]": -9.48, "y[mm]": 11.61, "z[mm]": 4.68},
            "EncapsulationLayer": {"Thickness[mm]": 0.0},
            "Contacts": [
                {
                    "Contact_ID": 1,
                    "Active": True,
                    "Voltage[V]": 1.0,
                },
                {
                    "Contact_ID": 3,
                    "Active": True,
                    "Voltage[V]": 0.0,
                },
            ],
        },
    ],
    "Mesh": {
        "LoadMesh": False,
        "SaveMesh": False,
        "AdaptiveMeshRefinement": {
            "Active": False,
            "MaxIterations": 2,
            "ErrorTolerance": 0.1,
        },
        "HPRefinement": {
            "Active": True,
            "Levels": 2,
            "Factor": 0.125,
        },
    },
    "ExportElectrode": False,
}

hp_mesh = ossdbs.generate_mesh(settings)
Draw(hp_mesh.ngsolvemesh)

settings["Mesh"]["HPRefinement"]["Active"] = False
mesh = ossdbs.generate_mesh(settings)

print("Number of elements before refinement:", mesh.ngsolvemesh.ne)
print("Number of elements after refinement:", hp_mesh.ngsolvemesh.ne)
