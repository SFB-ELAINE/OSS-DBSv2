"""
Example of a full model geometry
with one Vercise electrode without
encapsulation layer.
The mesh is generated and
local mesh refinements are taken
into account during the meshing
process.
"""

import ossdbs

settings = {
    "Electrodes": [
        {
            "Name": "BostonScientificVerciseDirected",
            "Rotation[Degrees]": 0,
            "Direction": {"x[mm]": 0, "y[mm]": 1, "z[mm]": -1},
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
                    "MaxMeshSizeEdge": 0.01,
                },
                {
                    "Contact_ID": 3,
                    "Active": True,
                    "Current[A]": 0.0,
                    "Voltage[V]": 0.0,
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
                    "MaxMeshSizeEdge": 0.01,
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
    "BrainRegion": {
        "Center": {"x[mm]": 5, "y[mm]": 14, "z[mm]": -4.5},
        "Dimension": {"x[mm]": 50.0, "y[mm]": 50.0, "z[mm]": 50.0},
        "Shape": "Ellipsoid",
    },
    "ExportElectrode": False,
}

electrodes = ossdbs.generate_electrodes(settings)
vercise = electrodes[0]
electrode_settings = settings["Electrodes"][0]


contact2 = []
# set edge mesh sizes


for edge in vercise.geometry.edges:
    if edge.name is not None:
        if edge.name in dict:
            dict[edge.name] += 1
        else:
            dict[edge.name] = 1

        if edge.name == "Contact_2":
            contact2.append((edge.center.x, edge.center.y, edge.center.z))

        contact_number = int(edge.name.replace("Contact_", ""))
        for contact_setting in electrode_settings["Contacts"]:
            if contact_setting["Contact_ID"] == contact_number:
                if "MaxMeshSizeEdge" in contact_setting:
                    edge.maxh = contact_setting["MaxMeshSizeEdge"]
                continue

vercise.export_electrode(
    output_path=".", brain_dict=settings["BrainRegion"], n_electrode=0
)
