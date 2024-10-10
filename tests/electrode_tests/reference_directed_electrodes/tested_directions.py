"""All tested directions."""

directions = [(-0.5775, -0.5775, -0.5775), (0.7075, 0.0, 0.7075), (1, 1, 1), (0, 0, 1)]

base_settings = {
    "Electrodes": [
        {
            "Rotation[Degrees]": 0,
            "TipPosition": {"x[mm]": 0, "y[mm]": 15, "z[mm]": -3},
        }
    ],
    "ExportElectrode": False,
    "OutputPath": ".",
    "BrainRegion": {
        "Center": {"x[mm]": 5, "y[mm]": 14, "z[mm]": -4.5},
        "Dimension": {"x[mm]": 50.0, "y[mm]": 50.0, "z[mm]": 50.0},
        "Shape": "Ellipsoid",
    },
    "Mesh": {
        "LoadMesh": False,
        "SaveMesh": False,
        "MeshingHypothesis": {"Type": "Coarse"},
    },
}


def get_direction_dict(direction: tuple) -> dict:
    return {"x[mm]": direction[0], "y[mm]": direction[1], "z[mm]": direction[2]}
