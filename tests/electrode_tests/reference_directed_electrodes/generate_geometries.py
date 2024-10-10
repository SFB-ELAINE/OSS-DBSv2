"""Generate geometries for all directed electrodes.
Four directions are tested to check if the spin of the
electrode agrees with the Lead-DBS convention.
"""

import ngsolve
from tested_directions import base_settings, directions, get_direction_dict

import ossdbs

directed_electrodes = [
    "BostonScientificVerciseDirected",
    "AbbottStJudeDirected6172",
    "AbbottStJudeDirected6173",
    "MedtronicSenSightB33015",
    "MedtronicSenSightB33005",
    "BostonScientificCartesiaX",
    "BostonScientificCartesiaHX",
]


with ngsolve.TaskManager():
    for electrode in directed_electrodes:
        settings = base_settings.copy()
        settings["Electrodes"][0]["Name"] = electrode
        for idx, direction in enumerate(directions):
            settings["Electrodes"][0]["Direction"] = get_direction_dict(direction)
            mesh = ossdbs.generate_mesh(settings)
            mesh.save(f"mesh_{electrode}_direction_{idx}.vol.gz")
