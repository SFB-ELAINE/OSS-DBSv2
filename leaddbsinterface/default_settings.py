# Copyright 2023, 2024 Konstantin Butenko, Jan Philipp Payonk, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later


def load_default_for_lead(settings):
    """Add parameters that are not defined in Lead-DBS GUI.

    Parameters
    ----------
    settings: dict
        OSS-DBS settings imported from Lead-DBS

    Returns
    -------
    settings: dict
        OSS-DBS compatible settings

    """
    settings["BrainRegion"]["Dimension"]["x[mm]"] = 60
    settings["BrainRegion"]["Dimension"]["y[mm]"] = 60
    settings["BrainRegion"]["Dimension"]["z[mm]"] = 80
    settings["BrainRegion"]["Shape"] = "Ellipsoid"
    settings["StimulationSignal"]["Type"] = "Multisine"
    settings["StimulationSignal"]["ListOfFrequencies"] = [10000]

    settings["PointModel"]["VoxelLattice"] = {
        "Active": False,
        "Shape": {"x": 31, "y": 31, "z": 41},
    }
    settings["Electrodes"][settings["ModelSide"]]["EncapsulationLayer"][
        "Thickness[mm]"
    ] = 0.1
    settings["Electrodes"][settings["ModelSide"]]["EncapsulationLayer"]["Material"] = (
        "White matter"
    )
    settings["ExportVTK"] = True
    settings["Mesh"]["MeshingHypothesis"]["Type"] = "Default"
    # material refinement
    settings["Mesh"]["MaterialRefinementSteps"] = 1
    # edge refinement is defined in other file depending on lead geometry
    settings["FEMOrder"] = 2
    settings["ComputeImpedance"] = False

    settings["Solver"]["MaximumSteps"] = 500
    settings["Solver"]["Precision"] = 1e-10

    return settings


def initialize_default_settings():
    """Initialize with necessary sub-dictionaries.

    Returns
    -------
    settings_template: dict

    """
    settings_template = {
        "ModelSide": 0,  # hardwired
        "BrainRegion": {
            "Dimension": {},
            "Shape": {},
        },
        "StimulationSignal": {},
        "PointModel": {
            "VoxelLattice": {
                "Shape": {},
            }
        },
        "Electrodes": [{"EncapsulationLayer": {}}],
        "Mesh": {"MeshingHypothesis": {}},
        "Solver": {},
    }

    return settings_template


def update_default_dict(default_settings: dict, custom_settings: dict) -> None:
    """Update default settings with custom input.

    Parameters
    ----------
    default_settings: dict
        predefined settings for some parameters
    custom_settings: dict
        provided settings, e.g. from Lead-DBS

    Returns
    -------
    updated_settings: dict

    """
    for key in custom_settings.keys():
        is_dict = False

        if key in default_settings.keys():
            if isinstance(default_settings[key], dict):
                # empty dicts yield False
                is_dict = bool(default_settings[key])

        if is_dict:
            update_default_dict(default_settings[key], custom_settings[key])
        else:
            default_settings[key] = custom_settings[key]

    updated_settings = default_settings
    return updated_settings
