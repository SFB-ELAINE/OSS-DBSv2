from ossdbs.point_analysis import imp_coord


def load_default_for_lead(settings):
    """Add parameters that are not defined in Lead-DBS GUI

    Parameters
    ----------

    settings: dict, OSS-DBS settings imported from Lead-DBS

    Returns
    -------

    settings: dict

    """

    settings["BrainRegion"]["Dimension"]["x[mm]"] = 60
    settings["BrainRegion"]["Dimension"]["y[mm]"] = 60
    settings["BrainRegion"]["Dimension"]["z[mm]"] = 80
    settings["BrainRegion"]["Shape"] = "Ellipsoid"
    settings["StimulationSignal"]["Type"] = "Multisine"
    settings["StimulationSignal"]["ListOfFrequencies"] = [10000]

    settings["PointModel"]["VoxelLattice"] = {"Active": False,
                                              "Shape": {'x': 31, 'y': 31, 'z': 41}}
    settings["Electrodes"][settings["ModelSide"]]["EncapsulationLayer"]["Thickness"] = 0.1
    settings["Electrodes"][settings["ModelSide"]]["EncapsulationLayer"]["Material"] = "White matter"
    settings["CalcAxonActivation"] = False
    settings["ExportVTK"] = True
    settings["Mesh"]["MeshingHypothesis"]["Type"] = "Fine"
    settings["FEMOrder"] = 2
    settings["ComputeImpedance"] = False
    settings["SaveImpedance"] = False

    first_contact = imp_coord(settings)
    settings["PointModel"]["Lattice"] = {"Active": True,
                                         "Center": {
                                             "x[mm]": first_contact[0],
                                             "y[mm]": first_contact[1] + 2.0,
                                             "z[mm]": first_contact[2] + 3.0
                                            },
                                         "Shape": {
                                             "x": 31,
                                             "y": 31,
                                             "z": 41
                                            },
                                         "Direction": {
                                             "x[mm]": settings["Electrodes"][0]["Direction"]["x[mm]"],
                                             "y[mm]": settings["Electrodes"][0]["Direction"]["y[mm]"],
                                             "z[mm]": settings["Electrodes"][0]["Direction"]["z[mm]"]
                                            },
                                         "PointDistance[mm]": 0.5
                                         }

    return settings
