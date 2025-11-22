import ossdbs

settings = {
    "MaterialDistribution": {
        "MRIPath": "",
        "MRIMapping": {
            "Unknown": 0,
            "CSF": 1,
            "White matter": 2,
            "Gray matter": 3,
            "Blood": 4,
        },
        "DiffusionTensorActive": False,
        "DTIPath": "",
    },
    "DielectricModel": {
        "Type": "Constant",
        "CustomParameters": {
            "Gray matter": {"permittivity": 2.22e4, "conductivity": 1.15e-1},
            "Unknown": {"permittivity": 2.22e4, "conductivity": 1.15e-1},
            "White matter": {"permittivity": 1.25e4, "conductivity": 6.95e-2},
            "CSF": {"permittivity": 1.09e2, "conductivity": 2.0},
            "Blood": {"permittivity": 5.25e3, "conductivity": 7e-1},
        },
    },
}

dielectric_properties = ossdbs.prepare_dielectric_properties(settings)

for material in settings["MaterialDistribution"]["MRIMapping"]:
    print(
        f"Constant complex permittivity of {material} "
        f"is {dielectric_properties[material].complex_permittivity(130)}"
    )
