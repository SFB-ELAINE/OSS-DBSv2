"""
Example of a full model geometry
with one Vercise electrode without
encapsulation layer.
The mesh is generated and adaptive
refinement is used.
"""

import logging

import ngsolve
from ngsolve import Draw

import ossdbs
from ossdbs.api import (
    create_bounding_box,
    generate_electrodes,
    load_images,
    prepare_dielectric_properties,
    prepare_solver,
    prepare_stimulation_signal,
    prepare_volume_conductor_model,
    run_volume_conductor_model,
    set_contact_and_encapsulation_layer_properties,
)
from ossdbs.fem import ConductivityCF
from ossdbs.model_geometry import BrainGeometry, ModelGeometry
from ossdbs.utils.settings import Settings

ossdbs.set_logger(logging.INFO)

input_settings = {
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
    "MaterialDistribution": {"MRIPath": "../../input_files/Butenko_segmask.nii.gz"},
    "Mesh": {
        "AdaptiveMeshRefinement": {
            "Active": True,
            "MaxIterations": 2,
            "ErrorTolerance": 0.1,
        },
        "HPRefinement": {
            "Active": False,
            "Levels": 2,
            "Factor": 0.125,
        },
    },
    "StimulationSignal": {
        "Type": "Multisine",
        "ListOfFrequencies": [1000.0],
    },
    "ComputeImpedance": True,
}

settings = Settings(input_settings).complete_settings()

mri_image, dti_image = load_images(settings)

electrodes = generate_electrodes(settings)

brain_region = create_bounding_box(settings["BrainRegion"])
brain_model = BrainGeometry(settings["BrainRegion"]["Shape"], brain_region)

geometry = ModelGeometry(brain_model, electrodes)

set_contact_and_encapsulation_layer_properties(settings, geometry)
dielectric_properties = prepare_dielectric_properties(settings)
materials = settings["MaterialDistribution"]["MRIMapping"]

conductivity = ConductivityCF(
    mri_image,
    brain_region,
    dielectric_properties,
    materials,
    geometry.encapsulation_layers,
    complex_data=settings["EQSMode"],
    dti_image=dti_image,
)

mesh = ossdbs.generate_mesh(settings)

with ngsolve.TaskManager():
    solver = prepare_solver(settings)
    volume_conductor = prepare_volume_conductor_model(
        settings, geometry, conductivity, solver
    )
    frequency_domain_signal = prepare_stimulation_signal(settings)

    run_volume_conductor_model(settings, volume_conductor, frequency_domain_signal)

print("Number of elements before refinement:", mesh.ngsolvemesh.ne)
print("Number of elements after refinement:", volume_conductor.mesh.ngsolvemesh.ne)
print("Number degrees of freedom:", ngsolve.H1(volume_conductor.mesh.ngsolvemesh).ndof)
Draw(volume_conductor.mesh.ngsolvemesh)
