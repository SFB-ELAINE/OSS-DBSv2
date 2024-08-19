# Copyright 2023, 2024 Konstantin Butenko, Jan Philipp Payonk
# Copyright 2023, 2024 Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import argparse
import json
import logging
import os
import pprint
import time

import ngsolve

from ossdbs import log_to_file, set_logger
from ossdbs.api import (
    create_bounding_box,
    generate_electrodes,
    load_images,
    prepare_dielectric_properties,
    prepare_solver,
    prepare_stimulation_signal,
    prepare_volume_conductor_model,
    run_stim_sets,
    run_volume_conductor_model,
    set_contact_and_encapsulation_layer_properties,
)
from ossdbs.fem import ConductivityCF
from ossdbs.model_geometry import BoundingBox, BrainGeometry, ModelGeometry
from ossdbs.utils.settings import Settings
from ossdbs.utils.type_check import TypeChecker

_logger = logging.getLogger(__name__)


def main() -> None:
    """Main function to run OSS-DBS in CLI mode."""
    parser = argparse.ArgumentParser(
        prog="OSS-DBS",
        description="Welcome to OSS-DBS v2.",
        epilog="Please report bugs and errors on GitHub",
    )
    parser.add_argument(
        "--loglevel", type=int, help="specify verbosity of logger", default=logging.INFO
    )
    parser.add_argument(
        "input_dictionary", type=str, help="input dictionary in JSON format"
    )
    args = parser.parse_args()

    timings = {}
    time_0 = time.time()
    set_logger(level=args.loglevel)

    _logger.info("Loading settings from input file")
    _logger.debug(f"Input file: {args.input_dictionary}")
    with open(args.input_dictionary) as json_file:
        input_settings = json.load(json_file)

    settings = Settings(input_settings).complete_settings()
    TypeChecker.check(settings)

    # add the stimulation folder (where input dict.json is stored, needed for Lead-DBS)
    settings["StimulationFolder"] = os.path.dirname(
        os.path.abspath(args.input_dictionary)
    )

    _logger.debug(f"Final settings:\\ {settings}")

    # create output path
    if not os.path.isdir(settings["OutputPath"]):
        os.mkdir(settings["OutputPath"])
    log_to_file(
        output_file=os.path.join(settings["OutputPath"], "ossdbs.log"),
        level=args.loglevel,
    )
    # create fail flag
    open(
        os.path.join(
            settings["StimulationFolder"], "fail_" + settings["FailFlag"] + ".txt"
        ),
        "w",
    ).close()

    time_1 = time.time()
    timings["Settings"] = time_1 - time_0
    time_0 = time_1

    mri_image, dti_image = load_images(settings)

    time_1 = time.time()
    timings["MRI"] = time_1 - time_0
    time_0 = time_1

    electrodes = generate_electrodes(settings)

    time_1 = time.time()
    timings["Electrodes"] = time_1 - time_0
    time_0 = time_1

    _logger.info("Generate full model geometry")
    # MRI image is default choice for brain construction
    if "BrainRegion" in settings:
        _logger.debug("Generating model geometry for fixed brain region")
        region_parameters = settings["BrainRegion"]
        brain_region = create_bounding_box(region_parameters)
        shape = settings["BrainRegion"]["Shape"]
        brain_model = BrainGeometry(shape, brain_region)
    else:
        _logger.debug("Generating model geometry from MRI image")
        # attention: bounding box is given in voxel space!
        brain_region = mri_image.bounding_box
        shape = "Ellipsoid"
        # transformation to real space in geometry creation
        _logger.debug(
            "Generate OCC model, passing transformation matrix from MRI image"
        )
        brain_model = BrainGeometry(
            shape,
            brain_region,
            trafo_matrix=mri_image.trafo_matrix,
            translation=mri_image.translation,
        )
        start, end = brain_model.geometry.bounding_box
        brain_region = BoundingBox(start, end)
        _logger.debug(
            f"Bounding box in real space: {brain_region.start}, {brain_region.end}"
        )

    geometry = ModelGeometry(brain_model, electrodes)

    time_1 = time.time()
    timings["ModelGeometry"] = time_1 - time_0
    time_0 = time_1

    set_contact_and_encapsulation_layer_properties(settings, geometry)

    time_1 = time.time()
    timings["ContactProperties"] = time_1 - time_0
    time_0 = time_1

    dielectric_properties = prepare_dielectric_properties(settings)

    time_1 = time.time()
    timings["DielectricModel"] = time_1 - time_0
    time_0 = time_1

    _logger.info("Prepare conductivity coefficient function")
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

    time_1 = time.time()
    timings["ConductivityCF"] = time_1 - time_0
    time_0 = time_1

    # save Mesh for StimSets
    if settings["StimSets"]["Active"]:
        settings["Mesh"]["SaveMesh"] = True
        settings["Mesh"]["SavePath"] = os.path.join(settings["OutputPath"], "tmp_mesh")
        settings["Mesh"]["LoadPath"] = os.path.join(
            settings["OutputPath"], "tmp_mesh.vol.gz"
        )
        settings["Mesh"]["LoadMesh"] = False
        # because of floating
        settings["Solver"]["PreconditionerKwargs"] = {"coarsetype": "local"}
    # run in parallel
    with ngsolve.TaskManager():
        solver = prepare_solver(settings)
        volume_conductor = prepare_volume_conductor_model(
            settings, geometry, conductivity, solver
        )
        frequency_domain_signal = prepare_stimulation_signal(settings)
        if not settings["StimSets"]["Active"]:
            vcm_timings = run_volume_conductor_model(
                settings, volume_conductor, frequency_domain_signal
            )
            _logger.info("Volume conductor timings:\n" f"{pprint.pformat(vcm_timings)}")
        else:
            # mesh was saved already
            settings["Mesh"]["SaveMesh"] = False
            settings["Mesh"]["LoadMesh"] = True
            run_stim_sets(
                settings, geometry, conductivity, solver, frequency_domain_signal
            )

    time_1 = time.time()
    timings["VolumeConductor"] = time_1 - time_0
    time_0 = time_1

    # run PAM
    if settings["PathwayFile"] is not None:
        _logger.info("Please compute the pathway activation separately.")
        # commented because of interaction with Lead-DBS
        """
        if settings["StimSets"]["Active"]:
            _logger.info(
                "No PAM run because you specified StimSets."
                "Compute the pathway activation separately."
            )
        elif settings["CalcAxonActivation"] is False:
            _logger.info("Axon activation is not computed.")
        else:
            run_PAM(settings)
            time_1 = time.time()
            timings["PAM"] = time_1 - time_0
        """

    _logger.info(f"Timings:\n {pprint.pformat(timings)}")

    # write success file
    open(
        os.path.join(
            settings["StimulationFolder"], "success_" + settings["FailFlag"] + ".txt"
        ),
        "w",
    ).close()
    os.remove(
        os.path.join(
            settings["StimulationFolder"], "fail_" + settings["FailFlag"] + ".txt"
        )
    )
    _logger.info("Process Completed")


if __name__ == "__main__":
    main()
