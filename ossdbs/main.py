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
    build_brain_model,
    generate_electrodes,
    generate_signal,
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
from ossdbs.model_geometry import ModelGeometry
from ossdbs.utils.settings import Settings
from ossdbs.utils.type_check import TypeChecker

_logger = logging.getLogger(__name__)


def main_run(input_settings: dict):
    """Run OSS-DBS from input dictionary.

    Parameters
    ----------
    input_settings: dict
        Input dictionary
    run_path: str
        Path where to run OSS-DBS
    """
    timings = {}
    time_0 = time.time()

    settings = Settings(input_settings).complete_settings()
    TypeChecker.check(settings)

    _logger.debug(f"Final settings:\\ {settings}")

    # create output path
    if not os.path.isdir(settings["OutputPath"]):
        os.mkdir(settings["OutputPath"])
    log_to_file(
        output_file=os.path.join(settings["OutputPath"], "ossdbs.log"),
        level=_logger.getEffectiveLevel(),
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
    brain_model = build_brain_model(settings, mri_image)
    try:
        geometry = ModelGeometry(brain_model, electrodes)
    except RuntimeError:
        _logger.warning(
            "Initial geometry failed, now building with rotated geometry."
            "If this fails, too, change the shape of the brain geometry."
        )
        brain_model = build_brain_model(settings, mri_image, rotate_initial_geo=True)
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
        brain_model.brain_region,
        dielectric_properties,
        materials,
        geometry.encapsulation_layers,
        complex_data=settings["EQSMode"],
        dti_image=dti_image,
    )

    time_1 = time.time()
    timings["ConductivityCF"] = time_1 - time_0
    time_0 = time_1

    # decide on truncation
    truncation_time = None
    if "TruncateAfterActivePartRatio" in settings:
        truncation_ratio = settings["TruncateAfterActivePartRatio"]
        if truncation_ratio is not None:
            if not isinstance(truncation_ratio, float):
                raise ValueError(
                    "Please provide the ratio to truncate the signal "
                    "as a floating-point number. "
                    "Set e.g. to 20 for 20 times pulse + counterpulse width."
                )
            time_domain_signal = generate_signal(settings)
            truncation_time = truncation_ratio * time_domain_signal.get_active_time()

    # save Mesh for StimSets
    if settings["StimSets"]["Active"]:
        settings["Mesh"]["SaveMesh"] = True
        settings["Mesh"]["SavePath"] = os.path.join(settings["OutputPath"], "tmp_mesh")
        settings["Mesh"]["LoadPath"] = os.path.join(
            settings["OutputPath"], "tmp_mesh.vol.gz"
        )
        settings["Mesh"]["LoadMesh"] = False
        # because of floating
        settings["Solver"]["Preconditioner"] = "local"
        settings["Solver"]["PreconditionerKwargs"] = {}
    # run in parallel
    with ngsolve.TaskManager():
        solver = prepare_solver(settings)
        volume_conductor = prepare_volume_conductor_model(
            settings, geometry, conductivity, solver
        )
        frequency_domain_signal = prepare_stimulation_signal(settings)
        if not settings["StimSets"]["Active"]:
            vcm_timings = run_volume_conductor_model(
                settings,
                volume_conductor,
                frequency_domain_signal,
                truncation_time=truncation_time,
            )
            _logger.info(f"Volume conductor timings:\n{pprint.pformat(vcm_timings)}")
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

    set_logger(level=args.loglevel)
    _logger.info("Loading settings from input file")
    _logger.debug(f"Input file: {args.input_dictionary}")
    with open(args.input_dictionary) as json_file:
        input_settings = json.load(json_file)

    # add the stimulation folder (where input dict.json is stored, needed for Lead-DBS)
    input_settings["StimulationFolder"] = os.path.dirname(
        os.path.abspath(args.input_dictionary)
    )

    main_run(input_settings)


if __name__ == "__main__":
    main()
