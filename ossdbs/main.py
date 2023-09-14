import os
import json
import ngsolve
from ossdbs import set_logger
from ossdbs.api import (load_images,
                        generate_electrodes,
                        prepare_dielectric_properties,
                        create_bounding_box,
                        prepare_solver,
                        prepare_volume_conductor_model,
                        run_volume_conductor_model,
                        set_contact_and_encapsulation_layer_properties,
                        )
from ossdbs.utils.settings import Settings
from ossdbs.utils.type_check import TypeChecker
from ossdbs.model_geometry import (ModelGeometry,
                                   BrainGeometry,
                                   BoundingBox)
from ossdbs.fem import ConductivityCF
import logging
import time
import pprint
import argparse

_logger = logging.getLogger(__name__)


def main() -> None:

    parser = argparse.ArgumentParser(prog="OSS-DBS",
                                     description="Welcome to OSS-DBS v2.",
                                     epilog="Please report bugs and errors on GitHub")
    parser.add_argument('--loglevel', type=int, help="specify verbosity of logger",
                        default=logging.INFO)
    parser.add_argument('input_dictionary', type=str,
                        help="input dictionary in JSON format"
                        )
    args = parser.parse_args()

    timings = {}
    time_0 = time.time()
    set_logger(level=args.loglevel)

    # navigate to the stim folder (where input dict.json is stored)
    stim_folder, _ = os.path.split(args.input_dictionary)
    os.chdir(stim_folder)

    _logger.info("Loading settings from input file")
    _logger.debug("Input file: {}".format(args.input_dictionary))
    with open(args.input_dictionary, 'r') as json_file:
        input_settings = json.load(json_file)

    settings = Settings(input_settings).complete_settings()
    TypeChecker.check(settings)
    _logger.debug("Final settings:\\ {}".format(settings))

    # create output path
    if not os.path.isdir(settings["OutputPath"]):
        os.mkdir(settings["OutputPath"])
    # create fail flag
    open("fail_" + settings["FailFlag"] + ".txt", 'w').close()

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
        region_parameters = settings['BrainRegion']
        brain_region = create_bounding_box(region_parameters)
        shape = settings['BrainRegion']['Shape']
        brain_model = BrainGeometry(shape, brain_region)
    else:
        _logger.debug("Generating model geometry from MRI image")
        # attention: bounding box is given in voxel space!
        brain_region = mri_image.bounding_box
        shape = "Ellipsoid"
        # transformation to real space in geometry creation
        _logger.debug("Generate OCC model, passing transformation matrix from MRI image")
        brain_model = BrainGeometry(shape, brain_region, trafo_matrix=mri_image.trafo_matrix, translation=mri_image.translation)
        start, end = brain_model.geometry.bounding_box
        brain_region = BoundingBox(start, end)
        _logger.debug("Bounding box in real space: {}, {}". format(brain_region.start, brain_region.end))

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
    conductivity = ConductivityCF(mri_image,
                                  brain_region,
                                  dielectric_properties,
                                  materials,
                                  geometry.encapsulation_layers,
                                  complex_data=settings["EQSMode"],
                                  dti_image=dti_image
                                  )

    time_1 = time.time()
    timings["ConductivityCF"] = time_1 - time_0
    time_0 = time_1

    # run in parallel
    with ngsolve.TaskManager():
        solver = prepare_solver(settings)
        volume_conductor = prepare_volume_conductor_model(settings, geometry, conductivity, solver)
        vcm_timings = run_volume_conductor_model(settings, volume_conductor)

    time_1 = time.time()
    timings["VolumeConductor"] = time_1 - time_0

    _logger.info("Timings:\n {}".format(pprint.pformat(timings)))
    _logger.info("Volume conductor timings:\n {}".format(pprint.pformat(vcm_timings)))

    # write success file
    open("success_" + settings["FailFlag"] + ".txt", 'w').close()
    os.remove("fail_" + settings["FailFlag"] + ".txt")
    _logger.info("Process Completed")


if __name__ == '__main__':
    main()
