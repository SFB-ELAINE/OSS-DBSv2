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
                        generate_neuron_grid
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
import pandas as pd
import numpy as np
import h5py
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

    _logger.info("Loading settings from input file")
    _logger.debug("Input file: {}".format(args.input_dictionary))
    with open(args.input_dictionary, 'r') as json_file:
        input_settings = json.load(json_file)

    settings = Settings(input_settings).complete_settings()
    TypeChecker.check(settings)

    # create fail flag
    open(os.path.join(settings["OutputPath"], "fail_" + settings["FailFlag"] + ".txt"), 'w').close()

    _logger.debug("Final settings:\\ {}".format(settings))

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
        run_volume_conductor_model(settings, volume_conductor)

    time_1 = time.time()
    timings["VolumeConductor"] = time_1 - time_0

    # TODO continue
    # Once the TimeResult object is implemented, you should no longer need to work
    # with the grid_pts, it should remain hidden behind the ossdbs interface
    # (probably implemented in a PointModel object)
    vl = generate_neuron_grid(settings)
    grid_pts = vl.coordinates()
    if settings["PointModel"]["Lattice"]["Active"] is True:
        grid_pts = vl.coordinates()
    else:
        grid_pts = vl.coordinates

    time_1 = time.time()
    timings["LatticeModel"] = time_1 - time_0
    time_0 = time_1

    potentials = volume_conductor.evaluate_potential_at_points(grid_pts)
    fields = volume_conductor.evaluate_field_at_points(grid_pts)
    field_mags = np.linalg.norm(fields, axis=1).reshape((fields.shape[0], 1))

    time_1 = time.time()
    timings["FieldProbing"] = time_1 - time_0
    time_0 = time_1

    # For high res use Latice
    # If computed in MNI, create affine
    if settings["PointModel"]['VoxelLattice']['Active']:
        vl.save_as_nifti(settings, field_mags, os.path.join(settings["OutputPath"], "E_field_solution.nii"))
        vl.save_as_nifti(settings, field_mags, os.path.join(settings["OutputPath"], "VTA_solution.nii"), binarize=True)

    # Save points
    h5f_pts = h5py.File(os.path.join(settings["OutputPath"], "oss_pts.h5"), 'w')
    h5f_pts.create_dataset("points", data=grid_pts)
    h5f_pts.close()

    # Save potential evaluation
    h5f_pot = h5py.File(os.path.join(settings["OutputPath"], 'oss_potentials.h5'), 'w')
    h5f_pot.create_dataset("points", data=grid_pts)
    h5f_pot.create_dataset("potentials", data=potentials)
    h5f_pot.close()
    df_pot = pd.DataFrame(np.concatenate([grid_pts, potentials.reshape((potentials.shape[0], 1))], axis=1),
                          columns=["x-pt", "y-pt", "z-pt", "potential"])
    df_pot.to_csv(os.path.join(settings["OutputPath"], "oss_potentials.csv"), index=False)

    # Save electric field evaluation
    h5f_field = h5py.File(os.path.join(settings["OutputPath"], "oss_field.h5"), 'w')
    h5f_field.create_dataset("points", data=grid_pts)
    h5f_field.create_dataset("field/field_vecs", data=fields)
    h5f_field.create_dataset("field/field_mags", data=field_mags)
    h5f_field.close()
    df_field = pd.DataFrame(np.concatenate([grid_pts, fields, field_mags], axis=1),
                            columns=["x-pt", "y-pt", "z-pt", "x-field", "y-field", "z-field", "magnitude"])
    if settings["TemplateSpace"]:
        df_field.to_csv(os.path.join(settings["OutputPath"], "E_field_Template_space.csv"), index=False)
    else:
        df_field.to_csv(os.path.join(settings["OutputPath"], "E_field_MRI_space.csv"), index=False)

    time_1 = time.time()
    timings["FieldExport"] = time_1 - time_0
    time_0 = time_1

    _logger.info("Timings:\n {}".format(pprint.pformat(timings)))

    # write success file
    open(os.path.join(settings["OutputPath"], "success_" + settings["FailFlag"] + ".txt"), 'w').close()
    os.remove(os.path.join(settings["OutputPath"], "fail_" + settings["FailFlag"] + ".txt"))
    _logger.info("Process Completed")


if __name__ == '__main__':
    main()
