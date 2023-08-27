from h5py._hl.base import default_lapl
import numpy as np
import pandas as pd
import os
import sys
import json
import h5py
import ossdbs
from ossdbs.leaddbs_interface.imp_coord import imp_coord
import logging
import time
#import ngsolve
#ngsolve.ngsglobals.msg_level = 10

_logger = logging.getLogger(__name__)

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
    settings["StimulationSignal"]["Type"] = "Multisine"
    settings["StimulationSignal"]["ListOfFrequencies"] = [10000]

    settings["PointModel"]["VoxelLattice"] = {"Active": False,
                                              "Shape": {'x':31,'y':31,'z':41}}
    settings["Electrodes"][settings["ModelSide"]]["EncapsulationLayer"]["Thickness"] = 0.1
    settings["Electrodes"][settings["ModelSide"]]["EncapsulationLayer"]["Material"] = "CSF"
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
                                             "x": 30,
                                             "y": 30,
                                             "z": 40
                                         },
                                         "Direction": {
                                             "x[mm]": settings["Electrodes"][0]["Direction"]["x[mm]"],
                                             "y[mm]": settings["Electrodes"][0]["Direction"]["y[mm]"],
                                             "z[mm]": settings["Electrodes"][0]["Direction"]["z[mm]"]
                                         },
                                         "PointDistance[mm]": 0.5
                                         }

    return settings

def main(lead_settings, hemi_side):

    timings = {}
    time_0 = time.time()
    hemi_side = int(hemi_side)

    ### Lead/OSS Interface ###
    output_path, tail = os.path.split(lead_settings)
    # sys.path.insert(output_path)
    os.environ["STIMFOLDER"] = output_path
    SIDES = ['rh', 'lh']
    side = SIDES[hemi_side]

    # get settings from oss-dbs_parameters.mat
    settings = ossdbs.load_from_lead_mat(lead_settings, hemi_side)

    # add default settings (alternatively, set using GUI)
    settings = load_default_for_lead(settings)

    # create the fail flag for Lead-DBS, will be removed if successful run
    with open(os.path.join(os.environ["STIMFOLDER"], "fail_" + side + ".txt"), 'w') as f:
        print("Settings imported")

    if not os.path.isdir(settings["OutputPath"]):
        os.mkdir(settings["OutputPath"])
    # Save settings as a json
    with open(os.path.join(settings["OutputPath"], "settings.json"), 'w') as f:
        json.dump(settings, f)

    time_1 = time.time()
    timings["Importing settings"] = time_1 - time_0
    time_0 = time_1

    ### OSS Processing ###

    electrodes = ossdbs.generate_electrodes(settings)
    electrode = electrodes[0]  # why 0 here?
    # Draw(electrode.geometry)
    time_1 = time.time()
    timings["Electrodes"] = time_1 - time_0
    time_0 = time_1

    mri_path = settings['MaterialDistribution']['MRIPath']

    mri_image = ossdbs.MagneticResonanceImage(os.path.join(os.environ["STIMFOLDER"], mri_path))
    dti_image = None
    if settings["MaterialDistribution"]["DiffusionTensorActive"]:
        dti_image = ossdbs.DiffusionTensorImage(os.path.join(os.environ["STIMFOLDER"],settings["MaterialDistribution"]["DTIPath"]))

    time_1 = time.time()
    timings["Importing Imaging Data"] = time_1 - time_0
    time_0 = time_1


    # brain_region = mri_image.bounding_box

    roi_bb = ossdbs.create_bounding_box(settings['BrainRegion'])

    brain = ossdbs.BrainGeometry("Ellipsoid", roi_bb)
    # Draw(brain.geometry)

    model_geometry = ossdbs.ModelGeometry(brain, electrodes)

    time_1 = time.time()
    timings["Model Geometry"] = time_1 - time_0
    time_0 = time_1

    ossdbs.set_contact_and_encapsulation_layer_properties(settings, model_geometry)

    time_1 = time.time()
    timings["ElectrodeProperties"] = time_1 - time_0
    time_0 = time_1

    dielectric_model = ossdbs.prepare_dielectric_properties(settings)

    time_1 = time.time()
    timings["DielectricModel"] = time_1 - time_0
    time_0 = time_1

    materials = settings["MaterialDistribution"]["MRIMapping"]
    conductivity = ossdbs.ConductivityCF(mri_image,
                                         roi_bb,
                                         dielectric_model,
                                         materials,
                                         model_geometry.encapsulation_layers,
                                         complex_data=settings["EQSMode"],
                                         dti_image=dti_image)

    time_1 = time.time()
    timings["ConductivityCF"] = time_1 - time_0
    time_0 = time_1


    solver = ossdbs.prepare_solver(settings)
    volume_conductor = ossdbs.prepare_volume_conductor_model(settings, model_geometry, conductivity, solver)
    ossdbs.run_volume_conductor_model(settings, volume_conductor)
    # DrawNG(volume_conductor.potential)
    # DrawNG(volume_conductor.electric_field, volume_conductor.mesh.ngsolvemesh)

    time_1 = time.time()
    timings["VolumeConductor"] = time_1 - time_0
    time_0 = time_1

    # Once the TimeResult object is implemented, you should no longer need to work
    # with the grid_pts, it should remain hidden behind the ossdbs interface
    # (probably implemented in a PointModel object)
    vl = ossdbs.generate_neuron_grid(settings)
    grid_pts = vl.coordinates()
    if settings["PointModel"]["Lattice"]["Active"] == True:
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
        vl.save_as_nifti(settings, field_mags, "E_field_solution.nii")
        vl.save_as_nifti(settings, field_mags, "VTA_solution.nii", binarize=True)


    #elif settings["PointModel"]["Lattice"]["Active"] and settings["TemplateSpace"]:


    # from ngsolve import VTKOutput
    # vtk = VTKOutput(ma=volume_conductor.mesh.ngsolvemesh,
    #            coefs=[conductivity.material_distribution],
    #            names = ["mat_distro"],
    #            filename=settings["OutputPath"]+"/mat_distro",
    #            subdivision=3)
    # vtk.Do()

    # Save points
    h5f_pts = h5py.File(settings["OutputPath"] + "/oss_pts.h5", 'w')
    h5f_pts.create_dataset("points", data=grid_pts)
    h5f_pts.close()
    df_pts = pd.DataFrame(grid_pts, columns=['x', 'y', 'z'])


        #settings.Estimate_In_Template

    # Save potential evaluation
    h5f_pot = h5py.File(settings["OutputPath"] + '/oss_potentials.h5', 'w')
    h5f_pot.create_dataset("points", data=grid_pts)
    h5f_pot.create_dataset("potentials", data=potentials)
    h5f_pot.close()
    df_pot = pd.DataFrame(np.concatenate([grid_pts, potentials.reshape((potentials.shape[0], 1))], axis=1),
                          columns=["x-pt", "y-pt", "z-pt", "potential"])
    df_pot.to_csv(settings["OutputPath"] + "/oss_potentials.csv", index=False)

    # Save electric field evaluation
    h5f_field = h5py.File(settings["OutputPath"] + "/oss_field.h5", 'w')
    h5f_field.create_dataset("points", data=grid_pts)
    h5f_field.create_dataset("field/field_vecs", data=fields)
    h5f_field.create_dataset("field/field_mags", data=field_mags)
    h5f_field.close()
    df_field = pd.DataFrame(np.concatenate([grid_pts, fields, field_mags], axis=1),
                            columns=["x-pt", "y-pt", "z-pt", "x-field", "y-field", "z-field", "magnitude"])
    if settings["TemplateSpace"]:
        df_field.to_csv(settings["OutputPath"] + "/E_field_Template_space.csv", index=False)
    else:
        df_field.to_csv(settings["OutputPath"] + "/E_field_MRI_space.csv", index=False)

    time_1 = time.time()
    timings["FieldExport"] = time_1 - time_0
    time_0 = time_1

    import pprint
    # Prints the nicely formatted dictionary
    pprint.pprint(timings)

    with open(os.path.join(os.environ["STIMFOLDER"], "success_" + side + ".txt"), 'w') as f:
        os.remove(os.path.join(os.environ["STIMFOLDER"], "fail_" + side + ".txt"))
        print("Process Completed")


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
