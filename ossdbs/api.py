# Copyright 2023, 2024 Konstantin Butenko, Shruthi Chakravarthy
# Copyright 2023, 2024 Jan Philipp Payonk, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import importlib
import json
import logging
import os
from typing import Optional

import numpy as np

from ossdbs.dielectric_model import (
    default_dielectric_parameters,
    dielectric_model_parameters,
    dielectric_models,
)
from ossdbs.electrodes import ELECTRODE_MODELS, ELECTRODE_PARAMETERS, ELECTRODES
from ossdbs.fem import (
    PRECONDITIONERS,
    SOLVERS,
    Mesh,
    VolumeConductor,
    VolumeConductorFloating,
    VolumeConductorFloatingImpedance,
    VolumeConductorNonFloating,
)
from ossdbs.model_geometry import BoundingBox, BrainGeometry, ModelGeometry
from ossdbs.point_analysis import Lattice, Pathway, VoxelLattice
from ossdbs.stimulation_signals import (
    FrequencyDomainSignal,
    RectangleSignal,
    TimeDomainSignal,
    TrapezoidSignal,
    TriangleSignal,
)
from ossdbs.stimulation_signals.utilities import get_positive_frequencies
from ossdbs.utils.nifti1image import DiffusionTensorImage, MagneticResonanceImage

_logger = logging.getLogger(__name__)

PAM_AVAILABLE = importlib.util.find_spec("neuron") is not None
if not PAM_AVAILABLE:
    _logger.warning("NEURON is not installed, disabling PAM analysis!")


def create_bounding_box(box_parameters: dict) -> BoundingBox:
    """Create a bounding box around a domain in space.

    Notes
    -----
    `box_parameters` need to contain the center (three coordinates,
    all following the style `x[mm]` (and similar for y and z-component))
    and the outer dimensions.
    """
    input_s = box_parameters["Dimension"]
    input_c = box_parameters["Center"]
    shape = (input_s["x[mm]"], input_s["y[mm]"], input_s["z[mm]"])
    center = (input_c["x[mm]"], input_c["y[mm]"], input_c["z[mm]"])
    start = center - np.divide(shape, 2)
    end = start + shape
    return BoundingBox(tuple(start), tuple(end))


def generate_electrodes(settings: dict):
    """Generate an OCC electrode model from the settings dict."""
    _logger.info("Generate electrode geometries")

    hp_refinement = False
    if "HPRefinement" in settings["Mesh"]:
        hp_refinement = settings["Mesh"]["HPRefinement"]["Active"]
    electrodes = []
    for electrode_parameters in settings["Electrodes"]:
        name = electrode_parameters["Name"]
        direction = (
            electrode_parameters["Direction"]["x[mm]"],
            electrode_parameters["Direction"]["y[mm]"],
            electrode_parameters["Direction"]["z[mm]"],
        )
        rotation = electrode_parameters["Rotation[Degrees]"]
        position = (
            electrode_parameters["TipPosition"]["x[mm]"],
            electrode_parameters["TipPosition"]["y[mm]"],
            electrode_parameters["TipPosition"]["z[mm]"],
        )

        # Implemented custom electrodes without using custom_electrodes.py
        if "Custom" in name:
            electrode_model = ELECTRODE_MODELS[name]
            parameter_class = ELECTRODE_PARAMETERS[electrode_model.__name__]
            custom_list = electrode_parameters["CustomParameters"]
            electrode = electrode_model(
                parameters=parameter_class(**custom_list),
                direction=direction,
                position=position,
                rotation=rotation,
            )

        else:
            electrode_type = ELECTRODES[name]
            electrode = electrode_type(
                direction=direction,
                position=position,
                rotation=rotation,
            )

        if hp_refinement:
            electrode.set_hp_flag(electrode_parameters=electrode_parameters)

        if "EncapsulationLayer" in electrode_parameters:
            electrode.encapsulation_thickness = electrode_parameters[
                "EncapsulationLayer"
            ]["Thickness[mm]"]
        electrodes.append(electrode)

    if settings["ExportElectrode"]:
        n_electrode = 0
        for electrode in electrodes:
            n_electrode = n_electrode + 1
            electrode.export_electrode(
                settings["OutputPath"], settings["BrainRegion"], n_electrode
            )

    return electrodes


def prepare_dielectric_properties(settings: dict) -> dict:
    """Return dictionary with dielectric properties for each tissue."""
    _logger.info("Prepare dielectric model")
    dielectric_settings = settings["DielectricModel"]
    model_type = dielectric_settings["Type"]
    custom_parameters = dielectric_settings["CustomParameters"]

    # create empty dict for collection of dielectric models
    dielectric_properties = {}
    dielectric_model = dielectric_models[model_type]
    parameter_template = dielectric_model_parameters[model_type]
    default_parameters = default_dielectric_parameters[model_type]
    for material in settings["MaterialDistribution"]["MRIMapping"]:
        if custom_parameters is not None:
            model_parameters = parameter_template(**custom_parameters[material])
        else:
            model_parameters = default_parameters[material]
        dielectric_properties[material] = dielectric_model(model_parameters)
    return dielectric_properties


def generate_brain_model(settings, rotate_initial_geo: bool = False):
    """Generate OCC brain model."""
    brain_region_parameters = settings["BrainRegion"]
    brain_shape = brain_region_parameters["Shape"]
    brain_region = create_bounding_box(brain_region_parameters)
    brain_model = BrainGeometry(
        brain_shape, brain_region, rotate_initial_geo=rotate_initial_geo
    )
    return brain_model


def generate_model_geometry(settings):
    """Generate a full geometry comprising brain and electrodes."""
    brain = generate_brain_model(settings)
    electrodes = generate_electrodes(settings)
    try:
        model_geometry = ModelGeometry(brain, electrodes)
    except RuntimeError:
        _logger.warning(
            "Could not build geometry, trying again after rotation of initial geometry"
        )
        brain = generate_brain_model(settings, rotate_initial_geo=True)
        model_geometry = ModelGeometry(brain, electrodes)
    return model_geometry


def build_brain_model(
    settings,
    mri_image: Optional[MagneticResonanceImage] = None,
    rotate_initial_geo: bool = False,
) -> BrainGeometry:
    """Build geometry model of brain."""
    # MRI image is default choice for brain construction
    if "BrainRegion" in settings:
        _logger.debug("Generating model geometry for fixed brain region")
        region_parameters = settings["BrainRegion"]
        brain_region = create_bounding_box(region_parameters)
        shape = settings["BrainRegion"]["Shape"]
        return BrainGeometry(shape, brain_region, rotate_initial_geo=rotate_initial_geo)
    else:
        _logger.debug("Generating model geometry from MRI image")
        if mri_image is None:
            raise ValueError("Need to provide MRI image to build geo.")
        # attention: bounding box is given in voxel space!
        brain_region = mri_image.bounding_box
        shape = "Ellipsoid"
        # transformation to real space in geometry creation
        _logger.debug(
            "Generate OCC model, passing transformation matrix from MRI image"
        )
        return BrainGeometry(
            shape,
            brain_region,
            trafo_matrix=mri_image.trafo_matrix,
            translation=mri_image.translation,
            rotate_intial_geo=rotate_initial_geo,
        )


def set_contact_and_encapsulation_layer_properties(settings, model_geometry):
    """Update boundary and material values on contacts and encapsulation layers."""
    _logger.info("Set values on contacts and encapsulation layers")
    electrode_settings = settings["Electrodes"]
    offset = 0
    for idx, new_parameters in enumerate(electrode_settings):
        _logger.debug(f"Update Electrode {idx} with settings {new_parameters}")
        if "Contacts" in new_parameters:
            for contact_info in new_parameters["Contacts"]:
                contact_idx = offset + contact_info["Contact_ID"]
                # contacts are zero-indexed in the model_geometry
                model_geometry.update_contact(contact_idx - 1, contact_info)
            offset += model_geometry.electrodes[idx].n_contacts
        if "EncapsulationLayer" in new_parameters:
            # encapsulation layer is one-indexed in the model_geometry
            _logger.debug(f"Updating encapsulation layer {idx + 1}")
            encap_idx = model_geometry.get_encapsulation_layer_index(
                f"EncapsulationLayer_{idx + 1}"
            )
            _logger.debug(f"Encapsulation layer has index {encap_idx}")
            if encap_idx != -1:
                _logger.info("Updating encapsulation layer properties")
                model_geometry.update_encapsulation_layer(
                    encap_idx, new_parameters["EncapsulationLayer"]
                )
    if "Surfaces" in settings:
        for surface in settings["Surfaces"]:
            idx = model_geometry.get_contact_index(surface["Name"])
            if idx == -1:
                raise ValueError(
                    "Surface {} not part of the geometry".format(surface["Name"])
                )
            model_geometry.update_contact(idx, surface)


def set_custom_mesh_sizes(settings, model_geometry):
    """Update the mesh sizes."""
    model_geometry.set_mesh_sizes(settings["Mesh"]["MeshSize"])


def generate_mesh(settings):
    """Generate a mesh from settings.

    Notes
    -----
    Attention! This mesh is not yet curved!
    """
    model_geometry = generate_model_geometry(settings)
    set_contact_and_encapsulation_layer_properties(settings, model_geometry)
    if "MeshSize" in settings["Mesh"]:
        set_custom_mesh_sizes(settings, model_geometry)

    mesh_settings = settings["Mesh"]
    mesh_order = 1

    mesh = Mesh(model_geometry.geometry, mesh_order)
    if mesh_settings["LoadMesh"]:
        mesh.load_mesh(mesh_settings["LoadPath"])
        return mesh

    if "MeshingHypothesis" not in mesh_settings:
        mesh_settings["MeshingHypothesis"] = {"Type": "Default"}
    mesh.generate_mesh(mesh_settings)
    if mesh_settings["SaveMesh"]:
        mesh.save(mesh_settings["SavePath"])
    return mesh


def prepare_solver(settings):
    """Set up solver and preconditioner."""
    _logger.info("Preparing solver")
    parameters = settings["Solver"]
    solver_type = parameters["Type"]
    solver = SOLVERS[solver_type]
    preconditioner_kwargs = parameters["PreconditionerKwargs"]
    preconditioner = PRECONDITIONERS[parameters["Preconditioner"]](
        **preconditioner_kwargs
    )

    return solver(
        precond_par=preconditioner,
        maxsteps=parameters["MaximumSteps"],
        precision=parameters["Precision"],
    )


def generate_point_models(settings: dict):
    """Generate a list of point models."""
    point_models = []
    if settings["PointModel"]["Pathway"]["Active"]:
        file_name = settings["PointModel"]["Pathway"]["FileName"]
        _logger.info(f"Import neuron geometries stored in {file_name}")
        export_field = settings["PointModel"]["Pathway"]["ExportField"]
        point_models.append(Pathway(file_name, export_field=export_field))
    if settings["PointModel"]["Lattice"]["Active"]:
        shape_par = settings["PointModel"]["Lattice"]["Shape"]
        shape = shape_par["x"], shape_par["y"], shape_par["z"]
        center_par = settings["PointModel"]["Lattice"]["Center"]
        center = center_par["x[mm]"], center_par["y[mm]"], center_par["z[mm]"]
        dir_par = settings["PointModel"]["Lattice"]["Direction"]
        direction = dir_par["x[mm]"], dir_par["y[mm]"], dir_par["z[mm]"]
        distance = settings["PointModel"]["Lattice"]["PointDistance[mm]"]
        collapse_vta = settings["PointModel"]["Lattice"]["CollapseVTA"]
        export_field = settings["PointModel"]["Lattice"]["ExportField"]

        point_models.append(
            Lattice(
                shape=shape,
                center=center,
                distance=distance,
                direction=direction,
                collapse_vta=collapse_vta,
                export_field=export_field,
            )
        )

    if settings["PointModel"]["VoxelLattice"]["Active"]:
        _logger.info("from voxel lattice")
        center_par = settings["PointModel"]["Lattice"]["Center"]
        center = center_par["x[mm]"], center_par["y[mm]"], center_par["z[mm]"]
        mri_image = MagneticResonanceImage(settings["MaterialDistribution"]["MRIPath"])
        affine = mri_image.affine
        header = mri_image.header
        shape_par = settings["PointModel"]["VoxelLattice"]["Shape"]
        shape = np.array([shape_par["x"], shape_par["y"], shape_par["z"]])
        export_field = settings["PointModel"]["VoxelLattice"]["ExportField"]
        point_models.append(
            VoxelLattice(center, affine, shape, header, export_field=export_field)
        )
    return point_models


def generate_meshsize_file_from_neuron_grid():
    """Use point grid to specify mesh sizes.

    TODO implement
    """
    raise NotImplementedError("Not yet supported")


def generate_signal(settings) -> TimeDomainSignal:
    """Generate a time-domain signal (waveform)."""
    signal_settings = settings["StimulationSignal"]
    signal_type = signal_settings["Type"]
    if signal_type == "Rectangle":
        signal = RectangleSignal(
            signal_settings["Frequency[Hz]"],
            1e-6 * signal_settings["PulseWidth[us]"],
            1e-6 * signal_settings["InterPulseWidth[us]"],
            1e-6 * signal_settings["CounterPulseWidth[us]"],
            signal_settings["CounterAmplitude"],
        )
    elif signal_type == "Triangle":
        signal = TriangleSignal(
            signal_settings["Frequency[Hz]"],
            1e-6 * signal_settings["PulseWidth[us]"],
            1e-6 * signal_settings["InterPulseWidth[us]"],
            1e-6 * signal_settings["CounterPulseWidth[us]"],
            signal_settings["CounterAmplitude"],
        )
    elif signal_type == "Trapezoid":
        signal = TrapezoidSignal(
            signal_settings["Frequency[Hz]"],
            1e-6 * signal_settings["PulseWidth[us]"],
            1e-6 * signal_settings["InterPulseWidth[us]"],
            1e-6 * signal_settings["CounterPulseWidth[us]"],
            1e-6 * signal_settings["PulseTopWidth[us]"],
            signal_settings["CounterAmplitude"],
        )
    signal.plot_time_domain_signal(
        signal_settings["CutoffFrequency"], settings["OutputPath"]
    )
    return signal


def prepare_volume_conductor_model(
    settings, model_geometry, conductivity, solver
) -> VolumeConductor:
    """Prepare the volume conductor model."""
    _logger.info("Generate volume conductor model")
    order = settings["FEMOrder"]

    mesh_parameters = settings["Mesh"]
    floating_mode = model_geometry.get_floating_mode()
    output_path = settings["OutputPath"]
    _logger.info(f"Output path set to: {output_path}")
    if floating_mode == "Floating":
        _logger.debug("Floating mode selected")
        return VolumeConductorFloating(
            model_geometry, conductivity, solver, order, mesh_parameters, output_path
        )

    elif floating_mode == "FloatingImpedance":
        _logger.debug("FloatingImpedance mode selected")
        return VolumeConductorFloatingImpedance(
            model_geometry, conductivity, solver, order, mesh_parameters, output_path
        )
    _logger.debug("Non floating mode selected")
    return VolumeConductorNonFloating(
        model_geometry, conductivity, solver, order, mesh_parameters, output_path
    )


def prepare_stimulation_signal(settings) -> FrequencyDomainSignal:
    """Prepare the frequency-domain representation of stimulation signal."""
    signal_settings = settings["StimulationSignal"]
    signal_type = signal_settings["Type"]
    current_controlled = signal_settings["CurrentControlled"]
    octave_band_approximation = False
    if signal_type == "Multisine":
        frequencies = signal_settings["ListOfFrequencies"]
        fourier_coefficients = np.ones(len(frequencies))
        base_frequency = frequencies[0]
        cutoff_frequency = frequencies[0]
        signal_length = len(frequencies)
    else:
        spectrum_mode = signal_settings["SpectrumMode"]
        if spectrum_mode == "OctaveBand":
            octave_band_approximation = True

        signal = generate_signal(settings)
        cutoff_frequency = signal_settings["CutoffFrequency"]
        base_frequency = signal.frequency
        fft_frequencies, fft_coefficients = signal.get_fft_spectrum(cutoff_frequency)
        signal_length = len(fft_coefficients)
        frequencies, fourier_coefficients = get_positive_frequencies(
            fft_frequencies, fft_coefficients
        )

    frequency_domain_signal = FrequencyDomainSignal(
        frequencies=frequencies,
        amplitudes=fourier_coefficients,
        current_controlled=current_controlled,
        base_frequency=base_frequency,
        cutoff_frequency=cutoff_frequency,
        signal_length=signal_length,
        octave_band_approximation=octave_band_approximation,
    )
    return frequency_domain_signal


def run_volume_conductor_model(
    settings, volume_conductor, frequency_domain_signal, truncation_time=None
):
    """TODO document.


    Notes
    -----
    Run at all frequencies.
    If the mode is multisine, a provided list of frequencies is used.
    """
    _logger.info("Run volume conductor model")

    out_of_core = settings["OutOfCore"]
    compute_impedance = False
    if "ComputeImpedance" in settings:
        if settings["ComputeImpedance"]:
            _logger.info("Will compute impedance at each frequency")
            compute_impedance = True
    if "ExportVTK" in settings:
        export_vtk = settings["ExportVTK"]
        if export_vtk:
            _logger.info("Will export solution to VTK")
    else:
        export_vtk = False
    if "ExportFrequency" in settings:
        export_frequency = settings["ExportFrequency"]
        if export_frequency is not None:
            _logger.info(f"Set custom export frequency to {export_frequency}.")

    point_models = generate_point_models(settings)

    vcm_timings = volume_conductor.run_full_analysis(
        frequency_domain_signal,
        compute_impedance,
        export_vtk,
        point_models=point_models,
        activation_threshold=settings["ActivationThresholdVTA[V-per-m]"],
        out_of_core=out_of_core,
        export_frequency=export_frequency,
        adaptive_mesh_refinement_settings=settings["Mesh"]["AdaptiveMeshRefinement"],
        material_mesh_refinement_steps=settings["Mesh"]["MaterialRefinementSteps"],
        truncation_time=truncation_time,
    )
    return vcm_timings


def run_stim_sets(settings, geometry, conductivity, solver, frequency_domain_signal):
    """TODO document.

    Notes
    -----
    Run at all frequencies.
    If the mode is multisine, a provided list of frequencies is used.
    """
    _logger.info("Run StimSets volume conductor model")

    out_of_core = settings["OutOfCore"]
    if not frequency_domain_signal.current_controlled:
        _logger.warning(
            "StimSets requires current-controlled stimulation"
            ", thus the setting was switched on"
        )
    # no vtk export
    export_vtk = settings["ExportVTK"]
    # no intermediate exports
    export_frequency = None
    # no VTA analysis
    activation_threshold = settings["ActivationThresholdVTA[V-per-m]"]
    # prepare point model
    point_models = generate_point_models(settings)

    ground_contact = None
    for contact in geometry.contacts:
        if np.isclose(contact.current, -1) and contact.active:
            ground_contact = contact.name
            _logger.info(f"Will skip ground contact {contact.name}")
    if ground_contact is None:
        raise ValueError(
            "No ground contact set. Choose one active contact with current -1."
        )
    for contact in geometry.contacts:
        if contact.name == ground_contact:
            continue
        # set current contact active, all other passive
        for upd_contact in geometry.contacts:
            # reset all voltages
            contact_idx = geometry.get_contact_index(upd_contact.name)
            geometry.update_contact(contact_idx, {"Voltage[V]": 0.0})
            # don't change ground
            if upd_contact.name == ground_contact:
                continue
            active = False
            floating = True
            current = 0.0
            voltage = False
            if contact.name == upd_contact.name:
                active = True
                floating = False
                current = 1.0
                voltage = 1.0
            # write new contact settings
            geometry.update_contact(
                contact_idx,
                {
                    "Floating": floating,
                    "Active": active,
                    "Current[A]": current,
                    "Voltage[V]": voltage,
                },
            )
        volume_conductor = prepare_volume_conductor_model(
            settings, geometry, conductivity, solver
        )
        _logger.info(f"Running with contacts:\n{volume_conductor.contacts}")

        volume_conductor.output_path = settings["OutputPath"] + contact.name
        vcm_timings = volume_conductor.run_full_analysis(
            frequency_domain_signal,
            export_vtk=export_vtk,
            point_models=point_models,
            activation_threshold=activation_threshold,
            out_of_core=out_of_core,
            export_frequency=export_frequency,
            adaptive_mesh_refinement_settings=settings["Mesh"][
                "AdaptiveMeshRefinement"
            ],
            material_mesh_refinement_steps=settings["Mesh"]["MaterialRefinementSteps"],
        )
        _logger.info(f"Timing for contact {contact.name}: {vcm_timings}")


def load_images(settings):
    """Load MRI and DTI images."""
    _logger.info("Load MRI image")
    mri_path = settings["MaterialDistribution"]["MRIPath"]
    _logger.debug(f"Input path: {mri_path}")
    mri_image = MagneticResonanceImage(mri_path)
    dti_image = None
    if settings["MaterialDistribution"]["DiffusionTensorActive"]:
        _logger.info("Load DTI image")
        dti_image = DiffusionTensorImage(settings["MaterialDistribution"]["DTIPath"])
    return mri_image, dti_image


def run_PAM(settings):
    """Run pathway activation analysis."""
    if not PAM_AVAILABLE:
        raise RuntimeError("PAM not available! Please install NEURON!")
    from ossdbs.axon_processing import get_neuron_model

    _logger.info("Running PAM")
    pathway_file = settings["PathwayFile"]
    pathway_solution_dir = settings["OutputPath"]
    time_domain_solution = os.path.join(
        settings["OutputPath"], "oss_time_result_PAM.h5"
    )
    with open(pathway_file) as fp:
        pathways_dict = json.load(fp)

    model_type = pathways_dict["Axon_Model_Type"]
    neuron_model = get_neuron_model(model_type, pathways_dict, pathway_solution_dir)

    if settings["StimSets"]["Active"]:
        if "CurrentVector" not in settings:
            settings["CurrentVector"] = None
        # files to load individual solutions from
        time_domain_solution_files = []

        if settings["StimSets"]["StimSetsFile"] is not None:
            _logger.info("Load current vectors form file.")
            stim_protocols = np.genfromtxt(
                settings["StimSets"]["StimSetsFile"],
                dtype=float,
                delimiter=",",
                names=True,
            )
            n_stim_protocols = stim_protocols.shape[0]
            n_contacts = len(list(stim_protocols[0]))
        else:
            if settings["CurrentVector"] is None:
                raise ValueError("Provide either a StimSetsFile or a CurrentVector")
            n_stim_protocols = 1
            # load current from input file
            stim_protocols = [settings["CurrentVector"]]
            # assign contacts
            n_contacts = len(stim_protocols[0])

        # load unit solutions once
        _logger.info("Load unit solutions")
        for contact_i in range(n_contacts):
            time_domain_solution_files.append(
                os.path.join(
                    settings["OutputPath"] + f"E1C{contact_i + 1}",
                    "oss_time_result_PAM.h5",
                )
            )

        td_unit_solutions = neuron_model.load_unit_solutions(time_domain_solution_files)

        # go through stimulation protocols
        _logger.info("Running stimulation protocols")
        for protocol_i in range(n_stim_protocols):
            # get the scaling vector for the current
            scaling_vector = list(stim_protocols[protocol_i])
            # swap NaNs to zero current and convert to A (StimSets in mA)
            scaling_vector = [0 if np.isnan(x) else 1e-3 * x for x in scaling_vector]

            td_solution = neuron_model.superimpose_unit_solutions(
                td_unit_solutions, scaling_vector
            )
            # when using optimizer, scaling_index is not used
            if (
                settings["CurrentVector"] is not None
                and settings["StimSets"]["StimSetsFile"] is None
            ):
                neuron_model.process_pathways(
                    td_solution, scaling=settings["Scaling"], scaling_index=None
                )
            else:
                neuron_model.process_pathways(
                    td_solution, scaling=settings["Scaling"], scaling_index=protocol_i
                )
    else:
        td_solution = neuron_model.load_solution(time_domain_solution)
        neuron_model.process_pathways(
            td_solution,
            scaling=settings["Scaling"],
            scaling_index=settings["ScalingIndex"],
        )
