import logging

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
from ossdbs.point_analysis import Lattice, Pathway, PointModel, VoxelLattice
from ossdbs.stimulation_signals import (
    FrequencyDomainSignal,
    RectangleSignal,
    TimeDomainSignal,
    TrapezoidSignal,
    TriangleSignal,
)
from ossdbs.utils.nifti1image import DiffusionTensorImage, MagneticResonanceImage

_logger = logging.getLogger(__name__)


def create_bounding_box(box_parameters: dict) -> BoundingBox:
    input_s = box_parameters["Dimension"]
    input_c = box_parameters["Center"]
    shape = (input_s["x[mm]"], input_s["y[mm]"], input_s["z[mm]"])
    center = (input_c["x[mm]"], input_c["y[mm]"], input_c["z[mm]"])
    start = center - np.divide(shape, 2)
    end = start + shape
    return BoundingBox(tuple(start), tuple(end))


def generate_electrodes(settings: dict):
    """Generate OCC electrode models.

    Notes
    -----
    """
    _logger.info("Generate electrode geometries")
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


def generate_brain_model(settings):
    """Generate OCC brain model.

    Notes
    -----
    TODO type checking

    """
    brain_region_parameters = settings["BrainRegion"]
    brain_shape = brain_region_parameters["Shape"]
    brain_region = create_bounding_box(brain_region_parameters)
    brain_model = BrainGeometry(brain_shape, brain_region)
    return brain_model


def generate_model_geometry(settings):
    brain = generate_brain_model(settings)
    electrodes = generate_electrodes(settings)
    model_geometry = ModelGeometry(brain, electrodes)
    return model_geometry


def set_contact_and_encapsulation_layer_properties(settings, model_geometry):
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
    model_geometry.set_mesh_sizes(settings["Mesh"]["MeshSize"])


def generate_mesh(settings):
    """

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

    if "MeshingHypothesis" in mesh_settings:
        mesh_hypothesis = mesh_settings["MeshingHypothesis"]
    else:
        mesh_hypothesis = {"Type": "Default"}
    mesh.generate_mesh(mesh_hypothesis)
    if mesh_settings["SaveMesh"]:
        mesh.save(mesh_settings["SavePath"])
    return mesh


def prepare_solver(settings):
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
        printrates=parameters["PrintRates"],
        maxsteps=parameters["MaximumSteps"],
        precision=parameters["Precision"],
    )


def generate_neuron_grid(settings: dict) -> PointModel:
    if settings["PointModel"]["Pathway"]["Active"]:
        file_name = settings["PointModel"]["Pathway"]["FileName"]
        _logger.info(f"Import neuron geometries stored in {file_name}")
        return Pathway(file_name)
    elif settings["PointModel"]["Lattice"]["Active"]:
        shape_par = settings["PointModel"]["Lattice"]["Shape"]
        shape = shape_par["x"], shape_par["y"], shape_par["z"]
        center_par = settings["PointModel"]["Lattice"]["Center"]
        center = center_par["x[mm]"], center_par["y[mm]"], center_par["z[mm]"]
        dir_par = settings["PointModel"]["Lattice"]["Direction"]
        direction = dir_par["x[mm]"], dir_par["y[mm]"], dir_par["z[mm]"]
        distance = settings["PointModel"]["Lattice"]["PointDistance[mm]"]
        collapse_vta = settings["PointModel"]["Lattice"]["CollapseVTA"]

        _logger.info("from lattice")
        return Lattice(
            shape=shape,
            center=center,
            distance=distance,
            direction=direction,
            collapse_vta=collapse_vta,
        )
    elif settings["PointModel"]["VoxelLattice"]["Active"]:
        _logger.info("from voxel lattice")
        center_par = settings["PointModel"]["Lattice"]["Center"]
        center = center_par["x[mm]"], center_par["y[mm]"], center_par["z[mm]"]
        mri_image = MagneticResonanceImage(settings["MaterialDistribution"]["MRIPath"])
        affine = mri_image.affine
        header = mri_image.header
        shape_par = settings["PointModel"]["VoxelLattice"]["Shape"]
        shape = np.array([shape_par["x"], shape_par["y"], shape_par["z"]])
        return VoxelLattice(center, affine, shape, header)
    return None


def filter_grid_points(electrodes, mesh, points, material_distribution):
    # TODO locations of points in relation to tissue. needs to be reworked
    raise NotImplementedError(
        "Grid points cannot yet be analysed w.r.t. the brain tissues"
    )


def generate_meshsize_file_from_neuron_grid():
    pass


def generate_signal(settings) -> TimeDomainSignal:
    signal_settings = settings["StimulationSignal"]
    signal_type = signal_settings["Type"]
    if signal_type == "Rectangle":
        signal = RectangleSignal(
            signal_settings["Frequency[Hz]"],
            signal_settings["PulseWidth[us]"],
            signal_settings["InterPulseWidth[us]"],
            signal_settings["CounterPulseWidth[us]"],
        )
    elif signal_type == "Triangle":
        signal = TriangleSignal(
            signal_settings["Frequency[Hz]"],
            signal_settings["PulseWidth[us]"],
            signal_settings["InterPulseWidth[us]"],
            signal_settings["CounterPulseWidth[us]"],
        )
    elif signal_type == "Trapezoid":
        signal = TrapezoidSignal(
            signal_settings["Frequency[Hz]"],
            signal_settings["PulseWidth[us]"],
            signal_settings["InterPulseWidth[us]"],
            signal_settings["CounterPulseWidth[us]"],
            signal_settings["PulseTopWidth[us]"],
        )
    return signal


def prepare_volume_conductor_model(
    settings, model_geometry, conductivity, solver
) -> VolumeConductor:
    _logger.info("Generate volume conductor model")
    order = settings["FEMOrder"]

    mesh_parameters = settings["Mesh"]
    floating_mode = model_geometry.get_floating_mode()
    frequency_domain_signal = prepare_stimulation_signal(settings)
    if floating_mode == "Floating":
        _logger.debug("Floating mode selected")
        return VolumeConductorFloating(
            model_geometry,
            conductivity,
            solver,
            order,
            mesh_parameters,
            frequency_domain_signal,
        )

    elif floating_mode == "FloatingImpedance":
        _logger.debug("FloatingImpedance mode selected")
        return VolumeConductorFloatingImpedance(
            model_geometry,
            conductivity,
            solver,
            order,
            mesh_parameters,
            frequency_domain_signal,
        )
    _logger.debug("Non floating mode selected")
    return VolumeConductorNonFloating(
        model_geometry,
        conductivity,
        solver,
        order,
        mesh_parameters,
        frequency_domain_signal,
    )


def prepare_stimulation_signal(settings) -> FrequencyDomainSignal:
    signal_settings = settings["StimulationSignal"]
    signal_type = signal_settings["Type"]
    current_controlled = signal_settings["CurrentControlled"]
    if signal_type == "Multisine":
        frequencies = signal_settings["ListOfFrequencies"]
        fourier_coefficients = np.ones(len(frequencies))
    else:
        spectrum_mode = signal_settings["SpectrumMode"]
        signal = generate_signal(settings)
        cutoff_frequency = signal_settings["CutoffFrequency"]

        if spectrum_mode == "OctaveBand":
            # TODO add cutoff?!
            frequencies, fourier_coefficients = signal.get_octave_band_spectrum(
                cutoff_frequency
            )
        elif spectrum_mode == "Truncation":
            frequencies, fourier_coefficients = signal.get_truncated_spectrum(
                cutoff_frequency
            )
        else:
            (
                frequencies,
                fourier_coefficients,
            ) = signal.get_frequencies_and_fourier_coefficients(cutoff_frequency)
    frequency_domain_signal = FrequencyDomainSignal(
        frequencies=frequencies,
        amplitudes=fourier_coefficients,
        current_controlled=current_controlled,
    )
    return frequency_domain_signal


def run_volume_conductor_model(settings, volume_conductor):
    """TODO document.


    Notes
    -----
    Run at all frequencies.
    If the mode is multisine, a provided list of frequencies is used.
    """
    _logger.info("Run volume conductor model")
    volume_conductor.output_path = settings["OutputPath"]
    _logger.info(f"Output path set to: {volume_conductor.output_path}")
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

    point_model = generate_neuron_grid(settings)
    template_space = settings["TemplateSpace"]
    vcm_timings = volume_conductor.run_full_analysis(
        compute_impedance,
        export_vtk,
        point_model=point_model,
        template_space=template_space,
        activation_threshold=settings["ActivationThresholdVTA"],
    )
    return vcm_timings


def load_images(settings):
    _logger.info("Load MRI image")
    mri_path = settings["MaterialDistribution"]["MRIPath"]
    _logger.debug(f"Input path: {mri_path}")
    mri_image = MagneticResonanceImage(mri_path)
    dti_image = None
    if settings["MaterialDistribution"]["DiffusionTensorActive"]:
        _logger.info("Load DTI image")
        dti_image = DiffusionTensorImage(settings["MaterialDistribution"]["DTIPath"])
    return mri_image, dti_image
