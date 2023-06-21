from ossdbs.electrodes import (ELECTRODES,
                               ELECTRODE_MODELS,
                               custom_parameters)
from ossdbs.dielectric_model import (create_cc4_model,
                                     create_constant_model,
                                     DIELECTRIC_MODELS)
from ossdbs.model_geometry import (BoundingBox,
                                   ModelGeometry,
                                   BrainGeometry)
from ossdbs.fem import (SOLVERS,
                        PRECONDITIONERS,
                        Mesh
                        )
from ossdbs.stimulation_signals import (TimeDomainSignal,
                                        FrequencyDomainSignal,
                                        RectangleSignal,
                                        TrapezoidSignal,
                                        TriangleSignal)
from ossdbs.point_analysis import (PointModel,
                                   Pathway,
                                   Lattice)
from ossdbs.fem import (VolumeConductor,
                        VolumeConductorNonFloating,
                        VolumeConductorFloating,
                        VolumeConductorFloatingImpedance)
import numpy as np
import logging

_logger = logging.getLogger(__name__)


def create_bounding_box(box_parameters: dict) -> BoundingBox:
    input_s = box_parameters['Dimension']
    input_c = box_parameters['Center']
    shape = (input_s['x[mm]'], input_s['y[mm]'], input_s['z[mm]'])
    center = (input_c['x[mm]'], input_c['y[mm]'], input_c['z[mm]'])
    start = center - np.divide(shape, 2)
    end = start + shape
    return BoundingBox(tuple(start), tuple(end))


def generate_electrodes(settings: dict):
    """Generate OCC electrode models

    Notes
    -----

    TODO custom electrodes
    """
    _logger.info("Generate electrode geometries")
    electrodes = []
    for electrode_parameters in settings["Electrodes"]:
        name = electrode_parameters["Name"]
        direction = (electrode_parameters["Direction"]["x[mm]"],
                     electrode_parameters["Direction"]["y[mm]"],
                     electrode_parameters["Direction"]["z[mm]"])
        rotation = electrode_parameters["Rotation[Degrees]"]
        position = (electrode_parameters["TipPosition"]["x[mm]"],
                    electrode_parameters["TipPosition"]["y[mm]"],
                    electrode_parameters["TipPosition"]["z[mm]"])

        if "Model" in name:
            electrode_model = ELECTRODE_MODELS[name]
            custom_electrode_parameters = custom_parameters(electrode_parameters["CustomParameters"])
            electrode = electrode_model(parameters=custom_electrode_parameters,
                                        direction=direction,
                                        position=position,
                                        rotation=rotation)

        else:
            electrode_type = ELECTRODES[name]
            electrode = electrode_type(direction=direction,
                                       position=position,
                                       rotation=rotation)

        if "EncapsulationLayer" in electrode_parameters:
            electrode.encapsulation_thickness = electrode_parameters["EncapsulationLayer"]["Thickness[mm]"]
        electrodes.append(electrode)
    return electrodes


def prepare_dielectric_properties(settings: dict):
    _logger.info("Prepare dielectric model")
    dielectric_parameters = settings["DielectricModel"]
    if 'Custom' in dielectric_parameters['Type']:
        model_parameters = dielectric_parameters['CustomParameters']
        if 'ColeCole4' in dielectric_parameters['Type']:
            return create_cc4_model(model_parameters)
        if 'Constant' in dielectric_parameters['Type']:
            return create_constant_model()

    return DIELECTRIC_MODELS[dielectric_parameters['Type']]


def generate_brain_model(settings):
    """Generate OCC brain model

    Notes
    -----

    TODO type checking

    """
    brain_region_parameters = settings['BrainRegion']
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
        if "Contacts" in new_parameters:
            for contact_info in new_parameters["Contacts"]:
                contact_idx = offset + contact_info["Contact_ID"]
                # contacts are zero-indexed in the model_geometry
                model_geometry.update_contact(contact_idx - 1, contact_info)
            offset += model_geometry.electrodes[idx].n_contacts
        if "EncapsulationLayer" in new_parameters:
            # encapsulation layer is one-indexed in the model_geometry
            encap_idx = model_geometry.get_encapsulation_layer_index("EncapsulationLayer_{}".format(idx + 1))
            if encap_idx != -1:
                _logger.info("Updating encapsulation layer properties")
                model_geometry.update_encapsulation_layer(idx, new_parameters["EncapsulationLayer"])
    if "Surfaces" in settings:
        for surface in settings["Surfaces"]:
            idx = model_geometry.get_contact_index(surface["Name"])
            if idx == -1:
                raise ValueError("Surface {} not part of the geometry".format(surface["Name"]))
            model_geometry.update_contact(idx, surface)


def set_custom_mesh_sizes(settings, model_geometry):
    model_geometry.set_volume_mesh_sizes(settings["VolumeMeshSizes"])


def generate_mesh(settings):
    """

    Notes
    -----

    Attention! This mesh is not yet curved!
    """
    model_geometry = generate_model_geometry(settings)
    set_contact_and_encapsulation_layer_properties(settings, model_geometry)
    if "VolumeMeshSizes" in settings:
        set_custom_mesh_sizes(settings, model_geometry)

    mesh_settings = settings["Mesh"]
    mesh_order = 1

    mesh = Mesh(model_geometry.geometry, mesh_order)
    if mesh_settings["LoadMesh"]:
        mesh.load_mesh(mesh_settings["LoadPath"])
        return mesh

    if 'MeshingHypothesis' in mesh_settings:
        mesh_hypothesis = mesh_settings['MeshingHypothesis']
    else:
        mesh_hypothesis = {"Type": "Default"}
    mesh.generate_mesh(mesh_hypothesis)
    if mesh_settings["SaveMesh"]:
        mesh.save(mesh_settings["SavePath"])
    return mesh


def prepare_solver(settings):
    _logger.info("Preparing solver")
    parameters = settings["Solver"]
    solver_type = parameters['Type']
    solver = SOLVERS[solver_type]
    preconditioner = PRECONDITIONERS[parameters['Preconditioner']]

    return solver(precond_par=preconditioner,
                  printrates=parameters['PrintRates'],
                  maxsteps=parameters['MaximumSteps'],
                  precision=parameters['Precision'])


def generate_neuron_grid(settings: dict) -> PointModel:
    if settings['Pathway']['Active']:
        file_name = settings['Pathway']['FileName']
        return Pathway(file_name)

    shape_par = settings['Lattice']['Shape']
    shape = shape_par['x'], shape_par['y'], shape_par['z']
    center_par = settings['Lattice']['Center']
    center = center_par['x[mm]'], center_par['y[mm]'], center_par['z[mm]']
    dir_par = settings['Lattice']['Direction']
    direction = dir_par['x[mm]'], dir_par['y[mm]'], dir_par['z[mm]']
    distance = settings['Lattice']['PointDistance[mm]']

    return Lattice(shape=shape,
                   center=center,
                   distance=distance,
                   direction=direction)


def filter_grid_points(electrodes, mesh, points, material_distribution):
    # TODO locations of points in relation to tissue. needs to be reworked
    raise NotImplementedError("Grid points cannot yet be analysed w.r.t. the brain tissues")


def generate_meshsize_file_from_neuron_grid():
    pass


def generate_signal(settings) -> TimeDomainSignal:
    signal_settings = settings['StimulationSignal']
    signal_type = signal_settings["Type"]
    if signal_type == "Rectangle":
        signal = RectangleSignal(
            signal_settings['Frequency[Hz]'],
            signal_settings['PulseWidth[us]'],
            signal_settings['InterPulseWidth[us]'],
            signal_settings['CounterPulseWidth[us]'])
    elif signal_type == "Triangle":
        signal = TriangleSignal(
            signal_settings['Frequency[Hz]'],
            signal_settings['PulseWidth[us]'],
            signal_settings['InterPulseWidth[us]'],
            signal_settings['CounterPulseWidth[us]'])
    elif signal_type == "Trapezoid":
        signal = TrapezoidSignal(
            signal_settings['Frequency[Hz]'],
            signal_settings['PulseWidth[us]'],
            signal_settings['InterPulseWidth[us]'],
            signal_settings['CounterPulseWidth[us]'],
            signal_settings['PulseTopWidth[us]'])
    return signal


def prepare_volume_conductor_model(settings, model_geometry, conductivity, solver) -> VolumeConductor:
    _logger.info("Generate volume conductor model")
    order = settings["FEMOrder"]

    mesh_parameters = settings["Mesh"]
    floating_mode = model_geometry.get_floating_mode()
    frequency_domain_signal = prepare_stimulation_signal(settings)
    if floating_mode == "Floating":
        return VolumeConductorFloating(model_geometry,
                                       conductivity,
                                       solver,
                                       order,
                                       mesh_parameters,
                                       frequency_domain_signal)

    elif floating_mode == "FloatingImpedance":
        return VolumeConductorFloatingImpedance(model_geometry,
                                                conductivity,
                                                solver,
                                                order,
                                                mesh_parameters,
                                                frequency_domain_signal)

    return VolumeConductorNonFloating(model_geometry,
                                      conductivity,
                                      solver,
                                      order,
                                      mesh_parameters,
                                      frequency_domain_signal)


def prepare_stimulation_signal(settings) -> FrequencyDomainSignal:
    signal_settings = settings['StimulationSignal']
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
            frequencies, fourier_coefficients = signal.get_octave_band_frequencies(cutoff_frequency)
        elif spectrum_mode == "Truncation":
            frequencies, fourier_coefficients = signal.get_truncated_spectrum(cutoff_frequency)
        else:
            frequencies, fourier_coefficients = signal.frequencies, signal.fourier_coefficients
    frequency_domain_signal = FrequencyDomainSignal(frequencies=frequencies,
                                                    amplitudes=fourier_coefficients,
                                                    current_controlled=current_controlled)
    return frequency_domain_signal