

from ossdbs.nifti1Image import Nifti1Image
from ossdbs.brain_geometry import BrainGeometry
from ossdbs.materials import Material
from ossdbs.output_directory import OutputDirectory
from ossdbs.settings import Settings
from ossdbs.type_check import TypeChecker
from ossdbs.factories import RegionOfInterestFactory
from ossdbs.factories import ConductivityFactory
from ossdbs.factories import CaseGroundContactFactory
from ossdbs.factories import DielectricModelFactory
from ossdbs.factories import ElectrodesFactory
from ossdbs.factories import MeshFactory
from ossdbs.factories import SignalFactory
from ossdbs.factories import SolverFactory
from ossdbs.factories import VolumeConductorFactory

from ossdbs.point_analysis.factories.pointmodel_factory import PointModelFactory
from ossdbs.point_analysis.factories import MaterialDistributionFactory
from ossdbs.point_analysis.factories import SpectrumFactory

import os
import numpy as np


def point_analysis(settings: dict) -> None:

    settings = Settings(settings).complete_settings()
    TypeChecker.check(settings)

    nifti_image = Nifti1Image(settings['MaterialDistribution']['MRIPath'])
    electrodes = ElectrodesFactory.create(settings['Electrodes'],
                                          settings['EncapsulatingLayer'])
    region_parameters = settings['RegionOfInterest']
    region_of_interest = RegionOfInterestFactory.create(region_parameters)
    dielectric_parameters = settings['DielectricModel']
    dielectric_model = DielectricModelFactory().create(dielectric_parameters)
    encapsulation_material = settings['EncapsulatingLayer']['Material']
    conductivity = ConductivityFactory(nifti_image,
                                       region_of_interest,
                                       dielectric_model,
                                       electrodes,
                                       encapsulation_material).create()

    geometry = BrainGeometry(region_of_interest, electrodes)
    mesh = MeshFactory(geometry).create(settings['Mesh'])
    mesh.datatype_complex(settings['EQSMode'])

    solver = SolverFactory.create(settings['Solver'])
    factory = VolumeConductorFactory(mesh, conductivity, solver)
    volume_conductor = factory.create(settings['Floating'])
    signal = SignalFactory.create(settings['StimulationSignal'])
    contacts = electrodes.contacts()
    contacts.append(CaseGroundContactFactory.create(settings['CaseGrounding']))
    spectrum_mode = SpectrumFactory.create(settings['SpectrumMode'],
                                           settings['CurrentControled'],
                                           len(contacts.active()))

    point_model = PointModelFactory().create(settings['PointModel'])
    points = point_model.coordinates()

    material_distribution = MaterialDistributionFactory(nifti_image,
                                                        region_of_interest
                                                        ).create()

    location = define_locations(electrodes, mesh, points, material_distribution)
    point_model.set_location_names(location)

    points = point_model.coordinates()
    result = spectrum_mode.compute(signal, volume_conductor, points, contacts)

    output = OutputDirectory(directory=settings['OutputPath'])
    if result.field_solution:
        filename = os.path.join(output.directory(), 'field_solution_130Hz')
        result.field_solution.save(filename=filename)

    point_model.save(result,
                     os.path.join(output.directory(), 'points_analysis.h5'))

    if settings['Mesh']['SaveMesh']:
        mesh_path = os.path.join(output.directory(), 'mesh')
        mesh.save(mesh_path)


def create_point_model(settings, mesh):
    nifti_image = Nifti1Image(settings['MaterialDistribution']['MRIPath'])
    region_parameters = settings['RegionOfInterest']
    region_of_interest = RegionOfInterestFactory.create(region_parameters)
    material_distribution = MaterialDistributionFactory(nifti_image,
                                                        region_of_interest).create()
    point_model = PointModelFactory().create(settings['PointModel'])
    points = point_model.coordinates()

    electrodes = ElectrodesFactory.create(settings['Electrodes'])
    capsule_d = settings['EncapsulatingLayer']['Thickness[mm]']
    capsule = electrodes.encapsulation(capsule_d)
    max_h = settings['EncapsulatingLayer']['MaxMeshSizeHeight']
    capsule.set_max_h(max_h=max_h if max_h else 0.1)

    location = define_locations(electrodes, capsule, mesh, points, material_distribution)
    point_model.set_location_names(location)
    return point_model


def define_locations(electrodes, mesh, points, material_distribution):
    location = np.full(len(points), '')
    in_mesh = mesh.is_included(points=points)
    caplsulation = electrodes.encapsulation()
    in_capsule = caplsulation.is_included(points=points)
    in_electrode = electrodes.is_included(points=points)
    is_tissue = np.logical_and(in_mesh, np.invert(in_capsule))
    is_outside = np.invert(np.logical_or(in_mesh, in_electrode))

    location[is_tissue] = 'Tissue'
    location[is_outside] = 'Outside'
    location[in_capsule] = 'Encapsulation'
    location[in_electrode] = 'Electrode'

    starts, ends = material_distribution.where(Material.CSF)
    for index, point in enumerate(points):
        in_voxel = np.logical_and(np.all(starts <= point, axis=1),
                                  np.all(point < ends, axis=1))
        if np.any(in_voxel):
            location[index] = 'CSF'

    return location
