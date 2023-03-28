

from ossdbs.nifti1ImageX import Nifti1Image
from ossdbs.brain_geometry import BrainGeometry
from ossdbs.materials import Material
from ossdbs.output import OutputDirectory

from ossdbs.factories import RegionOfInterestFactory
from ossdbs.factories import ConductivityFactory
from ossdbs.factories import ContactsFactory
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

    bounding_box = RegionOfInterestFactory.create(settings['RegionOfInterest'])
    nifti_image = Nifti1Image(settings['MaterialDistribution']['MRIPath'])
    dielectric_model = DielectricModelFactory().create(settings['DielectricModel'])

    electrodes = ElectrodesFactory.create(settings['Electrodes'])
    capsule_d = settings['EncapsulatingLayer']['Thickness[mm]']
    capsule = electrodes.encapsulating_layer(capsule_d)

    max_h = settings['Contacts']['MaxMeshSizeHeight']
    electrodes.set_max_h(max_h=max_h if max_h else 0.1)

    max_h = settings['EncapsulatingLayer']['MaxMeshSizeHeight']
    capsule.set_max_h(max_h=max_h if max_h else 0.1)

    conductivity = ConductivityFactory(nifti_image,
                                       bounding_box,
                                       dielectric_model,
                                       capsule).create()

    geometry = BrainGeometry(bounding_box, electrodes, capsule)
    mesh = MeshFactory(geometry).create(settings['Mesh'])
    mesh.datatype_complex(settings['EQSMode'])
    solver = SolverFactory.create(settings['Solver'])
    factory = VolumeConductorFactory(mesh, conductivity, solver)
    volume_conductor = factory.create(settings['Floating'])

    signal = SignalFactory.create(settings['StimulationSignal'])

    contacts = electrodes.contacts()
    contacts = ContactsFactory(electrodes).create(settings['CaseGrounding'])
    mode = SpectrumFactory.create(settings['SpectrumMode'],
                                  settings['CurrentControled'],
                                  len(contacts.active()))
    output = OutputDirectory(directory=settings['OutputPath'])

    point_model = PointModelFactory().create(settings['PointModel'])
    points = point_model.coordinates()
    result = mode.compute(signal, volume_conductor, points, contacts)

    if result.field_solution:
        filename = os.path.join(output.directory(), 'field_solution_130Hz')
        result.field_solution.save(filename=filename)

    material_distribution = MaterialDistributionFactory(nifti_image,
                                                        bounding_box).create()
    location = define_locations(electrodes, capsule, mesh, points, material_distribution)
    point_model.set_location_names(location)
    point_model.save(result, os.path.join(output.directory(),
                                          'points_analysis.h5'))

    if settings['Mesh']['SaveMesh']:
        mesh_path = os.path.join(output.directory(), 'mesh')
        mesh.save(mesh_path)


def define_locations(electrodes, capsule, mesh, points, material_distribution):
    location = np.full(len(points), '')
    in_mesh = mesh.is_included(points=points)
    in_capsule = capsule.is_included(points=points)
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
