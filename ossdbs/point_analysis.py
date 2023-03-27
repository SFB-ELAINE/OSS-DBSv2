

from ossdbs.nifti1ImageX import Nifti1Image
from ossdbs.brain_geometry import BrainGeometry
from ossdbs.materials import Material
from ossdbs.output import OutputDirectory

from ossdbs.factories import ActivationModelFactory
from ossdbs.factories import BoundingBoxFactory
from ossdbs.factories import ConductivityFactory
from ossdbs.factories import ContactsFactory
from ossdbs.factories import DielectricModelFactory
from ossdbs.factories import ElectrodesFactory
from ossdbs.factories import MeshFactory
from ossdbs.factories import MaterialDistributionFactory
from ossdbs.factories import SignalFactory
from ossdbs.factories import SolverFactory
from ossdbs.factories import SpectrumFactory
from ossdbs.factories import VolumeConductorFactory

import os
import numpy as np


def point_analysis(input: dict) -> None:

    bounding_box = BoundingBoxFactory.create(input['RegionOfInterest'])
    nifti_image = Nifti1Image(input['MaterialDistribution']['MRIPath'])

    material_distribution = MaterialDistributionFactory(nifti_image,
                                                        bounding_box).create()
    dielectric_model = DielectricModelFactory().create(input['DielectricModel'])

    electrodes = ElectrodesFactory.create(input['Electrodes'])
    capsule_d = input['EncapsulatingLayer']['Thickness[mm]']
    capsule = electrodes.encapsulating_layer(capsule_d)
    max_h = input['EncapsulatingLayer']['MaxMeshSizeHeight']
    capsule.set_max_h(max_h=max_h)
    conductivity = ConductivityFactory(nifti_image,
                                       bounding_box,
                                       dielectric_model,
                                       capsule).create()

    geometry = BrainGeometry(bounding_box, electrodes, capsule)
    mesh = MeshFactory(geometry).create_mesh(input['Mesh'])
    mesh.datatype_complex(input['EQSMode'])

    solver = SolverFactory.create(input['Solver'])

    factory = VolumeConductorFactory(mesh, conductivity, solver)
    volume_conductor = factory.create(input['Floating'])

    signal = SignalFactory.create(input['StimulationSignal'])

    acivation_model = ActivationModelFactory().create(input['ActivationModel'])
    points = acivation_model.coordinates()

    contacts = electrodes.contacts()
    contacts = ContactsFactory(electrodes).create(input['CaseGrounding'])
    mode = SpectrumFactory.create(input['SpectrumMode'],
                                  input['CurrentControled'],
                                  len(contacts.active()))
    output = OutputDirectory(directory=input['OutputPath'])

    result = mode.compute(signal, volume_conductor, points, contacts, output.directory())
    location = define_locations(electrodes, capsule, mesh, points, material_distribution)
    acivation_model.set_location_names(location)
    acivation_model.save(result, os.path.join(output.directory(), 'vta.h5'))

    if input['Mesh']['SaveMesh']:
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
        location[index] = np.any(in_voxel)

    return location
