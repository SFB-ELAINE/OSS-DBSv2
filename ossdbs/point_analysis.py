

from ossdbs.Nifti1Image import Nifti1Image
from ossdbs.brain_geometry import BrainGeometry
from ossdbs.outputX import OutputDirectory

from ossdbs.factories import BoundingBoxFactory
from ossdbs.factories import ConductivityFactory
from ossdbs.factories import ContactsFactory
from ossdbs.factories import DielectricModelFactory
from ossdbs.factories import ElectrodesFactory
from ossdbs.factories import MeshFactory
from ossdbs.factories import PointsFactory
from ossdbs.factories import SignalFactory
from ossdbs.factories import SolverFactory
from ossdbs.factories import SpectrumFactory
from ossdbs.factories import VolumeConductorFactory
from ossdbs.factories import VTAPointsFactory

import os


def point_analysis(input: dict) -> None:

    bounding_box = BoundingBoxFactory.create(input['RegionOfInterest'])
    nifti_image = Nifti1Image(input['MaterialDistribution']['MRIPath'])
    electrodes = ElectrodesFactory.create(input['Electrodes'])
    contacts = electrodes.contacts()
    contacts = ContactsFactory(electrodes).create(input['CaseGrounding'])
    output = OutputDirectory(directory=input['OutputPath'])

    capsule_d = input['EncapsulatingLayer']['Thickness[mm]']
    capsule = electrodes.encapsulating_layer(capsule_d)
    max_h = input['EncapsulatingLayer']['MaxMeshSizeHeight']
    capsule.set_max_h(max_h=max_h)

    dielectric_model = DielectricModelFactory().create(input['DielectricModel'])
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

    # points = PointsFactory.create(input['Points'])

    points = VTAPointsFactory.create(input['VTA'])

    mode = SpectrumFactory.create(input['SpectrumMode'], False, len(contacts.active()))

    result = mode.compute(signal, volume_conductor, points, contacts, output.output_directory())

    # categories = PointsFactory.categories(input['Points'])
    # result.save_by_categories("test_result.hdf5", categories)

    if input['Mesh']['SaveMesh']:
        mesh_path = os.path.join(output.output_directory(), 'mesh')
        mesh.save(mesh_path)
