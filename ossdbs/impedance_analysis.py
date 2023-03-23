import os
from ossdbs.Nifti1Image import Nifti1Image
from ossdbs.brain_geometry import BrainGeometry
from ossdbs.output import OutputDirectory

from ossdbs.factories import BoundingBoxFactory
from ossdbs.factories import ConductivityFactory
from ossdbs.factories import ContactsFactory
from ossdbs.factories import DielectricModelFactory
from ossdbs.factories import ElectrodesFactory
from ossdbs.factories import MeshFactory
from ossdbs.factories import SolverFactory
from ossdbs.factories import SignalFactory
from ossdbs.factories import VolumeConductorFactory
from ossdbs.factories import SpectrumImpedanceFactory


def impedance_analysis(input: dict) -> None:

    bounding_box = BoundingBoxFactory.create(input['RegionOfInterest'])
    nifti_image = Nifti1Image(input['MaterialDistribution']['MRIPath'])
    electrodes = ElectrodesFactory.create(input['Electrodes'])
    contacts = ContactsFactory(electrodes).create(input['CaseGrounding'])
    output = OutputDirectory(directory=input['OutputPath'])

    capsule_d = input['EncapsulatingLayer']['Thickness[mm]']
    capsule = electrodes.encapsulating_layer(capsule_d)
    electrodes.set_max_h(max_h=input['Contacts']['MaxMeshSizeHeight'])
    capsule.set_max_h(max_h=input['EncapsulatingLayer']['MaxMeshSizeHeight'])

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

    mode = SpectrumImpedanceFactory.create(input['SpectrumMode'])
    impedances = mode.compute(signal, volume_conductor, contacts)
    impedances.save(os.path.join(output.output_directory(), 'impedances.csv'))

    if input['Mesh']['SaveMesh']:
        mesh_path = os.path.join(output.output_directory(), 'mesh')
        mesh.save(mesh_path)
