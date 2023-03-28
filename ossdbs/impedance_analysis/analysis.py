from ossdbs.brain_geometry import BrainGeometry
from ossdbs.nifti1Image import Nifti1Image
from ossdbs.output import OutputDirectory
from ossdbs.factories import RegionOfInterestFactory
from ossdbs.factories import ConductivityFactory
from ossdbs.factories import ContactsFactory
from ossdbs.factories import DielectricModelFactory
from ossdbs.factories import ElectrodesFactory
from ossdbs.factories import MeshFactory
from ossdbs.factories import SolverFactory
from ossdbs.factories import SignalFactory
from ossdbs.factories import VolumeConductorFactory
from ossdbs.impedance_analysis.factories import SpectrumFactory

import os


def impedance_analysis(setting: dict) -> None:

    bounding_box = RegionOfInterestFactory.create(setting['RegionOfInterest'])
    nifti_image = Nifti1Image(setting['MaterialDistribution']['MRIPath'])
    electrodes = ElectrodesFactory.create(setting['Electrodes'])
    contacts = ContactsFactory(electrodes).create(setting['CaseGrounding'])
    output = OutputDirectory(directory=setting['OutputPath'])

    capsule_d = setting['EncapsulatingLayer']['Thickness[mm]']
    capsule = electrodes.encapsulating_layer(capsule_d)

    max_h = setting['Contacts']['MaxMeshSizeHeight']
    electrodes.set_max_h(max_h=max_h if max_h else 0.1)

    max_h = setting['EncapsulatingLayer']['MaxMeshSizeHeight']
    capsule.set_max_h(max_h=max_h if max_h else 0.1)

    dielectric_model = DielectricModelFactory().create(setting['DielectricModel'])
    conductivity = ConductivityFactory(nifti_image,
                                       bounding_box,
                                       dielectric_model,
                                       capsule).create()

    geometry = BrainGeometry(bounding_box, electrodes, capsule)
    mesh = MeshFactory(geometry).create(setting['Mesh'])
    mesh.datatype_complex(setting['EQSMode'])

    solver = SolverFactory.create(setting['Solver'])

    factory = VolumeConductorFactory(mesh, conductivity, solver)
    volume_conductor = factory.create(setting['Floating'])

    signal = SignalFactory.create(setting['StimulationSignal'])

    mode = SpectrumFactory.create(setting['SpectrumMode'])
    impedances = mode.compute(signal, volume_conductor, contacts)
    impedances.save(os.path.join(output.directory(), 'impedances.csv'))

    if setting['Mesh']['SaveMesh']:
        mesh_path = os.path.join(output.directory(), 'mesh')
        mesh.save(mesh_path)
