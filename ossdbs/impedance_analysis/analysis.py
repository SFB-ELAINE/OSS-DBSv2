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

from ossdbs.settings import Settings


def impedance_analysis(settings: dict) -> None:

    settings = Settings(settings).complete_settings()
    region_of_interest = RegionOfInterestFactory.create(
                                                settings['RegionOfInterest'])
    nifti_image = Nifti1Image(settings['MaterialDistribution']['MRIPath'])

    electrodes = ElectrodesFactory.create(settings['Electrodes'])

    capsule_d = settings['EncapsulatingLayer']['Thickness[mm]']
    capsule = electrodes.encapsulating_layer(capsule_d)

    max_h = settings['Contacts']['MaxMeshSizeHeight']
    electrodes.set_max_h(max_h=max_h if max_h else 0.1)

    max_h = settings['EncapsulatingLayer']['MaxMeshSizeHeight']
    capsule.set_max_h(max_h=max_h if max_h else 0.1)

    dielectric_model = DielectricModelFactory().create(
                                                settings['DielectricModel'])
    conductivity = ConductivityFactory(nifti_image,
                                       region_of_interest,
                                       dielectric_model,
                                       capsule).create()

    geometry = BrainGeometry(region_of_interest, electrodes, capsule)
    mesh = MeshFactory(geometry).create(settings['Mesh'])
    mesh.datatype_complex(settings['EQSMode'])

    solver = SolverFactory.create(settings['Solver'])
    factory = VolumeConductorFactory(mesh, conductivity, solver)
#########################################################
    volume_conductor = factory.create(settings['Floating'])
    signal = SignalFactory.create(settings['StimulationSignal'])
    contacts = ContactsFactory(electrodes).create(settings['CaseGrounding'])

    mode = SpectrumFactory.create(settings['SpectrumMode'])
    impedances = mode.compute(signal, volume_conductor, contacts)

    output = OutputDirectory(directory=settings['OutputPath'])
    impedances.save(os.path.join(output.directory(), 'impedances.csv'))

    if settings['Mesh']['SaveMesh']:
        mesh_path = os.path.join(output.directory(), 'mesh')
        mesh.save(mesh_path)
