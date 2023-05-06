from ossdbs.model_geometry import BrainGeometry
from ossdbs.nifti1Image import Nifti1Image
from ossdbs.output_directory import OutputDirectory
from ossdbs.factories import BrainRegionFactory
from ossdbs.factories import ConductivityFactory
from ossdbs.factories import CaseGroundContactFactory
from ossdbs.factories import DielectricModelFactory
from ossdbs.factories import ElectrodesFactory
from ossdbs.factories import MeshFactory
from ossdbs.factories import SolverFactory
from ossdbs.factories import SignalFactory
from ossdbs.factories import VolumeConductorFactory
from ossdbs.settings import Settings
from ossdbs.type_check import TypeChecker
from ossdbs.impedance_analysis.factories import SpectrumFactory

import os

import logging
_logger = logging.getLogger(__name__)


def impedance_analysis(settings: dict) -> None:
    _logger.info("Loading settings from input file")
    settings = Settings(settings).complete_settings()
    TypeChecker.check(settings)

    nifti_image = Nifti1Image(settings['MaterialDistribution']['MRIPath'])
    electrodes = ElectrodesFactory.create(settings['Electrodes'],
                                          settings['EncapsulatingLayer'])
    # TODO set MRI bounding box as default choice
    region_parameters = settings['BrainRegion']
    brain_region = BrainRegionFactory.create(region_parameters)
    dielectric_parameters = settings['DielectricModel']
    dielectric_model = DielectricModelFactory().create(dielectric_parameters)
    encapsulation_material = settings['EncapsulatingLayer']['Material']
    conductivity = ConductivityFactory(nifti_image,
                                       brain_region,
                                       dielectric_model,
                                       electrodes,
                                       encapsulation_material).create()

    geometry = BrainGeometry(brain_region, electrodes)
    mesh = MeshFactory(geometry).create(settings['Mesh'])
    mesh.datatype_complex(settings['EQSMode'])

    solver = SolverFactory.create(settings['Solver'])
    factory = VolumeConductorFactory(mesh, conductivity, solver)
    volume_conductor = factory.create(settings['Floating'])
    signal = SignalFactory.create(settings['StimulationSignal'])

    contacts = electrodes.contacts()
    contacts.append(CaseGroundContactFactory.create(settings['CaseGrounding']))

    # From here on different from point analysis
    spectrum_mode = SpectrumFactory.create(settings['SpectrumMode'])

    impedances = spectrum_mode.compute(signal, volume_conductor, contacts)

    output = OutputDirectory(directory=settings['OutputPath'])
    impedances.save(os.path.join(output.directory(), 'impedances.csv'))

    if settings['Mesh']['SaveMesh']:
        mesh_path = os.path.join(output.directory(), 'mesh')
        mesh.save(mesh_path)
