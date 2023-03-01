

from ossdbs.Nifti1Image import Nifti1Image
from ossdbs.brain_geometry import BrainGeometry
from ossdbs.factories.bounding_box_creation import BoundingBoxFactory
from ossdbs.factories.create_conductivity import ConductivityFactory
from ossdbs.factories.electrode_creation import ElectrodesFactory
from ossdbs.factories.mesh_creation import MeshFactory
from ossdbs.factories.points_creation import PointsFactory
from ossdbs.factories.signal_creation import SignalFactory
from ossdbs.factories.solver_creation import SolverFactory
from ossdbs.factories.spectrum_creation import SpectrumFactory
from ossdbs.factories.volume_conductor_creation import VolumeConductorFactory
from ossdbs.point_analysis.spectrum import FullSpectrum


def point_analysis(input: dict) -> None:

    bounding_box = BoundingBoxFactory.create(input['RegionOfInterest'])
    nifti_image = Nifti1Image(input['MaterialDistribution']['MRIPath'])
    electrodes = ElectrodesFactory.create(input['Electrodes'])

    capsule_d = input['EncapsulatingLayer']['Thickness[mm]']
    capsule = electrodes.encapsulating_layer(capsule_d)
    max_h = input['EncapsulatingLayer']['MaxMeshSizeHeight']
    capsule.set_max_h(max_h=max_h)

    conductivity = ConductivityFactory(nifti_image,
                                       bounding_box,
                                       capsule).create()
    conductivity.datatype_complex(input['EQSMode'])

    geometry = BrainGeometry(bounding_box, electrodes, capsule)
    mesh = MeshFactory(geometry).create_mesh(input['Mesh'])
    mesh.datatype_complex(input['EQSMode'])

    solver = SolverFactory.create(input['Solver'])

    factory = VolumeConductorFactory(mesh, conductivity, electrodes, solver)
    volume_conductor = factory.create(input['Floating'])

    signal = SignalFactory.create(input['StimulationSignal'])

    points = PointsFactory.create(input['Points'])

    mode = SpectrumFactory.create(input['SpectrumMode'])
    result = mode.compute(signal, volume_conductor, points)

    categories = PointsFactory.categories(input['Points'])
    result.save_by_categories("test_result.hdf5", categories)
