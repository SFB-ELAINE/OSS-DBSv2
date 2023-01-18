
from ossdbs.electrodes import AbbottStjudeActiveTip6142_6145
from ossdbs.electrodes import AbbottStjudeActiveTip6146_6149
from ossdbs.electrodes import AbbottStjudeDirected6172
from ossdbs.electrodes import BostonScientificVercise
from ossdbs.electrodes import BostonScientificVerciseDirected
from ossdbs.electrodes import Medtronic3387, Medtronic3389, Medtronic3391
from ossdbs.electrodes import MicroProbesSNEX_100
from ossdbs.electrodes import PINSMedicalL301
from ossdbs.electrodes import PINSMedicalL302
from ossdbs.electrodes import PINSMedicalL303
from ossdbs.electrodes import MicroProbesCustomRodent
from tests.geometry_converter import GeometryConverter
import ngsolve
import netgen


ELECTRODES = {'AbbottStjudeActiveTip6142_6145':
              AbbottStjudeActiveTip6142_6145,
              'AbbottStjudeActiveTip6146_6149':
              AbbottStjudeActiveTip6146_6149,
              'AbbottStjudeDirected6172':
              AbbottStjudeDirected6172,
              'BostonScientificVercise':
              BostonScientificVercise,
              'BostonScientificVerciseDirected':
              BostonScientificVerciseDirected,
              'Medtronic3387':
              Medtronic3387,
              'Medtronic3389':
              Medtronic3389,
              'Medtronic3391':
              Medtronic3391,
              'MicroProbesSNEX_100':
              MicroProbesSNEX_100,
              'PINSMedicalL301':
              PINSMedicalL301,
              'PINSMedicalL302':
              PINSMedicalL302,
              'PINSMedicalL303':
              PINSMedicalL303,
              'MicroProbesCustomRodent':
              MicroProbesCustomRodent
              }


def create_json(electrode_parameters, electrode_type, path):

    electrode_parameters = (0.0, (0, 0, 0), (0, 0, 1))
    electrode_type = AbbottStjudeActiveTip6142_6145
    path = "tests/test_data/AbbottStjudeActiveTip6142_6145_0.json"

    rotation, translation, direction = electrode_parameters
    electrode = electrode_type(rotation, direction, translation)
    GeometryConverter(electrode.generate_geometry()).to_json(path)


def draw_electrode(electrode):

    electrode = AbbottStjudeActiveTip6142_6145()

    geometry = electrode.generate_geometry()
    with ngsolve.TaskManager():
        mesh = ngsolve.Mesh(netgen.occ.OCCGeometry(geometry).GenerateMesh())
    bnd_dict = {"Contact_{}".format(i): float(i) for i in range(1, 9)}
    bndcf = mesh.BoundaryCF(bnd_dict, default=-1)
    ngsolve.Draw(bndcf, mesh, "BND")
