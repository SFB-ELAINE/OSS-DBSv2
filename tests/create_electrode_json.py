
from ossdbs.electrode_models import AbbottStjudeActiveTip6142_6145
from ossdbs.electrode_models import AbbottStjudeActiveTip6146_6149
from ossdbs.electrode_models import AbbottStjudeDirected6172
from ossdbs.electrode_models import AbbottStjudeDirected6173
from ossdbs.electrode_models import BostonScientificVercise
from ossdbs.electrode_models import BostonScientificVerciseDirected
from ossdbs.electrode_models import Medtronic3387, Medtronic3389, Medtronic3391
from ossdbs.electrode_models import MicroProbesSNEX_100
from ossdbs.electrode_models import PINSMedicalL301
from ossdbs.electrode_models import PINSMedicalL302
from ossdbs.electrode_models import PINSMedicalL303
from ossdbs.electrode_models import MicroProbesCustomRodent
#from tests.geometry_converter import GeometryConverter
import ngsolve
import netgen


ELECTRODES = {'AbbottStjudeActiveTip6142_6145':
              AbbottStjudeActiveTip6142_6145,
              'AbbottStjudeActiveTip6146_6149':
              AbbottStjudeActiveTip6146_6149,
              'AbbottStjudeDirected6172':
              AbbottStjudeDirected6172,
              'AbbottStjudeDirected6173':
              AbbottStjudeDirected6173,
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
    # GeometryConverter(electrode.generate_geometry()).to_json(path)


def draw_electrode():

    electrode = BostonScientificVerciseDirected(direction=(1,1,0))

    electrode_geo = electrode.geometry()
    capsule = electrode.capsule_geometry(0.1, max_h=1)

    box = netgen.occ.Box((-5, -5, -5,), (5, 5, 5))
    cut = capsule - box
    geometry = netgen.occ.Glue([box - capsule, capsule - electrode_geo - cut]) 

    print(set([edge.name for edge in geometry.edges]))
    geometry = capsule - electrode_geo
   # geometry = electrode_geo

    with ngsolve.TaskManager():
        mesh = ngsolve.Mesh(netgen.occ.OCCGeometry(geometry).GenerateMesh())

    print(set(mesh.GetBoundaries()))
    bnd_dict = {"Contact_{}".format(i): i for i in range(1, 9)}
    bnd_dict.update({"Capsule": 9})
    bndcf = mesh.BoundaryCF(bnd_dict, default=-1)
    ngsolve.Draw(bndcf, mesh, "BND")


if __name__ == '__main__':
    draw_electrode()
