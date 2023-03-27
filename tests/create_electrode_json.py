
from ossdbs.electrodes.contacts import Contacts
from ossdbs.electrodes.electrode_models import AbbottStjudeActiveTip6142_6145
from ossdbs.electrodes.electrode_models import AbbottStjudeActiveTip6146_6149
from ossdbs.electrodes.electrode_models import AbbottStjudeDirected6172
from ossdbs.electrodes.electrode_models import AbbottStjudeDirected6173
from ossdbs.electrodes.electrode_models import BostonScientificVercise
from ossdbs.electrodes.electrode_models import BostonScientificVerciseDirected
from ossdbs.electrodes.electrode_models import Medtronic3387, Medtronic3389, Medtronic3391
from ossdbs.electrodes.electrode_models import MicroProbesSNEX_100
from ossdbs.electrodes.electrode_models import PINSMedicalL301
from ossdbs.electrodes.electrode_models import PINSMedicalL302
from ossdbs.electrodes.electrode_models import PINSMedicalL303
from ossdbs.electrodes.electrode_models import MicroProbesCustomRodent
from tests.geometry_converter import GeometryConverter
import ngsolve
import netgen
from ossdbs.electrodes import Electrodes


class TestFileCreator:

    ELECTRODES = {
                #   'AbbottStjudeActiveTip6142_6145':
                #   AbbottStjudeActiveTip6142_6145,
                #   'AbbottStjudeActiveTip6146_6149':
                #   AbbottStjudeActiveTip6146_6149,
                #   'AbbottStjudeDirected6172':
                #   AbbottStjudeDirected6172,
                #   'AbbottStjudeDirected6173':
                #   AbbottStjudeDirected6173,
                #   'BostonScientificVercise':
                #   BostonScientificVercise,
                  'BostonScientificVerciseDirected':
                  BostonScientificVerciseDirected,
                #   'Medtronic3387':
                #   Medtronic3387,
                #   'Medtronic3389':
                #   Medtronic3389,
                #   'Medtronic3391':
                #   Medtronic3391,
                #   'MicroProbesSNEX_100':
                #   MicroProbesSNEX_100,
                #   'PINSMedicalL301':
                #   PINSMedicalL301,
                #   'PINSMedicalL302':
                #   PINSMedicalL302,
                #   'PINSMedicalL303':
                #   PINSMedicalL303,
                #   'MicroProbesCustomRodent':
                #   MicroProbesCustomRodent
                  }

    TESTDATA = [
        # electrode_parameters (Rotation, Position, Direction), file_path
        ((0.0, (0, 0, 0), (0, 0, 1))),
        ((0.0, (0, 0, 0), (0, 0, 0))),
        ((30.0, (0, 0, 0), (0, 0, 1))),
        ((0.0, (1, -2, 3), (0, 0, 1))),
        ((0.0, (1, -2, 3), (0, 0, 0))),
        ((30.0, (1, -2, 3), (0, 0, 0))),
        ((0.0, (1, -2, 3), (2.0, 0, 1.0))),
        ((0.0, (1, -2, 3), (2.0/3.0, 0, 1.0/3.0))),
        ]

    TESTDATA_CAPSULE = [
        # electrode_parameters (Thickness, Position, Direction), file_path
        ((1.0, (0, 0, 0), (0, 0, 1))),
        ((2.0, (0, 0, 0), (0, 0, 1))),
        ((1.0, (1, -2, 3), (0, 0, 1))),
        ]

    def create_json(self):

        for key, electrode_type in self.ELECTRODES.items():
            for index, parameters in enumerate(self.TESTDATA):
                path = key + '_' + str(index) + '.json'
                rotation, translation, direction = parameters
                electrode = electrode_type(rotation, direction, translation)
                GeometryConverter(electrode.geometry()).to_json(path)

        for key, electrode_type in self.ELECTRODES.items():
            for index, parameters in enumerate(self.TESTDATA_CAPSULE):
                path = key + '_' + 'Capsule_' + str(index) + '.json'
                thickness, translation, direction = parameters
                electrode = electrode_type(0.0, direction, translation)
                GeometryConverter(electrode.capsule_geometry(thickness)
                                  ).to_json(path)


def create_json(electrode_parameters, electrode_type, path):

    electrode_parameters = (0.0, (0, 0, 0), (0, 0, 1))
    electrode_type = AbbottStjudeActiveTip6142_6145
    path = "test_data/AbbottStjudeActiveTip6142_6145_0.json"

    rotation, translation, direction = electrode_parameters
    electrode = electrode_type(rotation, direction, translation)
    GeometryConverter(electrode.geometry()).to_json(path)


def draw_electrode():

    electrode = BostonScientificVerciseDirected(direction=(0, 0, 1))

    electrodes = Electrodes([electrode], Contacts([]))
    electrode_geo = electrodes.geometry()
    capsule = electrodes.encapsulating_layer(0.1).geometry()

    box = netgen.occ.Box((-10, -10, -10,), (10, 10, 10))
    cut_capsule = capsule - box
    cut_capsule.bc('Capsule')
    cut_electrode = electrode_geo - box
    # geometry = netgen.occ.Glue([box - capsule - electrode_geo, capsule - cut]) 

    capsule_part = capsule - electrode_geo - cut_capsule

    part_1 = box - (capsule - electrode_geo) - electrode_geo
    geometry = netgen.occ.Glue([part_1, capsule_part])


    print(set([edge.name for edge in geometry.edges]))

    with ngsolve.TaskManager():
        mesh = ngsolve.Mesh(netgen.occ.OCCGeometry(geometry).GenerateMesh())

    print(set(mesh.GetBoundaries()))
    bnd_dict = {"Contact_{}".format(i): i for i in range(1, 9)}
    bnd_dict.update({"Capsule": 9})
    bndcf = mesh.BoundaryCF(bnd_dict, default=-1)
    ngsolve.Draw(bndcf, mesh, "BND")


if __name__ == '__main__':
    draw_electrode()
    # TestFileCreator().create_json()
