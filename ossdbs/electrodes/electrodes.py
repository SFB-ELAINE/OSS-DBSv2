from dataclasses import dataclass
from . import AbstractElectrode
from . import AbbottStjudeActiveTip6142_6145
from . import AbbottStjudeActiveTip6146_6149
from . import AbbottStjudeDirected6172
from . import BostonScientificVercise
from . import Medtronic3387, Medtronic3389, Medtronic3391
from . import MicroProbesSNEX_100
from . import PINSMedicalL301
from . import PINSMedicalL302
from . import PINSMedicalL303
from . import MicroProbesCustomRodent


@dataclass
class ElectrodeParameters:
    name: str = 'Rodden'
    translation: tuple = (0., 0., 0.)
    direction: tuple = (0., 0., 0.)
    rotation: float = 0.0


class ElectrodeCreator:

    ELECTRODES = {'AbbottStjudeActiveTip6142_6145':
                  AbbottStjudeActiveTip6142_6145,
                  'AbbottStjudeActiveTip6146_6149':
                  AbbottStjudeActiveTip6146_6149,
                  'AbbottStjudeDirected6172':
                  AbbottStjudeDirected6172,
                  'BostonScientificVercise':
                  BostonScientificVercise,
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

    @classmethod
    def create(cls, parameters: ElectrodeParameters) -> AbstractElectrode:
        trans = parameters.translation
        rot = parameters.rotation
        dir = parameters.direction
        return cls.ELECTRODES[parameters.name](direction=dir,
                                               translation=trans,
                                               rotation=rot)
