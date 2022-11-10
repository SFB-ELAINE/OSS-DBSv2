from dataclasses import dataclass
from src.electrodes import AbstractElectrode
from src.electrodes import AbbottStjudeActiveTip6142_6145
from src.electrodes import AbbottStjudeActiveTip6146_6149
from src.electrodes import AbbottStjudeDirected6172
from src.electrodes import BostonScientificVercise
from src.electrodes import Medtronic3387, Medtronic3389, Medtronic3391
from src.electrodes import MicroProbesSNEX_100
from src.electrodes import PINSMedicalL301
from src.electrodes import PINSMedicalL302
from src.electrodes import PINSMedicalL303
from src.electrodes import MicroProbesCustomRodent


@dataclass
class ElectrodeParameters:
    name: str = 'Rodden'
    translation: tuple = (0., 0., 0.)
    direction: tuple = (0., 0., 0.)
    rotation: float = 0.0
    contact_values: list = None


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
        contacts = parameters.contact_values
        return cls.ELECTRODES[parameters.name](direction=dir,
                                               translation=trans,
                                               rotation=rot,
                                               boundaries=contacts)
