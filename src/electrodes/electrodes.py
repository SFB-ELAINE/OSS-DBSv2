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
    def create(cls,
               name: str,
               translation: tuple = (0, 0., 0.),
               direction: tuple = (0, 0., 0.),
               rotation: float = 0.0) -> AbstractElectrode:
        return cls.ELECTRODES[name](direction=direction,
                                    translation=translation,
                                    rotation=rotation)
