
from ossdbs.electrodes import Electrode
from ossdbs.electrodes import AbbottStjudeActiveTip6142_6145
from ossdbs.electrodes import AbbottStjudeActiveTip6146_6149
from ossdbs.electrodes import AbbottStjudeDirected6172
from ossdbs.electrodes import BostonScientificVercise
from ossdbs.electrodes import Medtronic3387, Medtronic3389, Medtronic3391
from ossdbs.electrodes import MicroProbesSNEX_100
from ossdbs.electrodes import PINSMedicalL301
from ossdbs.electrodes import PINSMedicalL302
from ossdbs.electrodes import PINSMedicalL303
from ossdbs.electrodes import MicroProbesCustomRodent


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
               rotation: float = 0.0) -> Electrode:
        return cls.ELECTRODES[name](direction=direction,
                                    translation=translation,
                                    rotation=rotation)
