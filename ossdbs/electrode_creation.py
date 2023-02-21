
from ossdbs.electrodes import Electrode
from ossdbs.electrodes import AbbottStjudeActiveTip6142_6145
from ossdbs.electrodes import AbbottStjudeActiveTip6146_6149
from ossdbs.electrodes import AbbottStjudeDirected6172
from ossdbs.electrodes import AbbottStjudeDirected6173
from ossdbs.electrodes import BostonScientificVercise
from ossdbs.electrodes import BostonScientificVerciseDirected
from ossdbs.electrodes import Medtronic3387
from ossdbs.electrodes import Medtronic3389
from ossdbs.electrodes import Medtronic3391
from ossdbs.electrodes import PINSMedicalL301
from ossdbs.electrodes import PINSMedicalL302
from ossdbs.electrodes import PINSMedicalL303
from ossdbs.electrodes import MicroProbesCustomRodent
from ossdbs.electrodes import MicroProbesSNEX_100

from ossdbs.electrodes import AbbottStjudeActiveCustom
from ossdbs.electrodes import AbbottStjudeDirectedCustom
from ossdbs.electrodes import BostonScientificVerciseCustom
from ossdbs.electrodes import BostonScientificVerciseDirectedCustom
from ossdbs.electrodes import MedtronicCustom
from ossdbs.electrodes import MicroProbesCustomRodentCustom
from ossdbs.electrodes import MicroProbesSNEX_100Custom
from ossdbs.electrodes import PINSMedicalCustom


class ElectrodeFactory:
    """Creates a list of Electrode objects."""

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
                  MicroProbesCustomRodent,
                  'AbbottStjudeActiveCustom':
                  AbbottStjudeActiveCustom,
                  'AbbottStjudeDirectedCustom':
                  AbbottStjudeDirectedCustom,
                  'BostonScientificVerciseCustom':
                  BostonScientificVerciseCustom,
                  'BostonScientificVerciseDirectedCustom':
                  BostonScientificVerciseDirectedCustom,
                  'MedtronicCustom':
                  MedtronicCustom,
                  'MicroProbesCustomRodentCustom':
                  MicroProbesCustomRodentCustom,
                  'MicroProbesSNEX_100Custom':
                  MicroProbesSNEX_100Custom,
                  'PINSMedicalCustom':
                  PINSMedicalCustom
                  }

    @classmethod
    def create(cls,
               name: str,
               direction: tuple,
               position: tuple,
               rotation: float) -> Electrode:
        """create a list of Electrode objects.

        Parameters
        ----------
        parameters : ElectrodeParameters

        Returns
        -------
        Electrode
            Electrode objects.
        """

        electrode_type = cls.ELECTRODES[name]
        return electrode_type(direction, position, rotation)
