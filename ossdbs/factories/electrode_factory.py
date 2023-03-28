

from ossdbs.electrodes.electrode_models import ElectrodeModel
from ossdbs.electrodes.electrode_models import AbbottStjudeActiveTip6142_6145
from ossdbs.electrodes.electrode_models import AbbottStjudeActiveTip6146_6149
from ossdbs.electrodes.electrode_models import AbbottStjudeDirected6172
from ossdbs.electrodes.electrode_models import AbbottStjudeDirected6173
from ossdbs.electrodes.electrode_models import BostonScientificVercise
from ossdbs.electrodes.electrode_models import BostonScientificVerciseDirected
from ossdbs.electrodes.electrode_models import Medtronic3387
from ossdbs.electrodes.electrode_models import Medtronic3389
from ossdbs.electrodes.electrode_models import Medtronic3391
from ossdbs.electrodes.electrode_models import PINSMedicalL301
from ossdbs.electrodes.electrode_models import PINSMedicalL302
from ossdbs.electrodes.electrode_models import PINSMedicalL303
from ossdbs.electrodes.electrode_models import MicroProbesRodentElectrode
from ossdbs.electrodes.electrode_models import MicroProbesSNEX100


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
                  'MicroProbesSNEX100':
                  MicroProbesSNEX100,
                  'PINSMedicalL301':
                  PINSMedicalL301,
                  'PINSMedicalL302':
                  PINSMedicalL302,
                  'PINSMedicalL303':
                  PINSMedicalL303,
                  'MicroProbesRodentElectrode':
                  MicroProbesRodentElectrode,
                  }

    @classmethod
    def create(cls,
               name: str,
               direction: tuple,
               position: tuple,
               rotation: float
               ) -> ElectrodeModel:
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
        return electrode_type(direction=direction,
                              position=position,
                              rotation=rotation)
