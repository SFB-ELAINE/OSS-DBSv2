

from ossdbs.electrodes.electrode_models import ElectrodeModel
from ossdbs.electrodes.electrode_models import AbbottStJudeActiveTip6142_6145
from ossdbs.electrodes.electrode_models import AbbottStJudeActiveTip6146_6149
from ossdbs.electrodes.electrode_models import AbbottStJudeDirected6172
from ossdbs.electrodes.electrode_models import AbbottStJudeDirected6173
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
    """Creates an Electrode object.

    See also
    --------

    :class:`ossdbs.electrodes.ElectrodeModel`

    """

    ELECTRODES = {'AbbottStJudeActiveTip6142_6145':
                  AbbottStJudeActiveTip6142_6145,
                  'AbbottStJudeActiveTip6146_6149':
                  AbbottStJudeActiveTip6146_6149,
                  'AbbottStJudeDirected6172':
                  AbbottStJudeDirected6172,
                  'AbbottStJudeDirected6173':
                  AbbottStJudeDirected6173,
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
        """Creates an Electrode object.

        Parameters
        ----------

        name: str
            Name of the electrode type

        Returns
        -------

        :class:`ossdbs.electrodes.ElectrodeModel`

        """

        electrode_type = cls.ELECTRODES[name]
        return electrode_type(direction=direction,
                              position=position,
                              rotation=rotation)
