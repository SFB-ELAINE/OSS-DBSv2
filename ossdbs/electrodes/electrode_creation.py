
from ossdbs.electrodes.electrode import Electrode
from ossdbs.electrodes.abbott_stjude_active_tip_6142_6145 \
    import AbbottStjudeActiveTip6142_6145
from ossdbs.electrodes.abbott_stjude_active_tip_6146_6149 \
    import AbbottStjudeActiveTip6146_6149
from ossdbs.electrodes.abbott_stjude_directed_6172 \
    import AbbottStjudeDirected6172
from ossdbs.electrodes.abbott_stjude_directed_6173 \
    import AbbottStjudeDirected6173
from ossdbs.electrodes.boston_scientific_vercise import BostonScientificVercise
from ossdbs.electrodes.boston_scientific_vercise_directed \
    import BostonScientificVerciseDirected
from ossdbs.electrodes.medtronic_3387 import Medtronic3387
from ossdbs.electrodes.medtronic_3389 import Medtronic3389
from ossdbs.electrodes.medtronic_3391 import Medtronic3391
from ossdbs.electrodes.pins_medical_L301 import PINSMedicalL301
from ossdbs.electrodes.pins_medical_L302 import PINSMedicalL302
from ossdbs.electrodes.pins_medical_L303 import PINSMedicalL303
from ossdbs.electrodes.micro_probes_custom_rodent \
    import MicroProbesCustomRodent
from ossdbs.electrodes.micro_probes_SNEX100 import MicroProbesSNEX_100

from ossdbs.electrodes.custom.abbott_stjude_active_tip_custom \
    import AbbottStjudeActiveCustom
from ossdbs.electrodes.custom.abbott_stjude_directed_custom \
    import AbbottStjudeDirectedCustom
from ossdbs.electrodes.custom.boston_scientific_vercise_custom \
    import BostonScientificVerciseCustom
from ossdbs.electrodes.custom.boston_scientific_vercise_directed_custom \
    import BostonScientificVerciseDirectedCustom
from ossdbs.electrodes.custom.medtronic_custom \
    import MedtronicCustom
from ossdbs.electrodes.custom.micro_probes_custom_rodent_custom \
    import MicroProbesCustomRodentCustom
from ossdbs.electrodes.custom.micro_probes_SNEX100_custom \
    import MicroProbesSNEX_100Custom
from ossdbs.electrodes.custom.pins_medical_custom \
    import PINSMedicalCustom

from dataclasses import dataclass


@dataclass
class ElectrodeParameters:
    name: str = 'BostonScientificVercise'
    position: tuple = (0, 0, 0)
    direction: tuple = (0, 0, 1)
    rotation: float = 0.0


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
    def create(cls, parameters: ElectrodeParameters) -> Electrode:
        """create a list of Electrode objects.

        Parameters
        ----------
        parameters : ElectrodeParameters

        Returns
        -------
        Electrode
            Electrode objects.
        """

        electrode_type = cls.ELECTRODES[parameters.name]
        return electrode_type(direction=parameters.direction,
                              position=parameters.position,
                              rotation=parameters.rotation)
