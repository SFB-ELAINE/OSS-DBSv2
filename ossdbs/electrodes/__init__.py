from .electrode_model_template import ElectrodeModel
from .abbott_stjude import (AbbottStJudeParameters,
                            AbbottStJudeActiveTipModel,
                            AbbottStJudeDirectedModel)
from .boston_scientific_vercise import (BostonScientificVerciseParameters,
                                        BostonScientificVerciseModel,
                                        BostonScientificVerciseDirectedModel)
from .medtronic import (MedtronicParameters,
                        MedtronicModel)
from .micro_probes import (MicroProbesSNEX100Parameters,
                           MicroProbesSNEX100Model,
                           MicroProbesRodentElectrodeParameters,
                           MicroProbesRodentElectrodeModel)
from .pins_medical import (PINSMedicalParameters,
                           PINSMedicalModel)
from .defaults import (default_electrode_parameters,
                       AbbottStJudeActiveTip6142_6145,
                       AbbottStJudeActiveTip6146_6149,
                       AbbottStJudeDirected6172,
                       AbbottStJudeDirected6173,
                       BostonScientificVercise,
                       BostonScientificVerciseDirected,
                       Medtronic3387,
                       Medtronic3389,
                       Medtronic3391,
                       MicroProbesSNEX100,
                       MicroProbesRodentElectrode,
                       PINSMedicalL301,
                       PINSMedicalL302,
                       PINSMedicalL303)
from .custom_electrodes import custom_parameters


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

ELECTRODE_MODELS = {'AbbottStJudeActiveTipCustom':
                    AbbottStJudeActiveTipModel,
                    'AbbottStJudeDirectedCustom':
                    AbbottStJudeDirectedModel,
                    'BostonScientificVerciseCustom':
                    BostonScientificVerciseModel,
                    'BostonScientificVerciseDirectedCustom':
                    BostonScientificVerciseDirectedModel,
                    'MedtronicCustom':
                    MedtronicModel,
                    'MicroProbesRodentElectrodeCustom':
                    MicroProbesRodentElectrodeModel,
                    'MicroProbesSNEX100Custom':
                    MicroProbesSNEX100Model,
                    'PINSMedicalCustom':
                    PINSMedicalModel
                    }


ELECTRODE_PARAMETERS = {'AbbottStJudeActiveTipModel':
                        AbbottStJudeParameters,
                        'AbbottStJudeDirectedModel':
                        AbbottStJudeParameters,
                        'BostonScientificVerciseModel':
                        BostonScientificVerciseParameters,
                        'BostonScientificVerciseDirectedModel':
                        BostonScientificVerciseParameters,
                        'MedtronicModel':
                        MedtronicParameters,
                        'MicroProbesRodentElectrodeModel':
                        MicroProbesRodentElectrodeParameters,
                        'MicroProbesSNEX100Model':
                        MicroProbesSNEX100Parameters,
                        'PINSMedicalModel':
                        PINSMedicalParameters
                        }

__all__ = ('ElectrodeModel',
           'AbbottStJudeActiveTip6142_6145',
           'AbbottStJudeActiveTip6146_6149',
           'AbbottStJudeActiveTipModel',
           'AbbottStJudeDirected6172',
           'AbbottStJudeDirected6173',
           'AbbottStJudeDirectedModel',
           'AbbottStJudeParameters',
           'BostonScientificVercise',
           'BostonScientificVerciseDirected',
           'BostonScientificVerciseModel',
           'BostonScientificVerciseDirectedModel',
           'BostonScientificVerciseParameters',
           'Medtronic3387',
           'Medtronic3389',
           'Medtronic3391',
           'MedtronicModel',
           'MedtronicParameters',
           'MicroProbesRodentElectrode',
           'MicroProbesRodentElectrodeModel',
           'MicroProbesRodentElectrodeParameters',
           'MicroProbesSNEX100',
           'MicroProbesSNEX100Model',
           'MicroProbesSNEX100Parameters',
           'PINSMedicalL301',
           'PINSMedicalL302',
           'PINSMedicalL303',
           'PINSMedicalModel',
           'PINSMedicalParameters',
           "default_electrode_parameters"
           )
