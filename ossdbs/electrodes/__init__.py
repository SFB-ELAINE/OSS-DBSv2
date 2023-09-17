from .abbott_stjude import (
    AbbottStJudeActiveTipModel,
    AbbottStJudeDirectedModel,
    AbbottStJudeParameters,
)
from .boston_scientific_vercise import (
    BostonScientificVerciseDirectedModel,
    BostonScientificVerciseModel,
    BostonScientificVerciseParameters,
)
from .defaults import (
    AbbottStJudeActiveTip6142_6145,
    AbbottStJudeActiveTip6146_6149,
    AbbottStJudeDirected6172,
    AbbottStJudeDirected6173,
    BostonScientificVercise,
    BostonScientificVerciseDirected,
    Medtronic3387,
    Medtronic3389,
    Medtronic3391,
    MedtronicSenSightB33005,
    MedtronicSenSightB33015,
    MicroProbesRodentElectrode,
    MicroProbesSNEX100,
    PINSMedicalL301,
    PINSMedicalL302,
    PINSMedicalL303,
    default_electrode_parameters,
)
from .electrode_model_template import ElectrodeModel
from .medtronic import MedtronicModel, MedtronicParameters, MedtronicSenSightModel
from .micro_probes import (
    MicroProbesRodentElectrodeModel,
    MicroProbesRodentElectrodeParameters,
    MicroProbesSNEX100Model,
    MicroProbesSNEX100Parameters,
)
from .pins_medical import PINSMedicalModel, PINSMedicalParameters

ELECTRODES = {
    "AbbottStJudeActiveTip6142_6145": AbbottStJudeActiveTip6142_6145,
    "AbbottStJudeActiveTip6146_6149": AbbottStJudeActiveTip6146_6149,
    "AbbottStJudeDirected6172": AbbottStJudeDirected6172,
    "AbbottStJudeDirected6173": AbbottStJudeDirected6173,
    "BostonScientificVercise": BostonScientificVercise,
    "BostonScientificVerciseDirected": BostonScientificVerciseDirected,
    "Medtronic3387": Medtronic3387,
    "Medtronic3389": Medtronic3389,
    "Medtronic3391": Medtronic3391,
    "MedtronicSenSightB33015": MedtronicSenSightB33015,
    "MedtronicSenSightB33005": MedtronicSenSightB33005,
    "MicroProbesSNEX100": MicroProbesSNEX100,
    "PINSMedicalL301": PINSMedicalL301,
    "PINSMedicalL302": PINSMedicalL302,
    "PINSMedicalL303": PINSMedicalL303,
    "MicroProbesRodentElectrode": MicroProbesRodentElectrode,
}

ELECTRODE_MODELS = {
    "AbbottStJudeActiveTipCustom": AbbottStJudeActiveTipModel,
    "AbbottStJudeDirectedCustom": AbbottStJudeDirectedModel,
    "BostonScientificVerciseCustom": BostonScientificVerciseModel,
    "BostonScientificVerciseDirectedCustom": BostonScientificVerciseDirectedModel,
    "MedtronicCustom": MedtronicModel,
    "MedtronicSensightCustom": MedtronicSenSightModel,
    "MicroProbesRodentElectrodeCustom": MicroProbesRodentElectrodeModel,
    "MicroProbesSNEX100Custom": MicroProbesSNEX100Model,
    "PINSMedicalCustom": PINSMedicalModel,
}


ELECTRODE_PARAMETERS = {
    "AbbottStJudeActiveTipModel": AbbottStJudeParameters,
    "AbbottStJudeDirectedModel": AbbottStJudeParameters,
    "BostonScientificVerciseModel": BostonScientificVerciseParameters,
    "BostonScientificVerciseDirectedModel": BostonScientificVerciseParameters,
    "MedtronicModel": MedtronicParameters,
    "MedtronicSenSightModel": MedtronicParameters,
    "MicroProbesRodentElectrodeModel": MicroProbesRodentElectrodeParameters,
    "MicroProbesSNEX100Model": MicroProbesSNEX100Parameters,
    "PINSMedicalModel": PINSMedicalParameters,
}

__all__ = (
    "ElectrodeModel",
    "AbbottStJudeActiveTip6142_6145",
    "AbbottStJudeActiveTip6146_6149",
    "AbbottStJudeActiveTipModel",
    "AbbottStJudeDirected6172",
    "AbbottStJudeDirected6173",
    "AbbottStJudeDirectedModel",
    "AbbottStJudeParameters",
    "BostonScientificVercise",
    "BostonScientificVerciseDirected",
    "BostonScientificVerciseModel",
    "BostonScientificVerciseDirectedModel",
    "BostonScientificVerciseParameters",
    "Medtronic3387",
    "Medtronic3389",
    "Medtronic3391",
    "MedtronicModel",
    "MedtronicSenSightB33015",
    "MedtronicSenSightB33005",
    "MedtronicParameters",
    "MicroProbesRodentElectrode",
    "MicroProbesRodentElectrodeModel",
    "MicroProbesRodentElectrodeParameters",
    "MicroProbesSNEX100",
    "MicroProbesSNEX100Model",
    "MicroProbesSNEX100Parameters",
    "PINSMedicalL301",
    "PINSMedicalL302",
    "PINSMedicalL303",
    "PINSMedicalModel",
    "PINSMedicalParameters",
    "default_electrode_parameters",
)
