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
from .defaults import default_electrode_parameters


def AbbottStJudeActiveTip6142_6145(rotation: float = 0,
                                   direction: tuple = (0, 0, 1),
                                   position: tuple = (0, 0, 0)):
    parameters = default_electrode_parameters["AbbottStJudeActiveTip6142_6145"]
    return AbbottStJudeActiveTipModel(parameters, rotation, direction, position)


def AbbottStJudeActiveTip6146_6149(rotation: float = 0,
                                   direction: tuple = (0, 0, 1),
                                   position: tuple = (0, 0, 0)):
    parameters = default_electrode_parameters["AbbottStJudeActiveTip6146_6149"]
    return AbbottStJudeActiveTipModel(parameters, rotation, direction, position)


def AbbottStJudeDirected6172(rotation: float = 0,
                             direction: tuple = (0, 0, 1),
                             position: tuple = (0, 0, 0)):
    parameters = default_electrode_parameters["AbbottStJudeDirected6173"]
    return AbbottStJudeDirectedModel(parameters, rotation, direction, position)


def AbbottStJudeDirected6173(rotation: float = 0,
                             direction: tuple = (0, 0, 1),
                             position: tuple = (0, 0, 0)):
    parameters = default_electrode_parameters["AbbottStJudeDirected6173"]
    return AbbottStJudeDirectedModel(parameters, rotation, direction, position)


def BostonScientificVercise(rotation: float = 0,
                            direction: tuple = (0, 0, 1),
                            position: tuple = (0, 0, 0)):
    parameters = default_electrode_parameters["BostonScientificVercise"]
    return BostonScientificVerciseModel(parameters, rotation, direction, position)


def BostonScientificVerciseDirected(rotation: float = 0,
                                    direction: tuple = (0, 0, 1),
                                    position: tuple = (0, 0, 0)):
    parameters = default_electrode_parameters["BostonScientificVerciseDirected"]
    return BostonScientificVerciseDirectedModel(parameters, rotation, direction, position)


def Medtronic3387(rotation: float = 0,
                  direction: tuple = (0, 0, 1),
                  position: tuple = (0, 0, 0)):
    parameters = default_electrode_parameters["Medtronic3387"]
    return MedtronicModel(parameters, rotation, direction, position)


def Medtronic3389(rotation: float = 0,
                  direction: tuple = (0, 0, 1),
                  position: tuple = (0, 0, 0)):
    parameters = default_electrode_parameters["Medtronic3389"]
    return MedtronicModel(parameters, rotation, direction, position)


def Medtronic3391(rotation: float = 0,
                  direction: tuple = (0, 0, 1),
                  position: tuple = (0, 0, 0)):
    parameters = default_electrode_parameters["Medtronic3391"]
    return MedtronicModel(parameters, rotation, direction, position)


def MicroProbesSNEX100(rotation: float = 0,
                       direction: tuple = (0, 0, 1),
                       position: tuple = (0, 0, 0)):
    parameters = default_electrode_parameters["MicroProbesSNEX100"]
    return MicroProbesSNEX100Model(parameters, rotation, direction, position)


def MicroProbesRodentElectrode(rotation: float = 0,
                               direction: tuple = (0, 0, 1),
                               position: tuple = (0, 0, 0)):
    parameters = default_electrode_parameters["MicroProbesRodentElectrode"]
    return MicroProbesRodentElectrodeModel(parameters, rotation, direction, position)


def PINSMedicalL301(rotation: float = 0,
                    direction: tuple = (0, 0, 1),
                    position: tuple = (0, 0, 0)):
    parameters = default_electrode_parameters["PINSMedicalL301"]
    return PINSMedicalModel(parameters, rotation, direction, position)


def PINSMedicalL302(rotation: float = 0,
                    direction: tuple = (0, 0, 1),
                    position: tuple = (0, 0, 0)):
    parameters = default_electrode_parameters["PINSMedicalL302"]
    return PINSMedicalModel(parameters, rotation, direction, position)


def PINSMedicalL303(rotation: float = 0,
                    direction: tuple = (0, 0, 1),
                    position: tuple = (0, 0, 0)):
    parameters = default_electrode_parameters["PINSMedicalL303"]
    return PINSMedicalModel(parameters, rotation, direction, position)


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
