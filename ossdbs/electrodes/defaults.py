from .abbott_stjude import (AbbottStJudeParameters,
                            AbbottStJudeActiveTipModel,
                            AbbottStJudeDirectedModel)
from .boston_scientific_vercise import (BostonScientificVerciseParameters,
                                        BostonScientificVerciseModel,
                                        BostonScientificVerciseDirectedModel)
from .medtronic import (MedtronicParameters,
                        MedtronicModel, MedtronicSenSightModel)
from .micro_probes import (MicroProbesSNEX100Parameters,
                           MicroProbesSNEX100Model,
                           MicroProbesRodentElectrodeParameters,
                           MicroProbesRodentElectrodeModel)
from .pins_medical import (PINSMedicalParameters,
                           PINSMedicalModel)


# TODO total_length might be too short
default_electrode_parameters = \
    {"AbbottStJudeActiveTip6146_6149":
     AbbottStJudeParameters(tip_length=3.0,
                            contact_length=1.5,
                            contact_spacing=0.5,
                            lead_diameter=1.4,
                            total_length=400.0,
                            offset=1.5),
     "AbbottStJudeActiveTip6142_6145":
     AbbottStJudeParameters(tip_length=3.0,
                            contact_length=1.5,
                            contact_spacing=1.5,
                            lead_diameter=1.4,
                            total_length=400.0,
                            offset=1.5),
     "AbbottStJudeDirected6172":
     AbbottStJudeParameters(tip_length=1.5,
                            contact_length=1.5,
                            contact_spacing=0.5,
                            lead_diameter=1.29,
                            total_length=400.0,
                            offset=2.25),
     "AbbottStJudeDirected6173":
     AbbottStJudeParameters(tip_length=1.5,
                            contact_length=1.5,
                            contact_spacing=1.5,
                            lead_diameter=1.29,
                            total_length=400.0,
                            offset=2.25),
     # The first version of ossdbs uses an offset of 1.85
     # for BostonScientificVericise, but 2.04 fits better here -
     # is it because there are some differences between
     # Boston Scientific Vercise models?
     "BostonScientificVercise":
     BostonScientificVerciseParameters(tip_length=1.29,
                                       contact_length=1.5,
                                       contact_spacing=0.5,
                                       lead_diameter=1.3,
                                       total_length=450.0,
                                       offset=2.04),
     "BostonScientificVerciseDirected":
     BostonScientificVerciseParameters(tip_length=1.5,
                                       contact_length=1.5,
                                       contact_spacing=0.5,
                                       lead_diameter=1.3,
                                       total_length=450.0,
                                       offset=0.75),
     "Medtronic3387":
     MedtronicParameters(tip_length=1.5,
                         contact_length=1.5,
                         contact_spacing=1.5,
                         lead_diameter=1.27,
                         total_length=400.0,
                         offset=2.25),
     "Medtronic3389":
     MedtronicParameters(tip_length=1.5,
                         contact_length=1.5,
                         contact_spacing=0.5,
                         lead_diameter=1.27,
                         total_length=400.0,
                         offset=2.25),
     "Medtronic3391":
     MedtronicParameters(tip_length=1.5,
                         contact_length=3.0,
                         contact_spacing=4.0,
                         lead_diameter=1.27,
                         total_length=400.0,
                         offset=3.0),
     "MedtronicSenSightB33015":
     MedtronicParameters(tip_length=1.0,
                         contact_length=1.5,
                         contact_spacing=1.5,
                         lead_diameter=1.36,
                         total_length=330,
                         offset=1.75),
     "MedtronicSenSightB33005":
     MedtronicParameters(tip_length=1.0,
                         contact_length=1.5,
                         contact_spacing=0.5,
                         lead_diameter=1.36,
                         total_length=330,
                         offset=1.75),
     "MicroProbesSNEX100":
     MicroProbesSNEX100Parameters(core_electrode_length=0.25,
                                  core_electrode_diameter=0.1,
                                  core_tubing_length=0.5,
                                  core_tubing_diameter=0.140,
                                  outer_electrode_length=0.25,
                                  outer_electrode_diameter=0.330,
                                  outer_tubing_diameter=0.411,
                                  total_length=100.0,
                                  offset=0),
     "MicroProbesRodentElectrode":
     MicroProbesRodentElectrodeParameters(exposed_wire=0,
                                          contact_radius=0.1125,
                                          lead_radius=0.1175,
                                          total_length=13.3,
                                          wire_radius=0.0,
                                          offset=0),
     "PINSMedicalL301":
     PINSMedicalParameters(tip_length=1.5,
                           contact_length=1.5,
                           contact_spacing=0.5,
                           lead_diameter=1.3,
                           total_length=400.0,
                           offset=2.25),
     "PINSMedicalL302":
     PINSMedicalParameters(tip_length=1.5,
                           contact_length=1.5,
                           contact_spacing=1.5,
                           lead_diameter=1.3,
                           total_length=400.0,
                           offset=2.25),
     "PINSMedicalL303":
     PINSMedicalParameters(tip_length=1.5,
                           contact_length=3.0,
                           contact_spacing=3.0,
                           lead_diameter=1.3,
                           total_length=400.0,
                           offset=2.25)

     }


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
    parameters = default_electrode_parameters["AbbottStJudeDirected6172"]
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


def MedtronicSenSightB33015(rotation: float = 0,
                            direction: tuple = (0, 0, 1),
                            position: tuple = (0, 0, 0)):
    parameters = default_electrode_parameters["MedtronicSenSightB33015"]
    return MedtronicSenSightModel(parameters, rotation, direction, position)


def MedtronicSenSightB33005(rotation: float = 0,
                            direction: tuple = (0, 0, 1),
                            position: tuple = (0, 0, 0)):
    parameters = default_electrode_parameters["MedtronicSenSightB33005"]
    return MedtronicSenSightModel(parameters, rotation, direction, position)


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