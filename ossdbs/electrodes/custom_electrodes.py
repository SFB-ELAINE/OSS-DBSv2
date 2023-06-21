from .abbott_stjude import AbbottStJudeParameters
from .boston_scientific_vercise import BostonScientificVerciseParameters
from .medtronic import MedtronicParameters
from .pins_medical import PINSMedicalParameters
from .micro_probes import (MicroProbesSNEX100Parameters,
                           MicroProbesRodentElectrodeParameters)


def custom_parameters(parameters):
    electrode_parameters = {
        "AbbottStJudeModel":
            AbbottStJudeParameters(parameters['TipLength[mm]'],
                                   parameters['ContactLength[mm]'],
                                   parameters['ContactSpacing[mm]'],
                                   parameters['LeadDiameter[mm]'],
                                   parameters['TotalLength[mm]']
                                   ),
        "BostonScientificVerciseModel":
            BostonScientificVerciseParameters(parameters['TipLength[mm]'],
                                              parameters['ContactLength[mm]'],
                                              parameters['ContactSpacing[mm]'],
                                              parameters['LeadDiameter[mm]'],
                                              parameters['TotalLength[mm]']
                                              ),
        "MedtronicModel":
            MedtronicParameters(parameters['TipLength[mm]'],
                                parameters['ContactLength[mm]'],
                                parameters['ContactSpacing[mm]'],
                                parameters['LeadDiameter[mm]'],
                                parameters['TotalLength[mm]']
                                ),
        "PINSMedicalModel":
            PINSMedicalParameters(parameters['TipLength[mm]'],
                                  parameters['ContactLength[mm]'],
                                  parameters['ContactSpacing[mm]'],
                                  parameters['LeadDiameter[mm]'],
                                  parameters['TotalLength[mm]']
                                  ),
        "MicroProbesRodentElectrodeModel":
            MicroProbesRodentElectrodeParameters(parameters['TubeThickness[mm]'],
                                                 parameters['ContactLength[mm]'],
                                                 parameters['LeadDiameter[mm]'],
                                                 parameters['TotalLength[mm]']
                                                 ),
        "MicroProbesSNEX100Model":
            MicroProbesSNEX100Parameters(parameters['CoreElectrodeDiameter[mm]'],
                                         parameters['CoreTubingDiameter[mm]'],
                                         parameters['CoreTubingLength[mm]'],
                                         parameters['CoreTubingDiameter[mm]'],
                                         parameters['OuterElectrodeLength[mm]'],
                                         parameters['OuterElectrodeDiameter[mm]'],
                                         parameters['OuterTubingDiameter[mm]'],
                                         parameters['TotalLength[mm]']
                                         )
                }

    electrode_parameters["AbbottStJudeDirectedModel"] = electrode_parameters["AbbottStJudeModel"]
    electrode_parameters["BostonScientificVerciseDirectedModel"] = electrode_parameters["BostonScientificVerciseModel"]
    return electrode_parameters
