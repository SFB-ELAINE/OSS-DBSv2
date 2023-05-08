from .abbott_stjude import AbbottStJudeParameters
from .boston_scientific_vercise import BostonScientificVerciseParameters
from .medtronic import MedtronicParameters
from .pins_medical import PINSMedicalParameters
from .micro_probes import (MicroProbesSNEX100Parameters,
                           MicroProbesRodentElectrodeParameters)

# TODO total_length might be too short
default_electrode_parameters = \
    {"AbbottStJudeActiveTip6146_6149":
     AbbottStJudeParameters(tip_length=1.1,
                            contact_length=1.5,
                            contact_spacing=1.5,
                            lead_diameter=1.3,
                            total_length=100.0),
     "AbbottStJudeActiveTip6142_6145":
     AbbottStJudeParameters(tip_length=2.6,
                            contact_length=1.5,
                            contact_spacing=0.5,
                            lead_diameter=1.3,
                            total_length=100.0),
     "AbbottStJudeDirected6172":
     AbbottStJudeParameters(tip_length=1.1,
                            contact_length=1.5,
                            contact_spacing=0.5,
                            lead_diameter=1.3,
                            total_length=100.0),
     "AbbottStJudeDirected6173":
     AbbottStJudeParameters(tip_length=1.1,
                            contact_length=1.5,
                            contact_spacing=1.5,
                            lead_diameter=1.3,
                            total_length=100.0),
     "BostonScientificVercise":
     BostonScientificVerciseParameters(tip_length=1.1,
                                       contact_length=1.5,
                                       contact_spacing=0.5,
                                       lead_diameter=1.3,
                                       total_length=100.0),
     "BostonScientificVerciseDirected":
     BostonScientificVerciseParameters(tip_length=1.5,
                                       contact_length=1.5,
                                       contact_spacing=0.5,
                                       lead_diameter=1.3,
                                       total_length=100.0),
     "Medtronic3387":
     MedtronicParameters(tip_length=1.5,
                         contact_length=1.5,
                         contact_spacing=0.5,
                         lead_diameter=1.27,
                         total_length=100.0),
     "Medtronic3389":
     MedtronicParameters(tip_length=1.5,
                         contact_length=1.5,
                         contact_spacing=1.5,
                         lead_diameter=1.27,
                         total_length=100.0),
     "Medtronic3391":
     MedtronicParameters(tip_length=1.5,
                         contact_length=3.0,
                         contact_spacing=3.0,
                         lead_diameter=1.27,
                         total_length=100.0),
     "MicroProbesSNEX100":
     MicroProbesSNEX100Parameters(core_electrode_length=0.25,
                                  core_electrode_diameter=0.1,
                                  core_tubing_length=0.5,
                                  core_tubing_diameter=0.140,
                                  outer_electrode_length=0.25,
                                  outer_electrode_diameter=0.330,
                                  outer_tubing_diameter=0.411,
                                  total_length=100.0),
     "MicroProbesRodentElectrode":
     MicroProbesRodentElectrodeParameters(tube_thickness=.01,
                                          contact_length=0.1125,
                                          lead_diameter=0.225,
                                          total_length=13.3),
     "PINSMedicalL301":
     PINSMedicalParameters(tip_length=1.1,
                           contact_length=1.5,
                           contact_spacing=0.5,
                           lead_diameter=1.3,
                           total_length=100.0),
     "PINSMedicalL302":
     PINSMedicalParameters(tip_length=1.1,
                           contact_length=1.5,
                           contact_spacing=1.5,
                           lead_diameter=1.3,
                           total_length=100.0),
     "PINSMedicalL303":
     PINSMedicalParameters(tip_length=1.1,
                           contact_length=3.0,
                           contact_spacing=3.0,
                           lead_diameter=1.3,
                           total_length=100.0)

     }
