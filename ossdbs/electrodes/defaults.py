# Copyright 2023, 2024 Konstantin Butenko, Shruthi Chakravarthy
# Copyright 2023, 2024 Jan Philipp Payonk, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

from .abbott_stjude import (
    AbbottStJudeActiveTipModel,
    AbbottStJudeActiveTipParameters,
    AbbottStJudeDirectedModel,
    AbbottStJudeParameters,
)
from .boston_scientific_cartesia import (
    BostonScientificCartesiaHXModel,
    BostonScientificCartesiaParameters,
    BostonScientificCartesiaXModel,
)
from .boston_scientific_vercise import (
    BostonScientificVerciseDirectedModel,
    BostonScientificVerciseDirectedParameters,
    BostonScientificVerciseModel,
    BostonScientificVerciseParameters,
)
from .dixi_microtechniques import (
    DixiSEEGModel,
    DixiSEEGParameters,
)
from .medtronic import (
    MedtronicModel,
    MedtronicParameters,
    MedtronicSenSightModel,
)
from .micro_probes import (
    MicroProbesRodentElectrodeModel,
    MicroProbesRodentElectrodeParameters,
    MicroProbesSNEX100Model,
    MicroProbesSNEX100Parameters,
)
from .microelectrode import (
    MicroElectrodeModel,
    MicroElectrodeParameters,
)
from .neuro_pace import (
    NeuroPaceModel,
    NeuroPaceParameters,
)
from .neuronexus import NeuroNexusElectrodeModel, NeuroNexusParameters
from .pins_medical import (
    PINSMedicalModel,
    PINSMedicalParameters,
)
from .sceneray import SceneRay1242Model, SceneRay1242Parameters

default_electrode_parameters = {
    "AbbottStJudeActiveTip6146_6149": AbbottStJudeActiveTipParameters(
        tip_length=3.0,
        contact_length=1.5,
        contact_spacing=0.5,
        lead_diameter=1.4,
        total_length=400.0,
    ),
    "AbbottStJudeActiveTip6142_6145": AbbottStJudeActiveTipParameters(
        tip_length=3.0,
        contact_length=1.5,
        contact_spacing=1.5,
        lead_diameter=1.4,
        total_length=400.0,
    ),
    "AbbottStJudeDirected6172": AbbottStJudeParameters(
        tip_length=1.0,
        contact_length=1.5,
        contact_spacing=0.5,
        lead_diameter=1.29,
        total_length=400.0,
    ),
    "AbbottStJudeDirected6173": AbbottStJudeParameters(
        tip_length=1.0,
        contact_length=1.5,
        contact_spacing=1.5,
        lead_diameter=1.29,
        total_length=400.0,
    ),
    "BostonScientificVercise": BostonScientificVerciseParameters(
        tip_length=1.1,
        contact_length=1.5,
        contact_spacing=0.5,
        lead_diameter=1.3,
        total_length=450.0,
    ),
    "BostonScientificVerciseDirected": BostonScientificVerciseDirectedParameters(
        tip_length=1.5,
        contact_length=1.5,
        contact_spacing=0.5,
        lead_diameter=1.3,
        total_length=450.0,
    ),
    "BostonScientificCartesiaX": BostonScientificCartesiaParameters(
        tip_length=1.1,
        contact_length=1.5,
        contact_spacing=0.5,
        lead_diameter=1.3,
        total_length=450.0,
        contacts_skipped=5.0,
    ),
    "BostonScientificCartesiaHX": BostonScientificCartesiaParameters(
        tip_length=1.1,
        contact_length=1.5,
        contact_spacing=0.5,
        lead_diameter=1.3,
        total_length=450.0,
        contacts_skipped=7.0,
    ),
    "Medtronic3387": MedtronicParameters(
        tip_length=1.5,
        contact_length=1.5,
        contact_spacing=1.5,
        lead_diameter=1.27,
        total_length=400.0,
    ),
    "Medtronic3389": MedtronicParameters(
        tip_length=1.5,
        contact_length=1.5,
        contact_spacing=0.5,
        lead_diameter=1.27,
        total_length=400.0,
    ),
    "Medtronic3391": MedtronicParameters(
        tip_length=1.5,
        contact_length=3.0,
        contact_spacing=4.0,
        lead_diameter=1.27,
        total_length=400.0,
    ),
    "MedtronicSenSightB33015": MedtronicParameters(
        tip_length=0.9,
        contact_length=1.5,
        contact_spacing=1.5,
        lead_diameter=1.36,
        total_length=330,
    ),
    "MedtronicSenSightB33005": MedtronicParameters(
        tip_length=0.9,
        contact_length=1.5,
        contact_spacing=0.5,
        lead_diameter=1.36,
        total_length=330,
    ),
    "MicroElectrode": MicroElectrodeParameters(
        tip_length=1.0,
        tip_diameter=0.7,
        contact_length=1.0,
        lead_diameter=1.0,
        total_length=200.0,
    ),
    "MicroProbesSNEX100": MicroProbesSNEX100Parameters(
        core_electrode_length=0.25,
        core_electrode_diameter=0.1,
        core_tubing_length=0.5,
        core_tubing_diameter=0.140,
        outer_electrode_length=0.25,
        outer_electrode_diameter=0.330,
        outer_tubing_diameter=0.411,
        total_length=100.0,
    ),
    "MicroProbesRodentElectrode": MicroProbesRodentElectrodeParameters(
        exposed_wire=0,
        contact_radius=0.1125,
        lead_radius=0.1175,
        total_length=100.0,
        wire_radius=0.1125,
    ),
    "NeuroPaceDL344_3_5": NeuroPaceParameters(
        tip_length=1.1,
        contact_length=2.0,
        contact_spacing=1.5,
        lead_diameter=1.27,
        total_length=450.0,
    ),
    "NeuroPaceDL344_10": NeuroPaceParameters(
        tip_length=1.1,
        contact_length=2.0,
        contact_spacing=8.0,
        lead_diameter=1.27,
        total_length=450.0,
    ),
    "PINSMedicalL301": PINSMedicalParameters(
        tip_length=1.5,
        contact_length=1.5,
        contact_spacing=0.5,
        lead_diameter=1.3,
        total_length=400.0,
    ),
    "PINSMedicalL302": PINSMedicalParameters(
        tip_length=1.5,
        contact_length=1.5,
        contact_spacing=1.5,
        lead_diameter=1.3,
        total_length=400.0,
    ),
    "PINSMedicalL303": PINSMedicalParameters(
        tip_length=1.5,
        contact_length=3.0,
        contact_spacing=3.0,
        lead_diameter=1.3,
        total_length=400.0,
    ),
    "DixiSEEG5": DixiSEEGParameters(
        tip_length=2.0,
        contact_length=2.0,
        contact_spacing=1.5,
        lead_diameter=0.8,
        total_length=400.0,
        n_contacts=5,
    ),
    "DixiSEEG8": DixiSEEGParameters(
        tip_length=2.0,
        contact_length=2.0,
        contact_spacing=1.5,
        lead_diameter=0.8,
        total_length=400.0,
        n_contacts=8,
    ),
    "DixiSEEG10": DixiSEEGParameters(
        tip_length=2.0,
        contact_length=2.0,
        contact_spacing=1.5,
        lead_diameter=0.8,
        total_length=400.0,
        n_contacts=10,
    ),
    "DixiSEEG12": DixiSEEGParameters(
        tip_length=2.0,
        contact_length=2.0,
        contact_spacing=1.5,
        lead_diameter=0.8,
        total_length=400.0,
        n_contacts=12,
    ),
    "DixiSEEG15": DixiSEEGParameters(
        tip_length=2.0,
        contact_length=2.0,
        contact_spacing=1.5,
        lead_diameter=0.8,
        total_length=400.0,
        n_contacts=15,
    ),
    "DixiSEEG18": DixiSEEGParameters(
        tip_length=2.0,
        contact_length=2.0,
        contact_spacing=1.5,
        lead_diameter=0.8,
        total_length=400.0,
        n_contacts=18,
    ),
    "PMTsEEG2102_08": DixiSEEGParameters(
        tip_length=2.0,
        contact_length=2.0,
        contact_spacing=1.5,
        lead_diameter=0.8,
        total_length=400.0,
        n_contacts=8,
    ),
    "PMTsEEG2102_10": DixiSEEGParameters(
        tip_length=2.0,
        contact_length=2.0,
        contact_spacing=1.5,
        lead_diameter=0.8,
        total_length=400.0,
        n_contacts=10,
    ),
    "PMTsEEG2102_12": DixiSEEGParameters(
        tip_length=2.0,
        contact_length=2.0,
        contact_spacing=1.5,
        lead_diameter=0.8,
        total_length=400.0,
        n_contacts=12,
    ),
    "PMTsEEG2102_14": DixiSEEGParameters(
        tip_length=2.0,
        contact_length=2.0,
        contact_spacing=1.5,
        lead_diameter=0.8,
        total_length=400.0,
        n_contacts=14,
    ),
    "PMTsEEG2102_16": DixiSEEGParameters(
        tip_length=2.0,
        contact_length=2.0,
        contact_spacing=1.5,
        lead_diameter=0.8,
        total_length=400.0,
        n_contacts=16,
    ),
    "NeuroNexusA1x16_5mm_50_177": NeuroNexusParameters(
        shank_thickness=20.0e-3,
        max_width=125.0e-3,
        min_width=33.0e-3,
        total_length=400.0,
        angle_tip=30.0,
        tip_length=50.0e-3,
        contact_spacing=50.0e-3,
        contact_diameter=30.0e-3,
    ),
    "SceneRay1242": SceneRay1242Parameters(
        tip_length=1.5,
        contact_length=3.0,
        first_contact_spacing=2.0,
        contact_spacing=4.0,
        lead_diameter=1.27,
        total_length=400.0,
    ),
    "SceneRay1202": BostonScientificVerciseParameters(
        tip_length=1.5,
        contact_length=1.5,
        contact_spacing=0.5,
        lead_diameter=1.27,
        total_length=450.0,
    ),
    "SceneRay1211": MedtronicParameters(
        tip_length=1.5,
        contact_length=1.5,
        contact_spacing=1.0,
        lead_diameter=1.27,
        total_length=400.0,
    ),
    "SceneRay1212": BostonScientificVerciseParameters(
        tip_length=1.5,
        contact_length=1.5,
        contact_spacing=1.5,
        lead_diameter=1.27,
        total_length=450.0,
    ),
}


def AbbottStJudeActiveTip6142_6145(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """Abbott St. Jude Active Tip 6142 6145 electrode."""
    parameters = default_electrode_parameters["AbbottStJudeActiveTip6142_6145"]
    return AbbottStJudeActiveTipModel(parameters, rotation, direction, position)


def AbbottStJudeActiveTip6146_6149(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """Abbott St. Jude Active Tip 6146 6149 electrode."""
    parameters = default_electrode_parameters["AbbottStJudeActiveTip6146_6149"]
    return AbbottStJudeActiveTipModel(parameters, rotation, direction, position)


def AbbottStJudeDirected6172(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """Abbott St. Jude Directed Tip 6172 electrode."""
    parameters = default_electrode_parameters["AbbottStJudeDirected6172"]
    return AbbottStJudeDirectedModel(parameters, rotation, direction, position)


def AbbottStJudeDirected6173(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """Abbott St. Jude Directed Tip 6173 electrode."""
    parameters = default_electrode_parameters["AbbottStJudeDirected6173"]
    return AbbottStJudeDirectedModel(parameters, rotation, direction, position)


def BostonScientificVercise(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """Boston Scientific Vercise electrode."""
    parameters = default_electrode_parameters["BostonScientificVercise"]
    return BostonScientificVerciseModel(parameters, rotation, direction, position)


def BostonScientificVerciseDirected(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """Boston Scientific Vercise Directed electrode."""
    parameters = default_electrode_parameters["BostonScientificVerciseDirected"]
    return BostonScientificVerciseDirectedModel(
        parameters, rotation, direction, position
    )


def BostonScientificCartesiaX(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """Boston Scientific Cartesia X electrode."""
    parameters = default_electrode_parameters["BostonScientificCartesiaX"]
    return BostonScientificCartesiaXModel(parameters, rotation, direction, position)


def BostonScientificCartesiaHX(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """Boston Scientific Cartesia HX electrode."""
    parameters = default_electrode_parameters["BostonScientificCartesiaHX"]
    return BostonScientificCartesiaHXModel(parameters, rotation, direction, position)


def Medtronic3387(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """Medtronic 3387 electrode."""
    parameters = default_electrode_parameters["Medtronic3387"]
    return MedtronicModel(parameters, rotation, direction, position)


def Medtronic3389(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """Medtronic 3389 electrode."""
    parameters = default_electrode_parameters["Medtronic3389"]
    return MedtronicModel(parameters, rotation, direction, position)


def Medtronic3391(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """Medtronic 3391 electrode."""
    parameters = default_electrode_parameters["Medtronic3391"]
    return MedtronicModel(parameters, rotation, direction, position)


def MedtronicSenSightB33015(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """Medtronic SenSight B33015 electrode."""
    parameters = default_electrode_parameters["MedtronicSenSightB33015"]
    return MedtronicSenSightModel(parameters, rotation, direction, position)


def MedtronicSenSightB33005(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """Medtronic SenSight B33005 electrode."""
    parameters = default_electrode_parameters["MedtronicSenSightB33005"]
    return MedtronicSenSightModel(parameters, rotation, direction, position)


def MicroElectrode(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """MicroElectrode electrode."""
    parameters = default_electrode_parameters["MicroElectrode"]
    return MicroElectrodeModel(parameters, rotation, direction, position)


def MicroProbesSNEX100(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """MicroProbes SNEX 100 electrode."""
    parameters = default_electrode_parameters["MicroProbesSNEX100"]
    return MicroProbesSNEX100Model(parameters, rotation, direction, position)


def MicroProbesRodentElectrode(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """MicroProbes rodent (unipolar) electrode."""
    parameters = default_electrode_parameters["MicroProbesRodentElectrode"]
    return MicroProbesRodentElectrodeModel(parameters, rotation, direction, position)


def NeuroPaceDL344_3_5(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """NeuroPace DL344_3_5 electrode."""
    parameters = default_electrode_parameters["NeuroPaceDL344_3_5"]
    return NeuroPaceModel(parameters, rotation, direction, position)


def NeuroPaceDL344_10(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """NeuroPace DL344_10 electrode."""
    parameters = default_electrode_parameters["NeuroPaceDL344_10"]
    return NeuroPaceModel(parameters, rotation, direction, position)


def PINSMedicalL301(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """PINSMedical L301 electrode."""
    parameters = default_electrode_parameters["PINSMedicalL301"]
    return PINSMedicalModel(parameters, rotation, direction, position)


def PINSMedicalL302(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """PINSMedical L302 electrode."""
    parameters = default_electrode_parameters["PINSMedicalL302"]
    return PINSMedicalModel(parameters, rotation, direction, position)


def PINSMedicalL303(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """PINSMedical L303 electrode."""
    parameters = default_electrode_parameters["PINSMedicalL303"]
    return PINSMedicalModel(parameters, rotation, direction, position)


def DixiSEEG5(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """DixiSEEG5 electrode."""
    parameters = default_electrode_parameters["DixiSEEG5"]
    return DixiSEEGModel(parameters, rotation, direction, position)


def DixiSEEG8(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """DixiSEEG8 electrode."""
    parameters = default_electrode_parameters["DixiSEEG8"]
    return DixiSEEGModel(parameters, rotation, direction, position)


def DixiSEEG10(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """DixiSEEG10 electrode."""
    parameters = default_electrode_parameters["DixiSEEG10"]
    return DixiSEEGModel(parameters, rotation, direction, position)


def DixiSEEG12(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """DixiSEEG12 electrode."""
    parameters = default_electrode_parameters["DixiSEEG12"]
    return DixiSEEGModel(parameters, rotation, direction, position)


def DixiSEEG15(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """DixiSEEG15 electrode."""
    parameters = default_electrode_parameters["DixiSEEG15"]
    return DixiSEEGModel(parameters, rotation, direction, position)


def DixiSEEG18(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """DixiSEEG18 electrode."""
    parameters = default_electrode_parameters["DixiSEEG18"]
    return DixiSEEGModel(parameters, rotation, direction, position)


def PMTsEEG2102_08(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """PMTsEEG2102_08 electrode."""
    parameters = default_electrode_parameters["PMTsEEG2102_08"]
    return DixiSEEGModel(parameters, rotation, direction, position)


def PMTsEEG2102_10(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """PMTsEEG2102_10 electrode."""
    parameters = default_electrode_parameters["PMTsEEG2102_10"]
    return DixiSEEGModel(parameters, rotation, direction, position)


def PMTsEEG2102_12(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """PMTsEEG2102_12 electrode."""
    parameters = default_electrode_parameters["PMTsEEG2102_12"]
    return DixiSEEGModel(parameters, rotation, direction, position)


def PMTsEEG2102_14(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """PMTsEEG2102_14 electrode."""
    parameters = default_electrode_parameters["PMTsEEG2102_14"]
    return DixiSEEGModel(parameters, rotation, direction, position)


def PMTsEEG2102_16(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """PMTsEEG2102_16 electrode."""
    parameters = default_electrode_parameters["PMTsEEG2102_16"]
    return DixiSEEGModel(parameters, rotation, direction, position)


def NeuroNexusA1x16_5mm_50_177(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """NeuroNexusA1x16_5mm_50_177 electrode."""
    parameters = default_electrode_parameters["NeuroNexusA1x16_5mm_50_177"]
    return NeuroNexusElectrodeModel(parameters, rotation, direction, position)


def SceneRay1242(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """SceneRay 1242 electrode."""
    parameters = default_electrode_parameters["SceneRay1242"]
    return SceneRay1242Model(parameters, rotation, direction, position)


def SceneRay1202(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """SceneRay 1202 electrode."""
    parameters = default_electrode_parameters["SceneRay1202"]
    return BostonScientificVerciseModel(parameters, rotation, direction, position)


def SceneRay1211(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """SceneRay 1211 electrode."""
    parameters = default_electrode_parameters["SceneRay1211"]
    return MedtronicModel(parameters, rotation, direction, position)


def SceneRay1212(
    rotation: float = 0, direction: tuple = (0, 0, 1), position: tuple = (0, 0, 0)
):
    """SceneRay 1212 electrode."""
    parameters = default_electrode_parameters["SceneRay1212"]
    return BostonScientificVerciseModel(parameters, rotation, direction, position)
