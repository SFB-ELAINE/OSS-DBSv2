# MicroProbes SNEX 100 Concentric Bipolar Electrodes
from .micro_probes_SNEX100_model import MicroProbesSNEX100Parameters
from .micro_probes_SNEX100_model import MicroProbesSNEX100Model
from .micro_probes_rodent_electrode_model \
                                    import MicroProbesRodentElectrodeModel
from .micro_probes_rodent_electrode_model \
                                    import MicroProbesRodentElectrodeParameters


class MicroProbesSNEX100(MicroProbesSNEX100Model):

    def __init__(self,
                 rotation: float = 0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        parameters = MicroProbesSNEX100Parameters(
                                                core_electrode_length=0.25,
                                                core_electrode_diameter=0.1,
                                                core_tubing_length=0.5,
                                                core_tubing_diameter=0.140,
                                                outer_electrode_length=0.25,
                                                outer_electrode_diameter=0.330,
                                                outer_tubing_diameter=0.411,
                                                total_length=100.0)
        super().__init__(parameters, rotation, direction, position)


class MicroProbesRodentElectrode(MicroProbesRodentElectrodeModel):

    def __init__(self,
                 rotation: float = 0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)
                 ) -> None:

        parameters = MicroProbesRodentElectrodeParameters(
                                                        tube_thickness=.01,
                                                        contact_length=0.1125,
                                                        lead_diameter=0.225,
                                                        total_length=13.3)
        super().__init__(parameters, rotation, direction, position)
