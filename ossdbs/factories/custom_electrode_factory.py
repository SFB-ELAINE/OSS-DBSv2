
from ossdbs.electrodes.electrode_models import ElectrodeModel

from ossdbs.electrodes.electrode_models import AbbottStJudeActiveTipModel
from ossdbs.electrodes.electrode_models import AbbottStJudeDirectedModel
from ossdbs.electrodes.electrode_models import BostonScientificVerciseModel
from ossdbs.electrodes.electrode_models \
                                    import BostonScientificVerciseDirectedModel
from ossdbs.electrodes.electrode_models import MedtronicModel
from ossdbs.electrodes.electrode_models import MicroProbesRodentElectrodeModel
from ossdbs.electrodes.electrode_models import MicroProbesSNEX100Model
from ossdbs.electrodes.electrode_models import PINSMedicalModel

from ossdbs.electrodes.electrode_models import AbbottStJudeParameters
from ossdbs.electrodes.electrode_models import BostonScientificVerciseParameters
from ossdbs.electrodes.electrode_models import MedtronicParameters
from ossdbs.electrodes.electrode_models \
                                    import MicroProbesRodentElectrodeParameters
from ossdbs.electrodes.electrode_models import MicroProbesSNEX100Parameters
from ossdbs.electrodes.electrode_models import PINSMedicalParameters


class CustomElectrodeFactory:
    """Creates a list of Electrode objects."""

    ELECTRODES = {'AbbottStJudeActiveTipCustom':
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

    @classmethod
    def create(cls,
               name: str,
               model_parameters: dict,
               direction: tuple,
               position: tuple,
               rotation: float
               ) -> ElectrodeModel:
        """create a list of Electrode objects.

        Parameters
        ----------
        parameters : ElectrodeParameters

        Returns
        -------
        Electrode
            Electrode objects.
        """

        electrode_type = cls.ELECTRODES[name]

        if electrode_type in (AbbottStJudeActiveTipModel,
                              AbbottStJudeDirectedModel):
            parameters = AbbottStJudeParameters(
                                        model_parameters['TipLength[mm]'],
                                        model_parameters['ContactLength[mm]'],
                                        model_parameters['ContactSpacing[mm]'],
                                        model_parameters['LeadDiameter[mm]'],
                                        model_parameters['TotalLength[mm]']
                                        )

        elif electrode_type in (BostonScientificVerciseModel,
                                BostonScientificVerciseDirectedModel):
            parameters = BostonScientificVerciseParameters(
                                        model_parameters['TipLength[mm]'],
                                        model_parameters['ContactLength[mm]'],
                                        model_parameters['ContactSpacing[mm]'],
                                        model_parameters['LeadDiameter[mm]'],
                                        model_parameters['TotalLength[mm]']
                                        )

        elif electrode_type == MedtronicModel:
            parameters = MedtronicParameters(
                                        model_parameters['TipLength[mm]'],
                                        model_parameters['ContactLength[mm]'],
                                        model_parameters['ContactSpacing[mm]'],
                                        model_parameters['LeadDiameter[mm]'],
                                        model_parameters['TotalLength[mm]']
                                        )

        elif electrode_type == PINSMedicalModel:
            parameters = PINSMedicalParameters(
                                        model_parameters['TipLength[mm]'],
                                        model_parameters['ContactLength[mm]'],
                                        model_parameters['ContactSpacing[mm]'],
                                        model_parameters['LeadDiameter[mm]'],
                                        model_parameters['TotalLength[mm]']
                                        )

        elif electrode_type == MicroProbesRodentElectrodeModel:
            parameters = MicroProbesRodentElectrodeParameters(
                                        model_parameters['TubeThickness[mm]'],
                                        model_parameters['ContactLength[mm]'],
                                        model_parameters['LeadDiameter[mm]'],
                                        model_parameters['TotalLength[mm]']
                                        )

        elif electrode_type == MicroProbesSNEX100Model:
            parameters = MicroProbesSNEX100Parameters(
                                model_parameters['CoreElectrodeDiameter[mm]'],
                                model_parameters['CoreTubingDiameter[mm]'],
                                model_parameters['CoreTubingLength[mm]'],
                                model_parameters['CoreTubingDiameter[mm]'],
                                model_parameters['OuterElectrodeLength[mm]'],
                                model_parameters['OuterElectrodeDiameter[mm]'],
                                model_parameters['OuterTubingDiameter[mm]'],
                                model_parameters['TotalLength[mm]']
                                )
        else:
            raise NotImplementedError("Electrode model not implemented.")

        return electrode_type(parameters=parameters,
                              direction=direction,
                              position=position,
                              rotation=rotation)
