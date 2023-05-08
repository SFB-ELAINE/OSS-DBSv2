from ossdbs.electrodes import (ElectrodeModel,
                               AbbottStJudeActiveTipModel,
                               AbbottStJudeDirectedModel,
                               BostonScientificVerciseModel,
                               BostonScientificVerciseDirectedModel,
                               MedtronicModel,
                               MicroProbesRodentElectrodeModel,
                               MicroProbesSNEX100Model,
                               PINSMedicalModel,
                               AbbottStJudeParameters,
                               BostonScientificVerciseParameters,
                               MedtronicParameters,
                               MicroProbesRodentElectrodeParameters,
                               MicroProbesSNEX100Parameters,
                               PINSMedicalParameters)


class CustomElectrodeFactory:
    """Create an Electrode object using custom parameters."""

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
               parameters: dict,
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
            parameters = \
                AbbottStJudeParameters(parameters['TipLength[mm]'],
                                       parameters['ContactLength[mm]'],
                                       parameters['ContactSpacing[mm]'],
                                       parameters['LeadDiameter[mm]'],
                                       parameters['TotalLength[mm]']
                                       )

        elif electrode_type in (BostonScientificVerciseModel,
                                BostonScientificVerciseDirectedModel):
            parameters = \
                BostonScientificVerciseParameters(parameters['TipLength[mm]'],
                                                  parameters['ContactLength[mm]'],
                                                  parameters['ContactSpacing[mm]'],
                                                  parameters['LeadDiameter[mm]'],
                                                  parameters['TotalLength[mm]']
                                                  )

        elif electrode_type == MedtronicModel:
            parameters = \
                MedtronicParameters(parameters['TipLength[mm]'],
                                    parameters['ContactLength[mm]'],
                                    parameters['ContactSpacing[mm]'],
                                    parameters['LeadDiameter[mm]'],
                                    parameters['TotalLength[mm]']
                                    )

        elif electrode_type == PINSMedicalModel:
            parameters = \
                PINSMedicalParameters(parameters['TipLength[mm]'],
                                      parameters['ContactLength[mm]'],
                                      parameters['ContactSpacing[mm]'],
                                      parameters['LeadDiameter[mm]'],
                                      parameters['TotalLength[mm]']
                                      )

        elif electrode_type == MicroProbesRodentElectrodeModel:
            parameters = \
                MicroProbesRodentElectrodeParameters(parameters['TubeThickness[mm]'],
                                                     parameters['ContactLength[mm]'],
                                                     parameters['LeadDiameter[mm]'],
                                                     parameters['TotalLength[mm]']
                                                     )

        elif electrode_type == MicroProbesSNEX100Model:
            parameters = \
                MicroProbesSNEX100Parameters(parameters['CoreElectrodeDiameter[mm]'],
                                             parameters['CoreTubingDiameter[mm]'],
                                             parameters['CoreTubingLength[mm]'],
                                             parameters['CoreTubingDiameter[mm]'],
                                             parameters['OuterElectrodeLength[mm]'],
                                             parameters['OuterElectrodeDiameter[mm]'],
                                             parameters['OuterTubingDiameter[mm]'],
                                             parameters['TotalLength[mm]']
                                             )
        else:
            raise NotImplementedError("Electrode model not implemented.")

        return electrode_type(parameters=parameters,
                              direction=direction,
                              position=position,
                              rotation=rotation)
