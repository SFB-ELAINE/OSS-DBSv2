from ossdbs.electrodes import (ElectrodeModel,
                               AbbottStJudeActiveTip6142_6145,
                               AbbottStJudeActiveTip6146_6149,
                               AbbottStJudeDirected6172,
                               AbbottStJudeDirected6173,
                               BostonScientificVercise,
                               BostonScientificVerciseDirected,
                               Medtronic3387,
                               Medtronic3389,
                               Medtronic3391,
                               PINSMedicalL301,
                               PINSMedicalL302,
                               PINSMedicalL303,
                               MicroProbesRodentElectrode,
                               MicroProbesSNEX100)


class ElectrodeFactory:
    """Creates an Electrode object using the default geometry.

    See also
    --------

    :class:`ossdbs.electrodes.ElectrodeModel`

    """

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

    @classmethod
    def create(cls,
               name: str,
               direction: tuple,
               position: tuple,
               rotation: float
               ) -> ElectrodeModel:
        """Creates an Electrode object.

        Parameters
        ----------

        name: str
            Name of the electrode type

        Returns
        -------

        :class:`ossdbs.electrodes.ElectrodeModel`

        Notes
        -----

        TODO define default choices?

        """

        electrode_type = cls.ELECTRODES[name]
        return electrode_type(direction=direction,
                              position=position,
                              rotation=rotation)
