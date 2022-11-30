from src.brain_imaging.diffusion_tensor_imaging import DiffusionTensorImage
from src.brain_imaging.magnetic_resonance_imaging import MagneticResonanceImage
from src.electrodes import AbbottStjudeActiveTip6142_6145
from src.electrodes import AbbottStjudeActiveTip6146_6149
from src.electrodes import AbbottStjudeDirected6172
from src.electrodes import BostonScientificVercise
from src.electrodes import Medtronic3387, Medtronic3389, Medtronic3391
from src.electrodes import MicroProbesSNEX_100
from src.electrodes import PINSMedicalL301
from src.electrodes import PINSMedicalL302
from src.electrodes import PINSMedicalL303
from src.electrodes import MicroProbesCustomRodent
from src.signals import SignalGenerator


class Configuration:

    ELECTRODES = {'AbbottStjudeActiveTip6142_6145':
                  AbbottStjudeActiveTip6142_6145,
                  'AbbottStjudeActiveTip6146_6149':
                  AbbottStjudeActiveTip6146_6149,
                  'AbbottStjudeDirected6172':
                  AbbottStjudeDirected6172,
                  'BostonScientificVercise':
                  BostonScientificVercise,
                  'Medtronic3387':
                  Medtronic3387,
                  'Medtronic3389':
                  Medtronic3389,
                  'Medtronic3391':
                  Medtronic3391,
                  'MicroProbesSNEX_100':
                  MicroProbesSNEX_100,
                  'PINSMedicalL301':
                  PINSMedicalL301,
                  'PINSMedicalL302':
                  PINSMedicalL302,
                  'PINSMedicalL303':
                  PINSMedicalL303,
                  'MicroProbesCustomRodent':
                  MicroProbesCustomRodent
                  }

    def __init__(self, input: dict) -> None:
        self.input = input

    def electrodes(self) -> list:
        return [self.__create_electrode(par, idx)
                for idx, par in enumerate(self.input['Electrodes'])]

    def __create_electrode(self, parameters, index):
        electrode_class = self.ELECTRODES[parameters['Name']]
        electrode = electrode_class(direction=parameters['Direction'],
                                    translation=parameters['Translation'],
                                    rotation=parameters['Rotation'])
        names = {'Contact_{}'.format(i+1): "E{}C{}".format(index, i)
                 for i in range(len(parameters['Contacts']['Active']))}
        names.update({'Body': 'E{}B'.format(index)})
        electrode.rename_boundaries(names)
        return electrode

    def stimulation_signal(self):
        signal_type = self.input['StimulationSignal']['Type']
        frequency = self.input['StimulationSignal']['Frequency']
        pulse_width = self.input['StimulationSignal']['PulseWidthPercentage']
        top_width = self.input['StimulationSignal']['TopWidthPercentage']
        return SignalGenerator(signal_type=signal_type,
                               frequency=frequency,
                               pulse_width=pulse_width,
                               top_width=top_width).generate()

    def boundary_values(self):
        boundary_values = {}
        for index, electrode in enumerate(self.input['Electrodes']):
            boundary_values.update(self.__electrode_values(index, electrode))

        if self.input['BrainSurface']['Active']:
            value = self.input['BrainSurface_value']['Value']
            boundary_values.update({'Brain': value})

        return boundary_values

    def __electrode_values(self, index, electrode):
        values = {'E{}C{}'.format(index, i): value
                  for i, value in enumerate(electrode['Contacts']['Value'])
                  if electrode['Contacts']['Active'][i]}

        if electrode['Body']['Active']:
            values.update({'E{}B'.format(index): electrode['Body']['Value']})
        return values

    def magnetic_resonance_image(self):
        path = self.input['MagneticResonanceImage']['Path']
        return MagneticResonanceImage(file_path=path)

    def diffusion_tensor_image(self):
        path = self.input['DiffusionTensorImage']['Path']
        return DiffusionTensorImage(file_path=path)

    def output_path(self):
        return self.input['Output_directoy']
