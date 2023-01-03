
from dataclasses import dataclass
import json
import numpy as np
from ossdbs.brain_imaging.magnetic_resonance_imaging import MagneticResonanceImage
from ossdbs.brainsubstance import Material
from ossdbs.electrodes import AbbottStjudeActiveTip6142_6145
from ossdbs.electrodes import AbbottStjudeActiveTip6146_6149
from ossdbs.electrodes import AbbottStjudeDirected6172
from ossdbs.electrodes import BostonScientificVercise
from ossdbs.electrodes import BostonScientificVerciseDirected
from ossdbs.electrodes import Medtronic3387, Medtronic3389, Medtronic3391
from ossdbs.electrodes import MicroProbesSNEX_100
from ossdbs.electrodes import PINSMedicalL301
from ossdbs.electrodes import PINSMedicalL302
from ossdbs.electrodes import PINSMedicalL303
from ossdbs.electrodes import MicroProbesCustomRodent
from ossdbs.signals import RectangleSignal, TrapzoidSignal, TriangleSignal, Signal
from ossdbs.volume_conductor_model import VolumeConductorQS, VolumeConductorEQS


@dataclass
class Region:
    start: tuple = (0, 0, 0)
    end: tuple = (0, 0, 0)


class Input:

    def __init__(self, json_path: str) -> None:
        self.__input = self.__load_json(path=json_path)
        mri_path = self.__input['MagneticResonanceImage']['Path']
        mri = MagneticResonanceImage(mri_path)
        self.__offset = np.multiply(mri.bounding_box()[0], -1)

    @staticmethod
    def __load_json(path):
        with open(path, 'r') as json_file:
            return json.load(json_file)

    def mesh_order(self):
        return self.__input['MeshElementOrder']

    def mri(self):
        coding = self.__input['MagneticResonanceImage']['MaterialCoding']
        mri_coding = {Material.GRAY_MATTER: coding['GrayMatter'],
                      Material.WHITE_MATTER: coding['WhiteMatter'],
                      Material.CSF: coding['CerebrospinalFluid'],
                      Material.UNKNOWN: coding['Unknown']}
        mri_path = self.__input['MagneticResonanceImage']['Path']
        mri = MagneticResonanceImage(mri_path, mri_coding)
        mri.set_offset(self.__offset)
        return mri

    def electrodes(self):
        self.__shift_electrodes()
        return ElectrodeGenerator.electrodes(self.__input['Electrodes'])

    def __shift_electrodes(self):
        start = self.__offset
        for index in range(len(self.__input['Electrodes'])):
            translation = self.__input['Electrodes'][index]['Translation']
            new_translation = np.add(translation, start)
            self.__input['Electrodes'][index]['Translation'] = new_translation

    def boundary_values(self):
        boundaries = {
            'Electrodes': [{'Contacts': electrode['Contacts'],
                            'Body': electrode['Body']}
                           for electrode in self.__input['Electrodes']],
            'BrainSurface': self.__input['BrainSurface']}
        return BoundaryGenerator.generate(boundaries)

    def stimulation_signal(self):
        return SignalGenerator.generate(self.__input['StimulationSignal'])

    def output_path(self):
        return self.__input['OutputPath']

    def volume_conductor_type(self):
        if self.__input['FEMMode'] == 'EQS':
            return VolumeConductorEQS
        return VolumeConductorQS

    def region_of_interest(self):

        if not self.__input['RegionOfInterest']['Active']:
            mri_start, mri_end = self.mri().bounding_box()
            return Region(start=mri_start, end=mri_end)

        shape = self.__input['RegionOfInterest']['Shape']
        center = self.__input['RegionOfInterest']['Center'] + self.__offset
        start = center - np.divide(shape, 2)
        end = start + shape
        return Region(start=tuple(start.astype(int)),
                      end=tuple(end.astype(int)))


class ElectrodeGenerator:

    ELECTRODES = {'AbbottStjudeActiveTip6142_6145':
                  AbbottStjudeActiveTip6142_6145,
                  'AbbottStjudeActiveTip6146_6149':
                  AbbottStjudeActiveTip6146_6149,
                  'AbbottStjudeDirected6172':
                  AbbottStjudeDirected6172,
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

    @classmethod
    def electrodes(cls, electrodes) -> list:
        return [cls.__create_electrode(parameters, idx)
                for idx, parameters in enumerate(electrodes)]

    @classmethod
    def __create_electrode(cls, parameters, index):
        electrode_class = cls.ELECTRODES[parameters['Name']]
        electrode = electrode_class(direction=tuple(parameters['Direction']),
                                    translation=tuple(parameters
                                    ['Translation']),
                                    rotation=parameters['Rotation'])
        names = {'Contact_{}'.format(i+1): "E{}C{}".format(index, i)
                 for i in range(len(parameters['Contacts']['Active']))}
        names.update({'Body': 'E{}B'.format(index)})
        electrode.rename_boundaries(names)
        return electrode


class BoundaryGenerator:

    @classmethod
    def generate(cls, boundaries) -> dict:
        boundary_values = {}

        for index, electrode in enumerate(boundaries['Electrodes']):
            boundary_values.update(cls.__electrode_values(index, electrode))

        if boundaries['BrainSurface']['Active']:
            value = boundaries['Brainsurface']['Value']
            boundary_values.update({'Brain': value})

        return boundary_values

    @staticmethod
    def __electrode_values(index, electrode):
        values = {'E{}C{}'.format(index, i): value
                  for i, value in enumerate(electrode['Contacts']['Value'])
                  if electrode['Contacts']['Active'][i]}
        if electrode['Body']['Active']:
            values.update({'E{}B'.format(index): electrode['Body']['Value']})
        return values


class SignalGenerator:

    SIGNALS = {'Rectangle': RectangleSignal,
               'Triangle': TriangleSignal,
               'Trapzoid': TrapzoidSignal
               }

    @classmethod
    def generate(cls, parameters) -> Signal:
        signal_type = parameters['Type']
        frequency = parameters['Frequency']
        pulse_width = parameters['PulseWidthPercentage']
        top_width = parameters['TopWidthPercentage']

        if signal_type == 'Trapzoid':
            return TrapzoidSignal(frequency, pulse_width, top_width)
        signal = cls.SIGNALS[signal_type]
        return signal(frequency, pulse_width)
