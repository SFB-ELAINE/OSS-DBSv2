from ossdbs.brain_geometry import BrainGeometry
from ossdbs.brain_imaging.mri import MagneticResonanceImage
from ossdbs.materials import Material
from ossdbs.conductivity import Conductivity
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
from ossdbs.electrodes import Electrode
from ossdbs.mesh import Mesh
from ossdbs.region import Region
from ossdbs.signals import Signal
from ossdbs.signals import RectangleSignal, TrapzoidSignal, TriangleSignal
from ossdbs.spectrum_modes import Octavevands, NoTruncationTest, SpectrumMode
from typing import List
import json
import numpy as np
import ngsolve
import netgen


class Input:
    """Transform the input json.

    Parameters
    ----------
    json_path : str
    """

    def __init__(self, json_path: str) -> None:
        self.__input = self.__load_json(path=json_path)
        mri_path = self.__input['MagneticResonanceImage']['Path']
        mri = MagneticResonanceImage(mri_path)
        self.__offset = np.multiply(mri.bounding_box()[0], -1)
        self.__shift_electrodes()

    def mesh(self):
        electrodes = ElectrodeFactory.create(self.__input['Electrodes'])
        geometry = BrainGeometry(region=self.__region_of_interest(),
                                 electrodes=electrodes)
        netgen_geometry = geometry.netgen_geometry()

        if self.__input["Mesh"]["LoadMesh"]:
            file_path = self.__input["Mesh"]["LoadPath"]
            ngsolve_mesh = ngsolve.Mesh(filename=file_path)
            ngsolve_mesh.ngmesh.SetGeometry(netgen_geometry)
        else:
            ngsolve_mesh = ngsolve.Mesh(ngmesh=netgen_geometry.GenerateMesh())

        order = self.__input["Mesh"]["MeshElementOrder"]
        mesh = Mesh(ngsolve_mesh=ngsolve_mesh, order=order)
        mesh.set_complex(state=self.__input['FEMMode'] == 'EQS')
        return mesh

    def conductivity(self) -> Conductivity:
        """Return the conductivity.

        Returns
        -------
        Conductivity
            Conductivity distribution in a given space.
        """
        coding = self.__input['MagneticResonanceImage']['MaterialCoding']
        mri_coding = {Material.GRAY_MATTER: coding['GrayMatter'],
                      Material.WHITE_MATTER: coding['WhiteMatter'],
                      Material.CSF: coding['CerebrospinalFluid'],
                      Material.UNKNOWN: coding['Unknown']}
        mri_path = self.__input['MagneticResonanceImage']['Path']
        mri = MagneticResonanceImage(mri_path, mri_coding)
        mri.set_offset(self.__offset)
        return Conductivity(mri)

    def boundary_values(self):
        """Return the boundary values.

        Returns
        -------
        dict
            Boundary names and the associated values.
        """

        boundaries = {
            'Electrodes': [{'Contacts': electrode['Contacts']}
                           for electrode in self.__input['Electrodes']],
            'BrainSurface': self.__input['BrainSurface']}
        return BoundaryFactory.create_boundaries(boundaries)

    def stimulation_signal(self):
        """Return stimulation signal.

        Returns
        -------
        Signal
        """
        return SignalFactory.create_signal(self.__input['StimulationSignal'])

    def output_path(self) -> str:
        """Return path for results.

        Returns
        -------
        str
            Directory for result files.
        """
        return self.__input['OutputPath']

    def __region_of_interest(self) -> Region:
        """Return the region of interest.

        Returns
        -------
        Region
        """
        if not self.__input['RegionOfInterest']['Active']:
            mri_start, mri_end = self.mri().bounding_box()
            return Region(start=mri_start, end=mri_end)

        shape = self.__input['RegionOfInterest']['Shape']
        center = self.__input['RegionOfInterest']['Center'] + self.__offset
        start = center - np.divide(shape, 2)
        end = start + shape
        return Region(start=tuple(start.astype(int)),
                      end=tuple(end.astype(int)))

    def spectrum_mode(self) -> SpectrumMode:
        """Return the spectrum mode for the FEM.

        Returns
        -------
        SpectrumMode
        """
        return {'NoTruncation': NoTruncationTest(),
                'OctaveBand': Octavevands()
                }[self.__input['SpectrumMode']]

    @staticmethod
    def __load_json(path) -> dict:
        with open(path, 'r') as json_file:
            return json.load(json_file)

    def __shift_electrodes(self) -> None:
        for index in range(len(self.__input['Electrodes'])):
            translation = self.__input['Electrodes'][index]['Translation']
            new_translation = np.add(translation, self.__offset)
            self.__input['Electrodes'][index]['Translation'] = new_translation


class ElectrodeFactory:
    """Creates a list of Electrode objects."""

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
    def create(cls, electrodes: dict) -> List[Electrode]:
        """create a list of Electrode objects.

        Parameters
        ----------
        electrodes : dict

        Returns
        -------
        list
            Collection of electrode objects.
        """
        return [cls.__create_electrode(parameters, idx)
                for idx, parameters in enumerate(electrodes)]

    @classmethod
    def __create_electrode(cls, parameters: dict, index: int) -> Electrode:
        elec_class = cls.ELECTRODES[parameters['Name']]
        electrode = elec_class(direction=tuple(parameters['Direction']),
                               translation=tuple(parameters['Translation']),
                               rotation=parameters['Rotation'])
        names = {'Contact_{}'.format(number+1): "E{}C{}".format(index, number)
                 for number in range(len(parameters['Contacts']['Active']))}
        names.update({'Body': 'E{}B'.format(index)})
        electrode.rename_boundaries(names)
        return electrode


class BoundaryFactory:
    """Creates a dictionary of boundaries and corresponding boundary values."""

    @classmethod
    def create_boundaries(cls, boundaries: dict) -> dict:
        """Create a dictionary of boundaries and corresponding boundary
        values.

        Parameters
        ----------
        boundaries : dict

        Returns
        -------
        dict
        """
        boundary_values = {}

        for index, electrode in enumerate(boundaries['Electrodes']):
            boundary_values.update(cls.__electrode_values(index, electrode))

        if boundaries['BrainSurface']['Active']:
            value = boundaries['Brainsurface']['Value']
            boundary_values.update({'Brain': value})

        return boundary_values

    @staticmethod
    def __electrode_values(electrode_index: int, electrode: dict):
        return {'E{}C{}'.format(electrode_index, contact_index): value
                for contact_index, value
                in enumerate(electrode['Contacts']['Value'])
                if electrode['Contacts']['Active'][contact_index]}


class SignalFactory:
    """Creates signal.

    Returns
    -------
    Signal
    """

    SIGNALS = {'Rectangle': RectangleSignal,
               'Triangle': TriangleSignal,
               'Trapzoid': TrapzoidSignal
               }

    @classmethod
    def create_signal(cls, parameters: dict) -> Signal:
        """Create a signal.

        Parameters
        ----------
        parameters : dict

        Returns
        -------
        dict
        """

        signal_type = parameters['Type']
        frequency = parameters['Frequency']
        pulse_width = parameters['PulseWidthMicroSeconds'] * frequency
        top_width = parameters['TopWidthMicroSeconds'] * frequency

        if signal_type == 'Trapzoid':
            return TrapzoidSignal(frequency, pulse_width, top_width)
        signal = cls.SIGNALS[signal_type]
        return signal(frequency, pulse_width)
