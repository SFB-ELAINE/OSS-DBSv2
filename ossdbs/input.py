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
from ossdbs.signals import RectangleSignal, TrapzoidSignal, TriangleSignal
from ossdbs.spectrum_modes import Octavevands, NoTruncationTest, SpectrumMode
from typing import List
import json
import numpy as np
import ngsolve


class Input:
    """Transform the input json.

    Parameters
    ----------
    json_path : str
    """

    def __init__(self, json_path: str) -> None:
        self.__input = self.__load_json(path=json_path)
        mri_path = self.__input['MaterialDistribution']['MRIPath']
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
        complex_datatpye = self.__input['EQSMode']
        return Mesh(ngsolve_mesh=ngsolve_mesh,
                    order=order,
                    complex_datatype=complex_datatpye)

    def conductivity(self) -> Conductivity:
        """Return the conductivity.

        Returns
        -------
        Conductivity
            Conductivity distribution in a given space.
        """
        coding = self.__input['MaterialDistribution']['MaterialCoding']
        mri_coding = {Material.GRAY_MATTER: coding['GrayMatter'],
                      Material.WHITE_MATTER: coding['WhiteMatter'],
                      Material.CSF: coding['CerebrospinalFluid'],
                      Material.UNKNOWN: coding['Unknown']}
        mri_path = self.__input['MaterialDistribution']['MRIPath']
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

        boundaries = {}
        for index, electrode in enumerate(self.__input['Electrodes']):
            boundaries.update(self.__active_electrode_values(index, electrode))

        case_grounding = self.__input['CaseGrounding']
        if case_grounding['Active']:
            boundaries.update({'BrainSurface': case_grounding['Value']})

        return boundaries

    def __active_electrode_values(self, electrode_index, electrode):
        return {'E{}C{}'.format(electrode_index, index):
                electrode['Contact_' + str(index+1)]['Value']
                for index in range(8)
                if electrode['Contact_' + str(index+1)]['Active']}

    def stimulation_signal(self):
        """Return stimulation signal.

        Returns
        -------
        Signal
        """

        parameters = self.__input['StimulationSignal']
        signal_type = parameters['Type']
        frequency = parameters['Frequency[Hz]']
        pulse_width = parameters['PulseWidth[µs]'] * frequency
        top_width = parameters['PulseTopWidth[µs]'] * frequency

        if signal_type == 'Trapzoid':
            return TrapzoidSignal(frequency, pulse_width, top_width)

        signal = {'Rectangle': RectangleSignal,
                  'Triangle': TriangleSignal,
                  'Trapzoid': TrapzoidSignal
                  }[signal_type]
        return signal(frequency, pulse_width)

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

        shape = (self.__input['RegionOfInterest']['Shape']['x[mm]'],
                 self.__input['RegionOfInterest']['Shape']['y[mm]'],
                 self.__input['RegionOfInterest']['Shape']['z[mm]'])
        center = (self.__input['RegionOfInterest']['Center']['x[mm]'],
                  self.__input['RegionOfInterest']['Center']['y[mm]'],
                  self.__input['RegionOfInterest']['Center']['z[mm]'])
        start = center - np.divide(shape, 2) + self.__offset
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
            par = self.__input['Electrodes'][index]['Translation']
            new_translation = {'x[mm]': par['x[mm]'] + self.__offset[0],
                               'y[mm]': par['y[mm]'] + self.__offset[1],
                               'z[mm]': par['z[mm]'] + self.__offset[2]}
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
        direction = (parameters['Direction']['x[mm]'],
                     parameters['Direction']['y[mm]'],
                     parameters['Direction']['z[mm]'])
        translation = (parameters['Translation']['x[mm]'],
                       parameters['Translation']['y[mm]'],
                       parameters['Translation']['z[mm]'])
        rotation = parameters['Rotation[Degrees]']
        electrode = elec_class(direction=direction,
                               translation=translation,
                               rotation=rotation)
        names = {'Contact_{}'.format(number+1): "E{}C{}".format(index, number)
                 for number in range(8)}
        names.update({'Body': 'E{}B'.format(index)})
        electrode.rename_boundaries(names)
        return electrode
