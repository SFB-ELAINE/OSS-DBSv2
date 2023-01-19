from ossdbs.brain_geometry import BrainGeometry
from ossdbs.brain_imaging.mri import MagneticResonanceImage
from ossdbs.materials import Material
from ossdbs.conductivity import Conductivity
from ossdbs.electrodes import Electrode, ElectrodeParameters, ElectrodeFactory
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
        electrodes = self.__create_electrodes()
        geometry = BrainGeometry(self.__region_of_interest(), electrodes)
        netgen_geometry = geometry.netgen_geometry()

        if self.__input["Mesh"]["LoadMesh"]:
            file_path = self.__input["Mesh"]["LoadPath"]
            ngsolve_mesh = ngsolve.Mesh(filename=file_path)
            ngsolve_mesh.ngmesh.SetGeometry(netgen_geometry)
        else:
            ngsolve_mesh = ngsolve.Mesh(ngmesh=netgen_geometry.GenerateMesh())

        order = self.__input["Mesh"]["MeshElementOrder"]
        complex_datatpye = self.__input['EQSMode']
        return Mesh(ngsolve_mesh, order, complex_datatpye)

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
        return Conductivity(mri=mri)

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
            boundaries.update({'BrainSurface': case_grounding['Voltage[V]']})

        return boundaries

    def __active_electrode_values(self, electrode_index, electrode):
        return {'E{}C{}'.format(electrode_index, index):
                electrode['Contact_' + str(index+1)]['Voltage[V]']
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
        mri_path = self.__input['MaterialDistribution']['MRIPath']
        mri = MagneticResonanceImage(mri_path)
        mri.set_offset(self.__offset)
        mri_start, mri_end = mri.bounding_box()

        if not self.__input['RegionOfInterest']['Active']:
            return Region(start=mri_start, end=mri_end)

        shape = (self.__input['RegionOfInterest']['Shape']['x[mm]'] * 1e-3,
                 self.__input['RegionOfInterest']['Shape']['y[mm]'] * 1e-3,
                 self.__input['RegionOfInterest']['Shape']['z[mm]'] * 1e-3)
        center = (self.__input['RegionOfInterest']['Center']['x[mm]'] * 1e-3,
                  self.__input['RegionOfInterest']['Center']['y[mm]'] * 1e-3,
                  self.__input['RegionOfInterest']['Center']['z[mm]'] * 1e-3)
        start = center - np.divide(shape, 2) + self.__offset
        end = start + shape
        return Region(start=tuple(start), end=tuple(end))

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

    def __create_electrodes(self) -> List[Electrode]:

        electrodes = []
        for input_par in self.__input['Electrodes']:
            direction = (input_par['Direction']['x[mm]'] * 1e-3,
                         input_par['Direction']['y[mm]'] * 1e-3,
                         input_par['Direction']['z[mm]'] * 1e-3)
            position = (input_par['Position']['x[mm]'] * 1e-3,
                        input_par['Position']['y[mm]'] * 1e-3,
                        input_par['Position']['z[mm]'] * 1e-3)
            rotation = input_par['Rotation[Degrees]']
            electrode_par = ElectrodeParameters(name=input_par['Name'],
                                                direction=direction,
                                                position=position,
                                                rotation=rotation)
            electrodes.append(ElectrodeFactory.create(electrode_par))

        for index, electrode in enumerate(electrodes):
            boundary_names = {'Contact_1': "E{}C0".format(index),
                              'Contact_2': "E{}C1".format(index),
                              'Contact_3': "E{}C2".format(index),
                              'Contact_4': "E{}C3".format(index),
                              'Contact_5': "E{}C4".format(index),
                              'Contact_6': "E{}C5".format(index),
                              'Contact_7': "E{}C6".format(index),
                              'Contact_8': "E{}C7".format(index),
                              'Body': 'E{}B'.format(index)}
            electrode.rename_boundaries(boundary_names)

        return electrodes

    def __shift_electrodes(self) -> None:
        x, y, z = np.multiply(self.__offset, 1e3)
        for index in range(len(self.__input['Electrodes'])):
            position = self.__input['Electrodes'][index]['Position']
            new_position = {'x[mm]': position['x[mm]'] + x,
                            'y[mm]': position['y[mm]'] + y,
                            'z[mm]': position['z[mm]'] + z}
            self.__input['Electrodes'][index]['Position'] = new_position
