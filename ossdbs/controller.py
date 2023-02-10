from ossdbs.electrode_contacts import ElectrodeContact, ContactCollection
from ossdbs.brain_geometry import BrainGeometry
from ossdbs.brain_imaging.mri import MagneticResonanceImage
from ossdbs.materials import Material
from ossdbs.conductivity import Conductivity
from ossdbs.electrodes import Electrode, ElectrodeParameters, ElectrodeFactory
from ossdbs.mesh import Mesh
from ossdbs.preconditioner import BDDCPreconditioner, LocalPreconditioner
from ossdbs.preconditioner import MultigridPreconditioner
from ossdbs.region import Region
from ossdbs.stimulation_signal import RectangleSignal, TrapzoidSignal, TriangleSignal
from ossdbs.solver import CGSolver, GMRESSolver
from typing import List
import json
import numpy as np
import ngsolve
from ossdbs.volume_conductor import VolumeConductor
from ossdbs.volume_conductor import VolumeConductorFloating
from ossdbs.volume_conductor import VolumeConductorFloatingImpedance
from ossdbs.volume_conductor import VolumeConductorNonFloating
from ossdbs.modes.spectrum import SpectrumMode
from ossdbs.modes.octave_band import OctaveBandModeImpedance, OctaveBandMode
from ossdbs.modes.no_truncation import NoTruncation, NoTruncationImpedance


class Controller:
    """Transform the input json.

    Parameters
    ----------
    json_path : str
    """

    def __init__(self, json_path: str) -> None:
        self.__input = self.__load_json(path=json_path)

    def mesh(self):
        electrodes = self.__create_electrodes()
        geometry = BrainGeometry(self.region_of_interest(), electrodes)
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
        return Conductivity(mri=mri, complex_datatype=self.__input['EQSMode'])

    def contacts(self) -> List[ElectrodeContact]:
        contacts = ContactCollection()
        for index, electrode in enumerate(self.__input['Electrodes']):
            for contact_index in range(0, 8):
                parameters = electrode['Contact_{}'.format(contact_index + 1)]
                contact = ElectrodeContact()
                contact.name = 'E{}C{}'.format(index, contact_index)
                contact.active = parameters['Active']
                contact.floating = parameters['Floating'] and not contact.active
                contact.current = parameters['Current[A]']
                contact.voltage = parameters['Voltage[V]']
                real = parameters['SurfaceImpedance[Ωm]']['real']
                imag = parameters['SurfaceImpedance[Ωm]']['imag']
                contact.surface_impedance = real + 1j * imag
                contacts.append(contact)

        contacts.append(contact)

        return contacts

    def output_path(self) -> str:
        """Return path for results.

        Returns
        -------
        str
            Directory for result files.
        """
        return self.__input['OutputPath']

    def volume_conductor(self) -> VolumeConductor:

        floating = self.__input['Floating']

        if not floating['Active']:
            return VolumeConductorNonFloating

        if not floating['FloatingImpedance']:
            return VolumeConductorFloating

        return VolumeConductorFloatingImpedance

    def region_of_interest(self) -> Region:
        """Return the region of interest.

        Returns
        -------
        Region
        """
        mri_path = self.__input['MaterialDistribution']['MRIPath']
        mri = MagneticResonanceImage(mri_path)
        shape = mri.xyz_shape()
        step_size = mri.voxel_size()

        size = (self.__input['RegionOfInterest']['Shape']['x[mm]'] * 1e-3,
                self.__input['RegionOfInterest']['Shape']['y[mm]'] * 1e-3,
                self.__input['RegionOfInterest']['Shape']['z[mm]'] * 1e-3)
        center = (self.__input['RegionOfInterest']['Center']['x[mm]'] * 1e-3,
                  self.__input['RegionOfInterest']['Center']['y[mm]'] * 1e-3,
                  self.__input['RegionOfInterest']['Center']['z[mm]'] * 1e-3)
        start = center - np.divide(size, 2)
        end = start + size

        x_step = step_size[0]
        y_step = step_size[1]
        z_step = step_size[2]

        x, y, z = size
        shape = int(x/x_step), int(y/y_step), int(z/z_step)
        return Region(tuple(start), tuple(end), shape, (x_step, y_step, z_step))

    def spectrum_mode(self) -> SpectrumMode:
        """Return the spectrum mode for the FEM.

        Returns
        -------
        SpectrumMode
        """

        return OctaveBandModeImpedance()

        if self.__input['SpectrumMode'] == 'OctaveBand':
            return OctaveBandMode()

        return NoTruncation()

    def stimulation_signal(self):
        """Return stimulation signal.

        Returns
        -------
        Signal
        """

        parameters = self.__input['StimulationSignal']
        frequency = parameters['Frequency[Hz]']
        pulse_width = parameters['PulseWidth[µs]'] * 1e-6 * frequency
        counter_width = parameters['CounterPulseWidth[µs]'] * 1e-6 * frequency
        inter_width = parameters['InterPulseWidth[µs]'] * 1e-6 * frequency
        top_width = parameters['PulseTopWidth[µs]'] * 1e-6 * frequency

        if parameters['Type'] == 'Trapzoid':
            return TrapzoidSignal(frequency=frequency,
                                  pulse_width=pulse_width,
                                  top_width=top_width,
                                  counter_pulse_width=counter_width,
                                  inter_pulse_width=inter_width)

        if parameters['Type'] == 'Triangle':
            return TriangleSignal(frequency=frequency,
                                  pulse_width=pulse_width,
                                  counter_pulse_width=counter_width,
                                  inter_pulse_width=inter_width)

        return RectangleSignal(frequency=frequency,
                               pulse_width=pulse_width,
                               counter_pulse_width=counter_width,
                               inter_pulse_width=inter_width)

    def coordinates(self):
        return np.array(self.__input['Points']) * 1e-3

    def __create_electrodes(self) -> List[Electrode]:

        electrodes = []
        for input_par in self.__input['Electrodes']:
            direction = (input_par['Direction']['x[mm]'] * 1e-3,
                         input_par['Direction']['y[mm]'] * 1e-3,
                         input_par['Direction']['z[mm]'] * 1e-3)
            position = (input_par['TipPosition']['x[mm]'] * 1e-3,
                        input_par['TipPosition']['y[mm]'] * 1e-3,
                        input_par['TipPosition']['z[mm]'] * 1e-3)
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

    @staticmethod
    def __load_json(path) -> dict:
        with open(path, 'r') as json_file:
            return json.load(json_file)

    def solver(self):
        solver_type = self.__input['Solver']['Type']
        solver = {'CG': CGSolver, 'GMRES': GMRESSolver}[solver_type]

        preconditioner_type = self.__input['Solver']['Preconditioner']
        preconditioner = {'bddc': BDDCPreconditioner(),
                          'local': LocalPreconditioner(),
                          'multigrid': MultigridPreconditioner()
                          }[preconditioner_type]

        return solver(precond_par=preconditioner,
                      printrates=self.__input['Solver']['PrintRates'],
                      maxsteps=self.__input['Solver']['MaximumSteps'],
                      precision=self.__input['Solver']['Precision'])
