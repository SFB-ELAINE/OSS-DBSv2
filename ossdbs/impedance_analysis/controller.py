from ossdbs.electrode_contacts import ElectrodeContact, ContactCollection
from ossdbs.brain_geometry import BrainGeometry
from ossdbs.brain_imaging.mri import MagneticResonanceImage
from ossdbs.materials import Material
from ossdbs.conductivity import Conductivity
from ossdbs.electrodes import Electrode, ElectrodeParameters, ElectrodeFactory
from ossdbs.mesh import Mesh
from ossdbs.preconditioner import BDDCPreconditioner, LocalPreconditioner
from ossdbs.preconditioner import MultigridPreconditioner
from ossdbs.region import BoundingBox
from ossdbs.stimulation_signal import RectangleSignal
from ossdbs.solver import CGSolver, GMRESSolver
from typing import List
import numpy as np
import ngsolve
from ossdbs.volume_conductor import VolumeConductor
from ossdbs.volume_conductor import VolumeConductorFloating
from ossdbs.volume_conductor import VolumeConductorFloatingImpedance
from ossdbs.volume_conductor import VolumeConductorNonFloating
from ossdbs.impedance_analysis.fourier_analysis.spectrum import SpectrumMode
from ossdbs.impedance_analysis.fourier_analysis.logarithm_scanning \
    import LogarithmScanning


class Controller:
    """Transform the input json.

    Parameters
    ----------
    json_path : str
    """

    def __init__(self, parameters: dict) -> None:
        self.__input = parameters

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
            for contact_par in electrode['Contacts']:
                name = 'E{}C{}'.format(index, contact_par['Contact_ID'] - 1)
                active = contact_par['Active']
                floating = contact_par['Floating'] and not active
                current = contact_par['Current[A]']
                voltage = contact_par['Voltage[V]']
                real = contact_par['SurfaceImpedance[Ωm]']['real']
                imag = contact_par['SurfaceImpedance[Ωm]']['imag']
                surface_impedance = real + 1j * imag

                contact = ElectrodeContact(name=name,
                                           active=active,
                                           floating=floating,
                                           current=current,
                                           voltage=voltage,
                                           surface_impedance=surface_impedance)

                contacts.append(contact)

        acvtive_contacts = contacts.active()
        for contact_name in acvtive_contacts[2:]:
            contacts.set_inactive(contact_name)

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

    def region_of_interest(self) -> BoundingBox:
        """Return the region of interest.

        Returns
        -------
        Region
        """
        size = (self.__input['RegionOfInterest']['Shape']['x[mm]'] * 1e-3,
                self.__input['RegionOfInterest']['Shape']['y[mm]'] * 1e-3,
                self.__input['RegionOfInterest']['Shape']['z[mm]'] * 1e-3)
        center = (self.__input['RegionOfInterest']['Center']['x[mm]'] * 1e-3,
                  self.__input['RegionOfInterest']['Center']['y[mm]'] * 1e-3,
                  self.__input['RegionOfInterest']['Center']['z[mm]'] * 1e-3)
        start = center - np.divide(size, 2)
        end = start + size
        return BoundingBox(tuple(start), tuple(end))

    def spectrum_mode(self) -> SpectrumMode:
        """Return the spectrum mode for the FEM.

        Returns
        -------
        SpectrumMode
        """

        return LogarithmScanning()

    def stimulation_signal(self):
        """Return stimulation signal.

        Returns
        -------
        Signal
        """

        parameters = self.__input['StimulationSignal']
        return RectangleSignal(frequency=parameters['Frequency[Hz]'],
                               pulse_width=0.0,
                               counter_pulse_width=0.0,
                               inter_pulse_width=0.0)

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
