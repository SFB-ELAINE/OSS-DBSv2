from ossdbs.Nifti1Image import Nifti1Image
from ossdbs.electrode_contacts import ElectrodeContact, ContactCollection
from ossdbs.brain_geometry import BrainGeometry
from ossdbs.materials import Material
from ossdbs.conductivity import Conductivity
from ossdbs.electrodes import Electrode
from ossdbs.electrode_creation import ElectrodeFactory
from ossdbs.mesh import Mesh
from ossdbs.preconditioner import BDDCPreconditioner, LocalPreconditioner
from ossdbs.preconditioner import MultigridPreconditioner
from ossdbs.region import BoundingBox
from ossdbs.stimulation_signal import RectangleSignal
from ossdbs.stimulation_signal import TrapzoidSignal
from ossdbs.stimulation_signal import TriangleSignal
from ossdbs.solver import CGSolver, GMRESSolver
from ossdbs.volume_conductor import VolumeConductor
from ossdbs.volume_conductor import VolumeConductorFloating
from ossdbs.volume_conductor import VolumeConductorFloatingImpedance
from ossdbs.volume_conductor import VolumeConductorNonFloating
from ossdbs.point_analysis.fourier_analysis.spectrum import SpectrumMode
from ossdbs.point_analysis.fourier_analysis.octave_band import OctaveBandMode
from ossdbs.point_analysis.fourier_analysis.full_spectrum import FullSpectrum
from typing import List
import numpy as np
import ngsolve
import netgen
import h5py


class Controller:
    """Transform the input json.

    Parameters
    ----------
    json_path : str
    """

    def __init__(self, paramters: dict) -> None:
        self.__input = paramters

    def mesh(self):
        electrodes = self.__create_electrodes()
        encapsulating_d = self.__input['EncapsulatingLayer']['Thickness[mm]']
        geometry = BrainGeometry(self.region_of_interest(),
                                 electrodes,
                                 encapsulating_d)
        contacts = self.contacts()
        edges = contacts.active() + contacts.floating()
        geometry.edges_for_finer_meshing(edges)
        netgen_geometry = geometry.netgen_geometry()

        if self.__input["Mesh"]["LoadMesh"]:
            file_path = self.__input["Mesh"]["LoadPath"]
            ngsolve_mesh = ngsolve.Mesh(filename=file_path)
            ngsolve_mesh.ngmesh.SetGeometry(netgen_geometry)
        else:
            parameters = self.meshing_parameters()
            ng_mesh = netgen_geometry.GenerateMesh(parameters)
            ngsolve_mesh = ngsolve.Mesh(ngmesh=ng_mesh)

        order = self.__input["Mesh"]["MeshElementOrder"]
        complex_datatpye = self.__input['EQSMode']
        return Mesh(ngsolve_mesh, order, complex_datatpye)

    def meshing_parameters(self) -> netgen.meshing.MeshingParameters:
        mesh_type = self.__input['Mesh']['MeshingHypothesis']['Type']
        max_h = self.__input['Mesh']['MeshingHypothesis']['MaxH']
        return {'Coarse': netgen.meshing.meshsize.coarse,
                'Fine': netgen.meshing.meshsize.fine,
                'Moderate': netgen.meshing.meshsize.moderate,
                'VeryCoarse': netgen.meshing.meshsize.very_coarse,
                'VeryFine': netgen.meshing.meshsize.very_fine,
                'Default': netgen.meshing.MeshingParameters(),
                'Custom': netgen.meshing.MeshingParameters(max_h=max_h)
                }[mesh_type]

    def conductivity(self) -> Conductivity:
        """Return the conductivity.

        Returns
        -------
        Conductivity
            Conductivity distribution in a given space.
        """
        mri_path = self.__input['MaterialDistribution']['MRIPath']

        nifti = Nifti1Image(file_path=mri_path)
        voxel_size = nifti.voxel_size()
        offset = nifti.offset()

        region = self.region_of_interest()
        start = (np.array(region.start) - offset) / voxel_size
        end = (np.array(region.end) - offset) / voxel_size

        start_indices = np.floor(start).astype(int)
        stop_indices = np.ceil(end).astype(int)

        x_s, y_s, z_s = start_indices
        x_e, y_e, z_e = stop_indices
        data = nifti.data_map()[x_s:x_e, y_s:y_e, z_s:z_e]
        bbox = (tuple(start_indices * voxel_size + offset),
                tuple(stop_indices * voxel_size + offset))

        encap_points = self.points_in_encapsulating_layer(voxel_size, bbox)
        point_indices = (encap_points - bbox[0]) / voxel_size

        for index in point_indices:
            x, y, z = index.astype(int)
            data[x, y, z] = Material(1)

        return Conductivity(material_distribution=data,
                            bounding_box=bbox,
                            complex_datatype=self.__input['EQSMode'])

    def contacts(self) -> ContactCollection:

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
        shape = self.__input['RegionOfInterest']['Shape']
        center = self.__input['RegionOfInterest']['Center']
        size = (shape['x[mm]'], shape['y[mm]'], shape['z[mm]'])
        center = (center['x[mm]'], center['y[mm]'], center['z[mm]'])
        start = center - np.divide(size, 2)
        end = start + size
        return BoundingBox(tuple(start), tuple(end))

    def spectrum_mode(self) -> SpectrumMode:
        """Return the spectrum mode for the FEM.

        Returns
        -------
        SpectrumMode
        """

        if self.__input['SpectrumMode'] == 'OctaveBand':
            return OctaveBandMode()

        return FullSpectrum()

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

    def points(self):
        with h5py.File(self.__input['Points'], "r") as file:
            points = [np.array(file[key]) for key in file.keys()]
        return np.concatenate(points, axis=0)

    def categories(self):
        with h5py.File(self.__input['Points'], "r") as file:
            categories = [(key, file[key].shape[0]) for key in file.keys()]
        return categories

    def __create_electrodes(self) -> List[Electrode]:

        electrodes = []
        for input_par in self.__input['Electrodes']:
            dir = input_par['Direction']
            pos = input_par['TipPosition']
            direction = (dir['x[mm]'], dir['y[mm]'], dir['z[mm]'])
            position = (pos['x[mm]'], pos['y[mm]'], pos['z[mm]'])
            rotation = input_par['Rotation[Degrees]']
            electrode = ElectrodeFactory.create(name=input_par['Name'],
                                                direction=direction,
                                                position=position,
                                                rotation=rotation)
            electrodes.append(electrode)

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

    def points_in_encapsulating_layer(self, voxel_size, bounding_box) -> np.ndarray:

        start_limit, end_limit = bounding_box
        electrodes = self.__create_electrodes()

        points = []

        for electrode in electrodes:

            geo = electrode.capsule_geometry(thickness=0.1)
            start, end = geo.bounding_box
            start = np.max([start_limit, tuple(start)], axis=0)
            end = np.min([end_limit, tuple(end)], axis=0)
            start = (start // voxel_size) * voxel_size
            end = ((end // voxel_size) + (end % voxel_size > 0)) * voxel_size

            x_coor = np.arange(start[0], end[0], voxel_size[0])
            y_coor = np.arange(start[1], end[1], voxel_size[1])
            z_coor = np.arange(start[2], end[2], voxel_size[2])

            points.append(np.array([(x, y, z)
                                    for x in x_coor
                                    for y in y_coor
                                    for z in z_coor]))

        points = np.concatenate(points)
        netgen_geometry = netgen.occ.OCCGeometry(geo)
        ng_mesh = netgen_geometry.GenerateMesh()
        ngsolve_mesh = ngsolve.Mesh(ngmesh=ng_mesh)
        mesh = Mesh(ngsolve_mesh=ngsolve_mesh, order=2)
        included = mesh.is_included(points)

        return points[included]
