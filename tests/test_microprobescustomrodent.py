from ossdbs.electrodes import MicroProbesRodentElectrode
import pytest
import netgen
import ngsolve
import numpy as np


class TestMicroProbesCustomRodent():

    """
    FILE_PREFIX = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'test_data',
                               'MicroProbes',
                               'MicroProbesCustomRodent')

    TESTDATA = [
        # electrode_parameters (Rotation, Position, Direction), file_path
        ((0.0, (0, 0, 0), (0, 0, 1)), FILE_PREFIX + '_0.json'),
        ((0.0, (0, 0, 0), (0, 0, 0)), FILE_PREFIX + '_0.json'),
        ((30.0, (0, 0, 0), (0, 0, 1)), FILE_PREFIX + '_0.json'),
        ((0.0, (1, -2, 3), (0, 0, 1)), FILE_PREFIX + '_1.json'),
        ((0.0, (1, -2, 3), (0, 0, 0)), FILE_PREFIX + '_1.json'),
        ((30.0, (1, -2, 3), (0, 0, 0)), FILE_PREFIX + '_1.json'),
        ((0.0, (1, -2, 3), (2.0, 0, 1.0)), FILE_PREFIX+'_2.json'),
        ((0.0, (1, -2, 3), (2.0/3.0, 0, 1.0/3.0)), FILE_PREFIX+'_2.json'),
        ]

    def load_geometry_data(self, path: str) -> dict:
        with open(path, "r") as file:
            geometry_data = json.load(file)
        return geometry_data

    @pytest.mark.parametrize('electrode_parameters, path', TESTDATA)
    def test_geometry(self, electrode_parameters, path) -> None:
        rotation, position, direction = electrode_parameters
        electrode = MicroProbesRodentElectrode(rotation=rotation,
                                                   direction=direction,
                                                   position=position)
        geometry = electrode.geometry
        desired = self.load_geometry_data(path=path)
        assert desired == GeometryConverter(geometry).to_dictionary()

    def test_geometry_default(self):
        electrode = MicroProbesRodentElectrode()
        geometry = electrode.geometry
        desired = self.load_geometry_data(path=self.FILE_PREFIX+'_0.json')
        assert desired == GeometryConverter(geometry).to_dictionary()
    """

    @pytest.fixture
    def MicroProbesRodentElectrode_electrode(self):
        return MicroProbesRodentElectrode()

    # Test whether set_contact_names() works
    def test_rename_boundaries(self, MicroProbesRodentElectrode_electrode):
        electrode = MicroProbesRodentElectrode_electrode
        electrode.set_contact_names({'Body': 'RenamedBody',
                                     'Contact_1': 'RenamedContact_1',
                                     'NonExistingPart': 'NonExistingPart'})
        geometry = electrode.geometry
        netgen_geometry = netgen.occ.OCCGeometry(geometry)
        with ngsolve.TaskManager():
            mesh = ngsolve.Mesh(netgen_geometry.GenerateMesh())
        desired = set(['RenamedBody',
                       'RenamedContact_1',
                       ])
        assert desired == set(mesh.GetBoundaries())

    # Test the number and names of contacts
    def test_contacts(self, MicroProbesRodentElectrode_electrode):
        electrode = MicroProbesRodentElectrode_electrode
        geometry = electrode.geometry
        netgen_geometry = netgen.occ.OCCGeometry(geometry)
        with ngsolve.TaskManager():
            mesh = ngsolve.Mesh(netgen_geometry.GenerateMesh())
        desired = set(['Body',
                       'Contact_1',
                       ])
        assert desired == set(mesh.GetBoundaries())

    # Test volume of the entire electrode
    def test_electrode_volume(self, MicroProbesRodentElectrode_electrode):
        electrode = MicroProbesRodentElectrode_electrode

        total_length = electrode._parameters.total_length
        lead_radius = electrode._parameters.lead_radius
        height = total_length - lead_radius

        desired = (height * lead_radius ** 2 * np.pi) + (4 / 3 * np.pi * lead_radius ** 3 * 0.5)
        actual = electrode.geometry.mass
        tolerance = 1e-3

        np.testing.assert_allclose(actual, desired, atol=tolerance)

    # Test volume of all the contacts
    def test_contacts_volume(self, MicroProbesRodentElectrode_electrode):
        electrode = MicroProbesRodentElectrode_electrode

        lead_radius = electrode._parameters.lead_radius
        contact_radius = electrode._parameters.contact_radius
        height = contact_radius - lead_radius

        desired = (height * lead_radius ** 2 * np.pi) + (4 / 3 * np.pi * lead_radius ** 3 * 0.5)
        actual = electrode._contacts().mass
        tolerance = 1e-3

        np.testing.assert_allclose(actual, desired, atol=tolerance)


"""
class TestMicroProbesCustomRodent_Capsule():

    FILE_PREFIX = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'test_data',
                               'MicroProbes',
                               'MicroProbesCustomRodent_Capsule')

    TESTDATA = [
        # electrode_parameters (Thickness, Position, Direction), file_path
        ((1.0, (0, 0, 0), (0, 0, 1)), FILE_PREFIX + '_0.json'),
        ((2.0, (0, 0, 0), (0, 0, 1)), FILE_PREFIX + '_1.json'),
        ((1.0, (1, -2, 3), (0, 0, 1)), FILE_PREFIX + '_2.json'),
        ]

    def load_geometry_data(self, path: str) -> dict:
        with open(path, "r") as file:
            geometry_data = json.load(file)
        return geometry_data

    @pytest.mark.parametrize('electrode_parameters, path', TESTDATA)
    def test_encapsulation_geometry(self, electrode_parameters, path) -> None:
        thickness, position, direction = electrode_parameters
        electrode = MicroProbesRodentElectrode(rotation=0.0,
                                                   direction=direction,
                                                   position=position)
        geometry = electrode.encapsulation_geometry(thickness=thickness)
        desired = self.load_geometry_data(path=path)
        assert desired == GeometryConverter(geometry).to_dictionary()
"""
