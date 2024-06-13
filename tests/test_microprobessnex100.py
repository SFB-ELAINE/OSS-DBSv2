import netgen
import ngsolve
import numpy as np
import pytest

from ossdbs.electrodes import MicroProbesSNEX100


class TestMicroProbesSNEX_100:

    """
    FILE_PREFIX = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'test_data',
                               'MicroProbes',
                               'MicroProbesSNEX_100')

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
        electrode = MicroProbesSNEX100(rotation=rotation,
                                                   direction=direction,
                                                   position=position)
        geometry = electrode.geometry
        desired = self.load_geometry_data(path=path)
        assert desired == GeometryConverter(geometry).to_dictionary()

    def test_geometry_default(self):
        electrode = MicroProbesSNEX100()
        geometry = electrode.geometry
        desired = self.load_geometry_data(path=self.FILE_PREFIX+'_0.json')
        assert desired == GeometryConverter(geometry).to_dictionary()
    """

    @pytest.fixture
    def MicroProbesSNEX100_electrode(self):
        return MicroProbesSNEX100()

    def test_rename_boundaries(self, MicroProbesSNEX100_electrode):
        """Test whether set_contact_names() works."""
        electrode = MicroProbesSNEX100_electrode
        electrode.set_contact_names(
            {
                "Body": "RenamedBody",
                "Contact_1": "RenamedContact_1",
                "NonExistingPart": "NonExistingPart",
            }
        )
        geometry = electrode.geometry
        netgen_geometry = netgen.occ.OCCGeometry(geometry)
        with ngsolve.TaskManager():
            mesh = ngsolve.Mesh(netgen_geometry.GenerateMesh())
        desired = {"RenamedBody", "RenamedContact_1", "Contact_2"}
        assert desired == set(mesh.GetBoundaries())

    def test_contacts(self, MicroProbesSNEX100_electrode):
        """Test the number and names of contacts."""
        electrode = MicroProbesSNEX100_electrode
        geometry = electrode.geometry
        netgen_geometry = netgen.occ.OCCGeometry(geometry)
        with ngsolve.TaskManager():
            mesh = ngsolve.Mesh(netgen_geometry.GenerateMesh())
        desired = {
            "Body",
            "Contact_1",
            "Contact_2",
        }
        assert desired == set(mesh.GetBoundaries())

    def test_electrode_volume(self, MicroProbesSNEX100_electrode):
        """Test volume of the entire electrode."""
        electrode = MicroProbesSNEX100_electrode

        radius_1 = electrode._parameters.core_tubing_diameter * 0.5
        height_1 = electrode._parameters.core_tubing_length

        distance_2 = (
            electrode._parameters.core_electrode_length
            + electrode._parameters.core_tubing_length
            + electrode._parameters.outer_electrode_length
        )

        radius_2 = electrode._parameters.outer_tubing_diameter * 0.5
        height_2 = electrode._parameters.total_length - distance_2

        body1 = np.pi * radius_1**2 * height_1
        body2 = np.pi * radius_2**2 * height_2

        outer_electrode_radius = electrode._parameters.outer_electrode_diameter * 0.5
        outer_electrode_length = electrode._parameters.outer_electrode_length

        core_electrode_radius = electrode._parameters.core_electrode_diameter * 0.5
        core_electrode_length = electrode._parameters.core_electrode_length

        contact_1 = np.pi * outer_electrode_radius**2 * outer_electrode_length
        contact_2 = (4 / 3 * np.pi * core_electrode_radius**3 * 0.5) + (
            np.pi
            * core_electrode_radius**2
            * (core_electrode_length - core_electrode_radius)
        )

        desired = body1 + body2 + contact_1 + contact_2
        actual = electrode.geometry.mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)

    def test_contacts_volume(self, MicroProbesSNEX100_electrode):
        """Test volume of all the contacts."""
        electrode = MicroProbesSNEX100_electrode

        outer_electrode_radius = electrode._parameters.outer_electrode_diameter * 0.5
        outer_electrode_length = electrode._parameters.outer_electrode_length

        core_electrode_radius = electrode._parameters.core_electrode_diameter * 0.5
        core_electrode_length = electrode._parameters.core_electrode_length

        contact_1 = np.pi * outer_electrode_radius**2 * outer_electrode_length
        contact_2 = (4 / 3 * np.pi * core_electrode_radius**3 * 0.5) + (
            np.pi
            * core_electrode_radius**2
            * (core_electrode_length - core_electrode_radius)
        )

        desired = contact_1 + contact_2
        actual = electrode._contacts().mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)


"""
class TestMicroProbesSNEX_100_Capsule():

    FILE_PREFIX = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'test_data',
                               'MicroProbes',
                               'MicroProbesSNEX_100_Capsule')

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
        electrode = MicroProbesSNEX100(rotation=0.0,
                                                   direction=direction,
                                                   position=position)
        geometry = electrode.encapsulation_geometry(thickness=thickness)
        desired = self.load_geometry_data(path=path)
        assert desired == GeometryConverter(geometry).to_dictionary()
"""
