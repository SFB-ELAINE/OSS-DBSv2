import netgen
import ngsolve
import numpy as np
import pytest

from ossdbs.electrodes import BostonScientificVerciseDirected


class TestBostonScientificVerciseDirected:

    """
    FILE_PREFIX = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'test_data',
                               'BostonScientific',
                               'BostonScientificVerciseDirected')

    TESTDATA = [
        # electrode_parameters (Rotation, Position, Direction), file_path
        ((0.0, (0, 0, 0), (0, 0, 1)), FILE_PREFIX + '_0.json'),
        ((0.0, (0, 0, 0), (0, 0, 0)), FILE_PREFIX + '_0.json'),
        ((30.0, (0, 0, 0), (0, 0, 1)), FILE_PREFIX + '_1.json'),
        ((0.0, (1, -2, 3), (0, 0, 1)), FILE_PREFIX + '_2.json'),
        ((0.0, (1, -2, 3), (0, 0, 0)), FILE_PREFIX + '_2.json'),
        ((30.0, (1, -2, 3), (0, 0, 0)), FILE_PREFIX + '_3.json'),
        ((0.0, (1, -2, 3), (2.0, 0, 1.0)), FILE_PREFIX+'_4.json'),
        ((0.0, (1, -2, 3), (2.0/3.0, 0, 1.0/3.0)), FILE_PREFIX+'_4.json'),
        ]

    def load_geometry_data(self, path: str) -> dict:
        with open(path, "r") as file:
            geometry_data = json.load(file)
        return geometry_data

    @pytest.mark.parametrize('electrode_parameters, path', TESTDATA)
    def test_geometry(self, electrode_parameters, path) -> None:
        rotation, position, direction = electrode_parameters
        electrode = BostonScientificVerciseDirected(rotation=rotation,
                                                   direction=direction,
                                                   position=position)
        geometry = electrode.geometry
        desired = self.load_geometry_data(path=path)
        assert desired == GeometryConverter(geometry).to_dictionary()

    def test_geometry_default(self):
        electrode = BostonScientificVerciseDirected()
        geometry = electrode.geometry
        desired = self.load_geometry_data(path=self.FILE_PREFIX+'_0.json')
        assert desired == GeometryConverter(geometry).to_dictionary()
    """

    @pytest.fixture
    def BostonScientificVerciseDirected_electrode(self):
        return BostonScientificVerciseDirected()

    def test_rename_boundaries(self, BostonScientificVerciseDirected_electrode):
        """Test whether set_contact_names() works."""
        electrode = BostonScientificVerciseDirected_electrode
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
        desired = {
            "RenamedBody",
            "RenamedContact_1",
            "Contact_2",
            "Contact_3",
            "Contact_4",
            "Contact_5",
            "Contact_6",
            "Contact_7",
            "Contact_8",
        }
        assert desired == set(mesh.GetBoundaries())

    def test_contacts(self, BostonScientificVerciseDirected_electrode):
        """Test the number and names of contacts."""
        electrode = BostonScientificVerciseDirected_electrode
        geometry = electrode.geometry
        netgen_geometry = netgen.occ.OCCGeometry(geometry)
        with ngsolve.TaskManager():
            mesh = ngsolve.Mesh(netgen_geometry.GenerateMesh())
        desired = {
            "Body",
            "Contact_1",
            "Contact_2",
            "Contact_3",
            "Contact_4",
            "Contact_5",
            "Contact_6",
            "Contact_7",
            "Contact_8",
        }
        assert desired == set(mesh.GetBoundaries())

    def test_electrode_volume(self, BostonScientificVerciseDirected_electrode):
        """Test volume of the entire electrode."""
        electrode = BostonScientificVerciseDirected_electrode

        total_length = electrode._parameters.total_length
        tip_length = electrode._parameters.tip_length
        radius = electrode._parameters.lead_diameter * 0.5
        height = total_length - tip_length

        desired = (np.pi * radius**2 * height) + (4 / 3 * np.pi * radius**3 * 0.5)
        actual = electrode.geometry.mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)

    def test_contacts_volume(self, BostonScientificVerciseDirected_electrode):
        """Test volume of all the contacts."""
        electrode = BostonScientificVerciseDirected_electrode

        contact_length = electrode._parameters.contact_length
        radius = electrode._parameters.lead_diameter * 0.5

        C1_height = electrode._parameters.tip_length - radius
        C1_volume = (4 / 3 * radius**3 * np.pi * 0.5) + (
            C1_height * radius**2 * np.pi
        )

        desired = (
            C1_volume
            + (np.pi * radius**2 * contact_length)
            + (np.pi * radius**2 * contact_length * 90 / 360) * 6
        )
        actual = electrode._contacts().mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)


"""
class TestBostonScientificVerciseDirected_Capsule():

    FILE_PREFIX = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'test_data',
                               'BostonScientific',
                               'BostonScientificVerciseDirected_Capsule')

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
        electrode = BostonScientificVerciseDirected(rotation=0.0,
                                                   direction=direction,
                                                   position=position)
        geometry = electrode.encapsulation_geometry(thickness=thickness)
        desired = self.load_geometry_data(path=path)
        assert desired == GeometryConverter(geometry).to_dictionary()
"""
