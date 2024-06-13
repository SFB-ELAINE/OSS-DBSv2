import netgen
import ngsolve
import numpy as np
import pytest

from ossdbs.electrodes import MicroElectrode


class TestMicroElectrode:
    @pytest.fixture
    def MicroElectrode_electrode(self):
        return MicroElectrode()

    def test_rename_boundaries(self, MicroElectrode_electrode):
        """Test whether set_contact_names() works."""
        electrode = MicroElectrode_electrode
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
        desired = {"RenamedBody", "RenamedContact_1", "fillet"}
        assert desired == set(mesh.GetBoundaries())

    def test_contacts(self, MicroElectrode_electrode):
        """Test the number and names of contacts."""
        electrode = MicroElectrode_electrode
        geometry = electrode.geometry
        netgen_geometry = netgen.occ.OCCGeometry(geometry)
        with ngsolve.TaskManager():
            mesh = ngsolve.Mesh(netgen_geometry.GenerateMesh())
        desired = {"Body", "Contact_1", "fillet"}
        assert desired == set(mesh.GetBoundaries())

    def test_electrode_volume(self, MicroElectrode_electrode):
        """Test volume of the entire electrode."""
        electrode = MicroElectrode_electrode

        total_length = electrode._parameters.total_length
        contact_length = electrode._parameters.contact_length

        tip_length = electrode._parameters.tip_length
        tip_radius = electrode._parameters.tip_diameter * 0.5

        lead_radius = electrode._parameters.lead_diameter * 0.5
        height = total_length - tip_length
        filet_val = 0.001142187023875807

        desired = (
            (contact_length * tip_radius**2 * np.pi)
            + (height * lead_radius**2 * np.pi)
            - filet_val
        )
        actual = electrode.geometry.mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)

    def test_contacts_volume(self, MicroElectrode_electrode):
        """Test volume of all the contacts."""
        electrode = MicroElectrode_electrode

        contact_length = electrode._parameters.contact_length
        tip_radius = electrode._parameters.tip_diameter * 0.5
        filet_val = 0.001142187023875807

        desired = (contact_length * tip_radius**2 * np.pi) - filet_val
        actual = electrode._contacts().mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)
