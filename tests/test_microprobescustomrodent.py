import netgen
import ngsolve
import numpy as np
import pytest

from ossdbs.electrodes import MicroProbesRodentElectrode


class TestMicroProbesCustomRodent:
    @pytest.fixture
    def MicroProbesRodentElectrode_electrode(self):
        return MicroProbesRodentElectrode()

    def test_rename_boundaries(self, MicroProbesRodentElectrode_electrode):
        """Test whether set_contact_names() works."""
        electrode = MicroProbesRodentElectrode_electrode
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
        }
        assert desired == set(mesh.GetBoundaries())

    def test_contacts(self, MicroProbesRodentElectrode_electrode):
        """Test the number and names of contacts."""
        electrode = MicroProbesRodentElectrode_electrode
        geometry = electrode.geometry
        netgen_geometry = netgen.occ.OCCGeometry(geometry)
        with ngsolve.TaskManager():
            mesh = ngsolve.Mesh(netgen_geometry.GenerateMesh())
        desired = {
            "Body",
            "Contact_1",
        }
        assert desired == set(mesh.GetBoundaries())

    def test_electrode_volume(self, MicroProbesRodentElectrode_electrode):
        """Test volume of the entire electrode."""
        electrode = MicroProbesRodentElectrode_electrode

        total_length = electrode._parameters.total_length
        lead_radius = electrode._parameters.lead_radius
        contact_radius = electrode._parameters.contact_radius
        height = total_length - contact_radius

        desired = (height * lead_radius**2 * np.pi) + (
            4 / 3 * np.pi * contact_radius**3 * 0.5
        )
        actual = electrode.geometry.mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)

    def test_contacts_volume(self, MicroProbesRodentElectrode_electrode):
        """Test volume of all the contacts."""
        electrode = MicroProbesRodentElectrode_electrode

        contact_radius = electrode._parameters.contact_radius

        desired = 4 / 3 * np.pi * contact_radius**3 * 0.5
        actual = electrode._contacts().mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)
