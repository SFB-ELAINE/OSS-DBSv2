import json
import math
import os

import ngsolve
import numpy as np
import pytest
from netgen.occ import Box, Cylinder, OCCGeometry, Pnt, Z

import ossdbs
from ossdbs.fem.mesh import Mesh
from ossdbs.fem.preconditioner import (
    AMGPreconditioner,
    BDDCPreconditioner,
    DirectPreconditioner,
    LocalPreconditioner,
    MultigridPreconditioner,
)
from ossdbs.fem.solver import CGSolver, DirectSolver, GMRESSolver
from ossdbs.fem.volume_conductor.floating import VolumeConductorFloating
from ossdbs.fem.volume_conductor.floating_impedance import (
    VolumeConductorFloatingImpedance,
)
from ossdbs.fem.volume_conductor.nonfloating import VolumeConductorNonFloating
from ossdbs.utils.settings import Settings


class TestPreconditioner:
    @pytest.mark.parametrize(
        "preconditioner_class",
        [
            AMGPreconditioner,
            BDDCPreconditioner,
            DirectPreconditioner,
            LocalPreconditioner,
            MultigridPreconditioner,
        ],
    )
    def test_preconditioner(self, preconditioner_class):
        try:
            solver = preconditioner_class()
            assert solver is not None
        except Exception:
            pytest.fail("Cannot be instantiated.")


class TestSolver:
    @pytest.mark.parametrize("solver_class", [CGSolver, GMRESSolver, DirectSolver])
    def test_solver(self, solver_class):
        try:
            solver = solver_class(
                precond_par=BDDCPreconditioner(), maxsteps=10000, precision=1e-12
            )
            assert solver is not None
        except Exception:
            pytest.fail("Cannot be instantiated.")


@pytest.fixture
def settings_fixture():
    json_path = os.path.join(
        os.getcwd(), "input_test_cases/input_case1/input_homogeneous.json"
    )
    with open(json_path) as file:
        settings = json.load(file)
    settings = Settings(settings).complete_settings()
    return settings


@pytest.fixture
def mri_fixture(settings_fixture):
    settings_fixture["MaterialDistribution"]["MRIPath"] = os.path.join(
        os.getcwd(), "input_files/sub-John_Doe/JD_segmask.nii.gz"
    )
    settings_fixture["MaterialDistribution"]["DiffusionTensorActive"] = True
    settings_fixture["MaterialDistribution"]["DTIPath"] = os.path.join(
        os.getcwd(), "input_files/sub-John_Doe/JD_DTI_NormMapping.nii.gz"
    )
    mri_image, dti_image = ossdbs.load_images(settings_fixture)
    return mri_image, dti_image


@pytest.fixture
def geometry_fixture(settings_fixture):
    settings = settings_fixture
    brain_region = ossdbs.create_bounding_box(settings["BrainRegion"])
    electrodes = ossdbs.generate_electrodes(settings)
    brain = ossdbs.BrainGeometry(settings["BrainRegion"]["Shape"], brain_region)
    geometry = ossdbs.ModelGeometry(brain, electrodes)
    return brain_region, electrodes, geometry


class TestMesh:
    def test_mesh(self, geometry_fixture, settings_fixture):
        try:
            geometry = geometry_fixture[2].geometry

            mesh = Mesh(geometry, settings_fixture["FEMOrder"])
            assert mesh is not None
        except Exception:
            pytest.fail("Cannot be instantiated.")


class TestConductivity:
    def test_conductivityCF(self, settings_fixture, mri_fixture, geometry_fixture):
        try:
            mri_image, _ = mri_fixture
            brain_region, _, geometry = geometry_fixture
            dielectric_model = ossdbs.prepare_dielectric_properties(settings_fixture)
            materials = settings_fixture["MaterialDistribution"]["MRIMapping"]

            conductivity = ossdbs.ConductivityCF(
                mri_image,
                brain_region,
                dielectric_model,
                materials,
                geometry.encapsulation_layers,
                complex_data=settings_fixture["EQSMode"],
            )
            assert conductivity is not None
        except Exception:
            pytest.fail("Cannot be instantiated.")


class TestDTIMasking:
    def test_DTImasking(self, settings_fixture, mri_fixture, geometry_fixture):
        try:
            mri_image, dti_image = mri_fixture
            brain_region, _, geometry = geometry_fixture
            dielectric_model = ossdbs.prepare_dielectric_properties(settings_fixture)
            materials = settings_fixture["MaterialDistribution"]["MRIMapping"]

            # Create ConductivityCF instances
            conductivity_unmasked = ossdbs.ConductivityCF(
                mri_image,
                brain_region,
                dielectric_model,
                materials,
                geometry.encapsulation_layers,
                complex_data=settings_fixture["EQSMode"],
                dti_image=dti_image,
                wm_masking=False,  # No masking
            )

            conductivity_masked = ossdbs.ConductivityCF(
                mri_image,
                brain_region,
                dielectric_model,
                materials,
                geometry.encapsulation_layers,
                complex_data=settings_fixture["EQSMode"],
                dti_image=dti_image,
                wm_masking=True,  # With masking
            )

            # Generate mesh and get underlying NGSolve mesh
            mesh = ossdbs.generate_mesh(settings_fixture)
            ngmesh = mesh.ngsolvemesh

            # Build tensor-valued conductivity fields
            sigma_unmasked_cf = conductivity_unmasked(mesh=mesh, frequency=10000.0)
            sigma_masked_cf = conductivity_masked(mesh=mesh, frequency=10000.0)

            # Helper to evaluate sigma at a physical point and return a 3x3 tensor
            def eval_sigma(cf, point):
                mp = ngmesh(*point)  # unpack (x, y, z) -> ngmesh(x, y, z)
                vals = cf(mp)  # 9 components (flattened 3x3)
                return np.array(vals, dtype=float).reshape((3, 3))

            # Tissue-specific test points (updated by you)
            test_points = {
                "GM": (-11.2, -2.5, 5.2),
                "WM": (-10.9, -2.6, -2.3),
            }

            # Evaluate tensors
            sigma_gm_unmasked = eval_sigma(sigma_unmasked_cf, test_points["GM"])
            sigma_gm_masked = eval_sigma(sigma_masked_cf, test_points["GM"])

            sigma_wm_unmasked = eval_sigma(sigma_unmasked_cf, test_points["WM"])
            sigma_wm_masked = eval_sigma(sigma_masked_cf, test_points["WM"])

            # WM: masking should not change conductivity (only CSF/GM)
            assert np.allclose(
                sigma_wm_unmasked,
                sigma_wm_masked,
            ), (
                "Conductivity tensors in white matter should be identical with and without masking."
            )

            # GM: masking should modify anisotropic GM tensors (become isotropic)
            assert not np.allclose(
                sigma_gm_unmasked,
                sigma_gm_masked,
            ), (
                "Conductivity tensors in anisotropic gray matter should differ with masking applied."
            )

        except Exception as e:
            pytest.fail(f"Test failed with exception: {e}")


class TestVolumeConductorModel:
    def test_volume_conductor_model(
        self, settings_fixture, mri_fixture, geometry_fixture
    ):
        try:
            mri_image, _ = mri_fixture
            brain_region, _, geometry = geometry_fixture
            dielectric_model = ossdbs.prepare_dielectric_properties(settings_fixture)
            materials = settings_fixture["MaterialDistribution"]["MRIMapping"]
            solver = ossdbs.prepare_solver(settings_fixture)

            conductivity = ossdbs.ConductivityCF(
                mri_image,
                brain_region,
                dielectric_model,
                materials,
                geometry.encapsulation_layers,
                complex_data=settings_fixture["EQSMode"],
            )

            floating_mode = geometry.get_floating_mode()
            volume_conductor_classes = {
                "Floating": VolumeConductorFloating,
                "FloatingImpedance": VolumeConductorFloatingImpedance,
                "NonFloating": VolumeConductorNonFloating,
            }

            # If floating_mode is None, call VolumeConductorNonFloating.
            VolumeConductorClass = volume_conductor_classes.get(
                floating_mode, VolumeConductorNonFloating
            )
            volume_conductor = VolumeConductorClass(
                geometry,
                conductivity,
                solver,
                settings_fixture["FEMOrder"],
                settings_fixture["Mesh"],
            )
            assert volume_conductor is not None
        except Exception:
            pytest.fail("Cannot be instantiated.")


class TestHPRefineCurvedBoundaryIntegration:
    """Regression test for NGSolve bug in 6.2.2602:
    Integrate on curved boundaries returns NaN inside TaskManager
    after RefineHP + Curve applied outside TaskManager.
    """

    @pytest.fixture
    def hp_refined_mesh(self):
        box = Box(Pnt(-2, -2, -2), Pnt(2, 2, 2))
        cyl = Cylinder(Pnt(0, 0, -1), Z, r=0.5, h=2)
        cyl.faces.Max(Z).name = "contact_top"
        cyl.faces.Min(Z).name = "contact_bottom"
        box.faces.name = "outer"
        geo = OCCGeometry(box - cyl)

        with ngsolve.TaskManager():
            mesh = ngsolve.Mesh(geo.GenerateMesh(maxh=0.5))
            mesh.Curve(2)

        mesh.RefineHP(levels=2, factor=0.125)
        mesh.Curve(2)
        return mesh

    def test_curved_boundary_integration_not_nan(self, hp_refined_mesh):
        """Curved boundary integration inside TaskManager must not return NaN."""
        with ngsolve.TaskManager():
            val = ngsolve.Integrate(
                ngsolve.CF(1.0) * ngsolve.ds("contact_top"), hp_refined_mesh
            )
        assert not math.isnan(val), (
            "Integrate on curved boundary 'contact_top' returned NaN inside TaskManager. "
            "This is a known NGSolve 6.2.2602 regression."
        )

    def test_curved_boundary_area_accuracy(self, hp_refined_mesh):
        """Curved boundary area must be close to the analytical value."""
        expected_area = math.pi * 0.5**2  # disk of radius 0.5
        with ngsolve.TaskManager():
            val = ngsolve.Integrate(
                ngsolve.CF(1.0) * ngsolve.ds("contact_top"), hp_refined_mesh
            )
        assert abs(val - expected_area) < 0.02, (
            f"Integrate on 'contact_top' = {val}, expected ~{expected_area:.4f}"
        )

    def test_flat_boundary_integration_inside_taskmanager(self, hp_refined_mesh):
        """Flat boundary integration inside TaskManager should work regardless."""
        expected_area = 4.0 * 4.0 * 6  # 6 faces of a 4x4x4 box
        with ngsolve.TaskManager():
            val = ngsolve.Integrate(
                ngsolve.CF(1.0) * ngsolve.ds("outer"), hp_refined_mesh
            )
        assert not math.isnan(val)
        assert abs(val - expected_area) < 0.1


class TestContactSurfaceImpedance:
    """Tests for Contact and Contacts surface impedance functionality."""

    def test_contact_get_surface_impedance_R_model(self):
        """Test Contact.get_surface_impedance with a resistor model."""
        from ossdbs.model_geometry.contacts import Contact

        contact = Contact(
            name="test_contact",
            area=1.0,  # 1 mm^2
            active=False,
            floating=True,
            surface_impedance_model="R",
            surface_impedance_parameters={"R": 100.0},
        )
        # R model: Z = R, independent of frequency
        z = contact.get_surface_impedance(frequency=1000.0, is_complex=False)
        assert np.isclose(z, 100.0, rtol=1e-6), f"Expected Z=100, got {z}"

    def test_contact_get_surface_impedance_scales_with_area(self):
        """Surface impedance should scale linearly with contact area."""
        from ossdbs.model_geometry.contacts import Contact

        R = 50.0
        area = 2.5
        contact = Contact(
            name="test",
            area=area,
            surface_impedance_model="R",
            surface_impedance_parameters={"R": R},
        )
        z = contact.get_surface_impedance(frequency=1000.0, is_complex=False)
        assert np.isclose(z, R * area, rtol=1e-6), f"Expected Z={R * area}, got {z}"

    def test_contact_get_surface_impedance_complex(self):
        """Test complex return for EQS mode with RC model."""
        from ossdbs.model_geometry.contacts import Contact

        contact = Contact(
            name="test",
            area=1.0,
            surface_impedance_model="R",
            surface_impedance_parameters={"R": 100.0},
        )
        z = contact.get_surface_impedance(frequency=1000.0, is_complex=True)
        # For pure R, imaginary part should be zero
        assert isinstance(z, complex)
        assert np.isclose(z.real, 100.0, rtol=1e-6)
        assert np.isclose(z.imag, 0.0, atol=1e-10)

    def test_contacts_get_surface_impedances(self):
        """Test Contacts.get_surface_impedances returns dict with correct keys."""
        from ossdbs.model_geometry.contacts import Contact, Contacts

        contacts = Contacts(
            [
                Contact(
                    name="c1",
                    active=True,
                    floating=False,
                    area=1.0,
                    surface_impedance_model="R",
                    surface_impedance_parameters={"R": 100.0},
                ),
                Contact(name="c2", active=False, floating=True, area=1.0),
            ]
        )
        z_dict = contacts.get_surface_impedances(frequency=1000.0, is_complex=False)
        assert z_dict["c1"] is not None
        assert np.isclose(z_dict["c1"], 100.0)
        assert z_dict["c2"] is None


class TestPrepareCurrentControlledMode:
    """Tests for prepare_current_controlled_mode with various contact configs."""

    @pytest.fixture
    def floating_impedance_vc(self, settings_fixture, mri_fixture, geometry_fixture):
        """Create a FloatingImpedance volume conductor."""
        mri_image, _ = mri_fixture
        brain_region, _, geometry = geometry_fixture
        dielectric_model = ossdbs.prepare_dielectric_properties(settings_fixture)
        materials = settings_fixture["MaterialDistribution"]["MRIMapping"]
        solver = ossdbs.prepare_solver(settings_fixture)

        conductivity = ossdbs.ConductivityCF(
            mri_image,
            brain_region,
            dielectric_model,
            materials,
            geometry.encapsulation_layers,
            complex_data=settings_fixture["EQSMode"],
        )

        # Set up contacts for floating impedance mode:
        # 1 active (ground) + 1 floating with surface impedance
        for contact in geometry.contacts:
            if contact.active:
                contact.voltage = 0.0
                break
        # Make the second contact floating with impedance
        c2 = geometry.contacts[1]
        c2.active = False
        c2.floating = True
        c2.current = 0.001
        c2.surface_impedance_model = "R"
        c2.surface_impedance_parameters = {"R": 100.0}

        vc = VolumeConductorFloatingImpedance(
            geometry,
            conductivity,
            solver,
            settings_fixture["FEMOrder"],
            settings_fixture["Mesh"],
        )
        return vc

    def test_prepare_cc_mode_1_active_0_floating(
        self, settings_fixture, mri_fixture, geometry_fixture
    ):
        """2 active contacts should work (standard CC mode)."""
        mri_image, _ = mri_fixture
        brain_region, _, geometry = geometry_fixture
        dielectric_model = ossdbs.prepare_dielectric_properties(settings_fixture)
        materials = settings_fixture["MaterialDistribution"]["MRIMapping"]
        solver = ossdbs.prepare_solver(settings_fixture)

        conductivity = ossdbs.ConductivityCF(
            mri_image,
            brain_region,
            dielectric_model,
            materials,
            geometry.encapsulation_layers,
            complex_data=settings_fixture["EQSMode"],
        )

        vc = VolumeConductorNonFloating(
            geometry,
            conductivity,
            solver,
            settings_fixture["FEMOrder"],
            settings_fixture["Mesh"],
        )
        # Should not raise for 2 active contacts
        vc.prepare_current_controlled_mode()

    def test_prepare_cc_mode_3_active_raises(
        self, settings_fixture, mri_fixture, geometry_fixture
    ):
        """3+ active contacts should raise ValueError."""
        mri_image, _ = mri_fixture
        brain_region, _, geometry = geometry_fixture
        dielectric_model = ossdbs.prepare_dielectric_properties(settings_fixture)
        materials = settings_fixture["MaterialDistribution"]["MRIMapping"]
        solver = ossdbs.prepare_solver(settings_fixture)

        conductivity = ossdbs.ConductivityCF(
            mri_image,
            brain_region,
            dielectric_model,
            materials,
            geometry.encapsulation_layers,
            complex_data=settings_fixture["EQSMode"],
        )

        # Set 3 contacts to active (if enough contacts exist)
        active_count = 0
        for contact in geometry.contacts:
            if active_count < 3:
                contact.active = True
                contact.floating = False
                contact.voltage = 0.0
                active_count += 1

        if active_count >= 3:
            vc = VolumeConductorNonFloating(
                geometry,
                conductivity,
                solver,
                settings_fixture["FEMOrder"],
                settings_fixture["Mesh"],
            )
            with pytest.raises(ValueError):
                vc.prepare_current_controlled_mode()
