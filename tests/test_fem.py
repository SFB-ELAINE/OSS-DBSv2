import json
import math
import os
from copy import deepcopy

import ngsolve
import numpy as np
import pytest
from netgen.occ import Box, Cylinder, OCCGeometry, Pnt, Z

import ossdbs
from ossdbs.fem.mesh import Mesh
from ossdbs.fem.preconditioner import (
    AMGPreconditioner,
    BDDCPreconditioner,
    CustomizedLocalPreconditioner,
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
            CustomizedLocalPreconditioner,
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


def _build_volume_conductor_with_solver(
    settings: dict,
    solver,
):
    """Build a volume conductor for a given solver on the stable case-1 setup."""
    local_settings = deepcopy(settings)
    local_settings["MaterialDistribution"]["MRIPath"] = os.path.join(
        os.getcwd(), "input_files/sub-John_Doe/JD_segmask.nii.gz"
    )
    local_settings["MaterialDistribution"]["DiffusionTensorActive"] = True
    local_settings["MaterialDistribution"]["DTIPath"] = os.path.join(
        os.getcwd(), "input_files/sub-John_Doe/JD_DTI_NormMapping.nii.gz"
    )
    mri_image, _ = ossdbs.load_images(local_settings)
    brain_region = ossdbs.create_bounding_box(local_settings["BrainRegion"])
    electrodes = ossdbs.generate_electrodes(local_settings)
    brain = ossdbs.BrainGeometry(local_settings["BrainRegion"]["Shape"], brain_region)
    geometry = ossdbs.ModelGeometry(brain, electrodes)
    ossdbs.set_contact_and_encapsulation_layer_properties(local_settings, geometry)

    dielectric_model = ossdbs.prepare_dielectric_properties(local_settings)
    materials = local_settings["MaterialDistribution"]["MRIMapping"]
    conductivity = ossdbs.ConductivityCF(
        mri_image,
        brain_region,
        dielectric_model,
        materials,
        geometry.encapsulation_layers,
        complex_data=local_settings["EQSMode"],
    )

    floating_mode = geometry.get_floating_mode()
    volume_conductor_classes = {
        "Floating": VolumeConductorFloating,
        "FloatingImpedance": VolumeConductorFloatingImpedance,
        "NonFloating": VolumeConductorNonFloating,
    }
    volume_conductor_class = volume_conductor_classes.get(
        floating_mode, VolumeConductorNonFloating
    )
    return volume_conductor_class(
        geometry,
        conductivity,
        solver,
        local_settings["FEMOrder"],
        local_settings["Mesh"],
    )


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


class TestCustomizedLocalPreconditioner:
    def test_customized_local_cg_solver_runs(
        self, settings_fixture, mri_fixture, geometry_fixture
    ):
        mri_image, _ = mri_fixture
        brain_region, _, geometry = geometry_fixture
        dielectric_model = ossdbs.prepare_dielectric_properties(settings_fixture)
        materials = settings_fixture["MaterialDistribution"]["MRIMapping"]
        ossdbs.set_contact_and_encapsulation_layer_properties(
            settings_fixture, geometry
        )

        solver = CGSolver(
            precond_par=CustomizedLocalPreconditioner(),
            maxsteps=2000,
            precision=1e-10,
        )

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

        volume_conductor.compute_solution(10000.0)

        sol = volume_conductor._potential.vec.FV().NumPy()
        assert np.all(np.isfinite(sol)), (
            "Customized local preconditioner produced NaN/Inf values."
        )
        assert not np.allclose(sol, 0.0), (
            "Customized local preconditioner produced a trivial zero solution."
        )

    def test_customized_local_impedance_matches_native_local(
        self, settings_fixture
    ) -> None:
        local_solver = CGSolver(
            precond_par=LocalPreconditioner(),
            maxsteps=2000,
            precision=1e-10,
        )
        customized_solver = CGSolver(
            precond_par=CustomizedLocalPreconditioner(),
            maxsteps=2000,
            precision=1e-10,
        )

        local_volume_conductor = _build_volume_conductor_with_solver(
            settings_fixture, local_solver
        )
        customized_volume_conductor = _build_volume_conductor_with_solver(
            settings_fixture, customized_solver
        )

        local_volume_conductor.compute_solution(10000.0)
        customized_volume_conductor.compute_solution(10000.0)

        local_impedance = local_volume_conductor.compute_impedance()
        customized_impedance = customized_volume_conductor.compute_impedance()

        assert np.isfinite(local_impedance), "Native local impedance is not finite."
        assert np.isfinite(customized_impedance), (
            "Customized local impedance is not finite."
        )

        relative_difference = abs(customized_impedance - local_impedance) / abs(
            local_impedance
        )
        assert relative_difference < 1e-2, (
            "Customized local impedance deviates too much from native local "
            f"({relative_difference:.3%})."
        )


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
