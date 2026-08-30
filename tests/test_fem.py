import json
import math
import os
import subprocess
import sys
import textwrap
from copy import deepcopy
from types import SimpleNamespace

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
from ossdbs.fem.solver import (
    CGSolver,
    DirectSolver,
    GMRESSolver,
    _finalize_krylov_solve,
    _krylov_convergence_status,
)
from ossdbs.fem.volume_conductor.floating import VolumeConductorFloating
from ossdbs.fem.volume_conductor.floating_impedance import (
    VolumeConductorFloatingImpedance,
)
from ossdbs.fem.volume_conductor.nonfloating import VolumeConductorNonFloating
from ossdbs.fem.volume_conductor.volume_conductor_model import VolumeConductor
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
                precond_par=BDDCPreconditioner(),
                maxsteps=10000,
                relative_tolerance=1e-12,
                absolute_tolerance=1e-14,
            )
            assert solver is not None
            assert solver._absolute_tolerance == 1e-14
        except Exception:
            pytest.fail("Cannot be instantiated.")

    def test_absolute_tolerance_from_settings(self):
        settings = Settings({"Solver": {"AbsoluteTolerance": 1e-8}}).complete_settings()

        solver = ossdbs.prepare_solver(settings)

        assert solver._absolute_tolerance == 1e-8

    def test_absolute_tolerance_default(self):
        settings = Settings({}).complete_settings()

        solver = ossdbs.prepare_solver(settings)

        assert np.isclose(solver._absolute_tolerance, 1e-12)

    def test_absolute_tolerance_must_be_positive(self):
        settings = Settings({"Solver": {"AbsoluteTolerance": 0.0}}).complete_settings()

        with pytest.raises(ValueError, match="AbsoluteTolerance must be positive"):
            ossdbs.prepare_solver(settings)

    @pytest.mark.parametrize(
        ("solver_class", "ngsolve_solver_name"),
        [(CGSolver, "CGSolver"), (GMRESSolver, "GMResSolver")],
    )
    def test_absolute_tolerance_forwarded_to_ngsolve(
        self, solver_class, ngsolve_solver_name, mocker
    ):
        ngsolve_solver = SimpleNamespace(Solve=mocker.Mock())
        constructor = mocker.patch(
            f"ossdbs.fem.solver.ngsolve.krylovspace.{ngsolve_solver_name}",
            return_value=ngsolve_solver,
        )
        mocker.patch("ossdbs.fem.solver.ngsolve.Preconditioner")
        residual = mocker.MagicMock()
        mocker.patch(
            "ossdbs.fem.solver._log_initial_system_state", return_value=residual
        )
        mocker.patch("ossdbs.fem.solver._get_debug_mode", return_value=(False, False))
        finalize = mocker.patch("ossdbs.fem.solver._finalize_krylov_solve")
        bilinear_form = mocker.MagicMock()
        linear_form = mocker.MagicMock()
        grid_function = mocker.MagicMock()

        solver = solver_class(
            precond_par=BDDCPreconditioner(),
            maxsteps=100,
            relative_tolerance=1e-8,
            absolute_tolerance=1e-10,
        )
        solver.bvp(bilinear_form, linear_form, grid_function)

        assert constructor.call_args.kwargs["tol"] == 1e-8
        assert constructor.call_args.kwargs["atol"] == 1e-10
        finalize.assert_called_once()

    @pytest.mark.parametrize(
        (
            "residuals",
            "relative_tolerance",
            "absolute_tolerance",
            "expected_target",
            "expected_converged",
        ),
        [
            ([2.0, 0.02], 1e-2, None, 0.02, True),
            ([2.0, 0.021], 1e-2, None, 0.02, False),
            ([2.0, 0.04], 1e-2, 0.05, 0.05, True),
            ([2.0, 0.06], None, 0.05, 0.05, False),
        ],
    )
    def test_krylov_convergence_status(
        self,
        residuals,
        relative_tolerance,
        absolute_tolerance,
        expected_target,
        expected_converged,
    ):
        solver = SimpleNamespace(
            residuals=residuals,
            tol=relative_tolerance,
            atol=absolute_tolerance,
        )

        initial, final, target, converged = _krylov_convergence_status(solver)

        assert initial == residuals[0]
        assert final == residuals[-1]
        assert target == pytest.approx(expected_target)
        assert converged is expected_converged

    def test_convergence_at_maximum_steps_does_not_raise(self, mocker):
        solver = SimpleNamespace(
            residuals=[1.0, 1e-8],
            tol=1e-8,
            atol=None,
            iterations=500,
        )
        mocker.patch("ossdbs.fem.solver._warn_local_preconditioner_issue")

        _finalize_krylov_solve(
            solver,
            mocker.MagicMock(),
            mocker.MagicMock(),
            mocker.MagicMock(),
            mocker.MagicMock(),
            mocker.MagicMock(),
            False,
            "CG",
            500,
        )

    def test_nonconvergence_message_reports_effective_target(self, mocker):
        solver = SimpleNamespace(
            residuals=[0.1, 4.448e-9],
            tol=1e-8,
            atol=None,
            iterations=500,
        )
        mocker.patch("ossdbs.fem.solver._warn_local_preconditioner_issue")

        with pytest.raises(RuntimeError) as error:
            _finalize_krylov_solve(
                solver,
                mocker.MagicMock(),
                mocker.MagicMock(),
                mocker.MagicMock(),
                mocker.MagicMock(),
                mocker.MagicMock(),
                False,
                "CG",
                500,
            )

        message = str(error.value)
        assert "Initial residual: 0.1" in message
        assert "relative tolerance: 1e-08" in message
        assert "absolute tolerance: None" in message
        assert "target residual: 1e-09" in message
        assert "final residual: 4.448e-09" in message

    def test_absolute_tolerance_accepts_previous_final_residual(self, mocker):
        solver = SimpleNamespace(
            residuals=[0.1, 4.448e-9],
            tol=1e-8,
            atol=1e-8,
            iterations=500,
        )
        mocker.patch("ossdbs.fem.solver._warn_local_preconditioner_issue")

        _finalize_krylov_solve(
            solver,
            mocker.MagicMock(),
            mocker.MagicMock(),
            mocker.MagicMock(),
            mocker.MagicMock(),
            mocker.MagicMock(),
            False,
            "CG",
            500,
        )


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

            # Tissue-specific test points in JD MRI frame
            test_points = {
                "GM": (1.9, -34.6, 4.4),
                "WM": (0.8, -33.7, -0.6),
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
            relative_tolerance=1e-10,
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
            relative_tolerance=1e-10,
        )
        customized_solver = CGSolver(
            precond_par=CustomizedLocalPreconditioner(),
            maxsteps=2000,
            relative_tolerance=1e-10,
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


class TestRefineAfterRefineHPIsBroken:
    """Document the NGSolve limitation that prevents combining AMR with
    hp-refinement.

    Conceptually we *want* to run AMR on top of an hp-refined mesh: hp
    handles geometry-driven refinement near electrode contacts, and AMR
    would then handle solution-driven refinement elsewhere. However,
    NGSolve's standard Mesh.Refine() (the primitive AMR drives) silently
    corrupts the mesh when called after RefineHP() — the integrated
    volume of the domain stops matching the true geometry. Because the
    failure is silent (no exception, no warning), VolumeConductor must
    proactively disable AMR whenever hp-refinement has been applied.
    See VolumeConductor._resolve_amr_active.
    """

    @staticmethod
    def _build_mesh():
        box = Box(Pnt(-2, -2, -2), Pnt(2, 2, 2))
        cyl = Cylinder(Pnt(0, 0, -1), Z, r=0.5, h=2)
        cyl.faces.Max(Z).name = "contact_top"
        cyl.faces.Min(Z).name = "contact_bottom"
        box.faces.name = "outer"
        geo = OCCGeometry(box - cyl)
        with ngsolve.TaskManager():
            mesh = ngsolve.Mesh(geo.GenerateMesh(maxh=1.0))
            mesh.Curve(2)
        return mesh

    @staticmethod
    def _true_volume():
        # 4x4x4 box minus a cylinder of radius 0.5 and height 2
        return 4.0 * 4.0 * 4.0 - math.pi * 0.5**2 * 2.0

    def test_refine_hp_alone_preserves_volume(self):
        """Sanity check: RefineHP on its own keeps the mesh volume correct."""
        mesh = self._build_mesh()
        mesh.RefineHP(levels=1, factor=0.25)
        mesh.Curve(2)
        with ngsolve.TaskManager():
            vol = ngsolve.Integrate(ngsolve.CF(1.0) * ngsolve.dx, mesh)
        assert abs(vol - self._true_volume()) < 1e-2, (
            f"RefineHP alone should preserve volume, got {vol} "
            f"vs expected {self._true_volume()}"
        )

    def test_standard_refine_after_refine_hp_corrupts_mesh(self):
        """Calling NGSolve's Refine() (the primitive used by AMR) after
        RefineHP() yields a mesh whose integrated volume no longer
        matches the true geometry. This is precisely the scenario we
        would need to support to combine AMR with hp-refinement, and it
        is the reason that combination is currently disabled in
        VolumeConductor._resolve_amr_active.

        On some platforms (observed on macOS) the corruption manifests
        as a hard crash inside ``mesh.Curve(2)`` rather than a silently
        broken volume. Either symptom proves the bug is still present,
        so we run the dangerous sequence in a subprocess: a non-zero
        exit code (segfault, abort, ...) and a finite-but-wrong volume
        both count as "bug still present"; only a correct volume should
        fail this test and prompt us to revisit the guard.
        """
        script = textwrap.dedent(
            """
            import math
            import ngsolve
            from netgen.occ import Box, Cylinder, OCCGeometry, Pnt, Z

            box = Box(Pnt(-2, -2, -2), Pnt(2, 2, 2))
            cyl = Cylinder(Pnt(0, 0, -1), Z, r=0.5, h=2)
            cyl.faces.Max(Z).name = "contact_top"
            cyl.faces.Min(Z).name = "contact_bottom"
            box.faces.name = "outer"
            geo = OCCGeometry(box - cyl)
            with ngsolve.TaskManager():
                mesh = ngsolve.Mesh(geo.GenerateMesh(maxh=1.0))
                mesh.Curve(2)
            mesh.RefineHP(levels=1, factor=0.25)
            mesh.Curve(2)
            for el in mesh.Elements():
                mesh.SetRefinementFlag(el, True)
            mesh.Refine()
            mesh.Curve(2)
            with ngsolve.TaskManager():
                vol = ngsolve.Integrate(ngsolve.CF(1.0) * ngsolve.dx, mesh)
            print(f"VOLUME={vol}")
            """
        )

        result = subprocess.run(
            [sys.executable, "-c", script],
            capture_output=True,
            text=True,
            timeout=300,
        )

        true_vol = self._true_volume()

        if result.returncode != 0:
            # Subprocess crashed (e.g. segfault) — bug is still present,
            # the AMR/hp guard is still required. Nothing more to check.
            return

        vol = None
        for line in result.stdout.splitlines():
            if line.startswith("VOLUME="):
                vol = float(line.split("=", 1)[1])
                break
        assert vol is not None, (
            "Could not parse volume from subprocess output:\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )

        # The mesh is broken: volume diverges by several percent.
        # If NGSolve ever fixes this, the assertion will fail and we
        # can revisit the AMR/hp exclusion in VolumeConductor.
        assert abs(vol - true_vol) > 1.0, (
            "NGSolve Refine() after RefineHP() unexpectedly preserved the "
            f"mesh volume ({vol} vs {true_vol}). The hp/AMR mutual-exclusion "
            "guard in VolumeConductor may no longer be necessary."
        )


class TestResolveAmrActiveGuard:
    """Verify the hp-refinement / AMR mutual-exclusion guard in
    VolumeConductor._resolve_amr_active.
    """

    @staticmethod
    def _call(hp_applied: bool, settings):
        """Invoke the unbound method on a minimal stub object.

        ``hp_refinement_applied`` is declared as a ``@property`` here to
        match the real ``ossdbs.fem.mesh.Mesh`` class — otherwise the
        stub would silently mask a regression where the production code
        calls the attribute as a method (``False()`` -> TypeError).
        """

        class _MeshStub:
            def __init__(self, hp):
                self._hp = hp

            @property
            def hp_refinement_applied(self) -> bool:
                return self._hp

        class _Stub:
            mesh = _MeshStub(hp_applied)

            # _resolve_amr_active calls self._check_AMR_settings; reuse
            # the real implementation so we exercise the same validation.
            _check_AMR_settings = VolumeConductorNonFloating._check_AMR_settings

        return VolumeConductorNonFloating._resolve_amr_active(_Stub(), settings)

    def test_mesh_exposes_hp_refinement_applied_as_property(self):
        """Guard against regressions that turn ``hp_refinement_applied``
        into a method. The production guard in ``_resolve_amr_active``
        reads it as an attribute, so the two must stay in sync.
        """
        assert isinstance(Mesh.hp_refinement_applied, property), (
            "Mesh.hp_refinement_applied must remain a @property — "
            "VolumeConductor._resolve_amr_active reads it as an attribute."
        )

    def test_settings_none_returns_false(self):
        assert self._call(hp_applied=False, settings=None) is False

    def test_active_true_without_hp_returns_true(self):
        settings = {"Active": True, "ErrorTolerance": 0.1, "MaxIterations": 1}
        assert self._call(hp_applied=False, settings=settings) is True

    def test_active_false_returns_false(self):
        settings = {"Active": False, "ErrorTolerance": 0.1, "MaxIterations": 1}
        assert self._call(hp_applied=False, settings=settings) is False

    def test_active_true_with_hp_is_disabled_and_warns(self, caplog):
        settings = {"Active": True, "ErrorTolerance": 0.1, "MaxIterations": 1}
        with caplog.at_level("WARNING", logger="ossdbs"):
            result = self._call(hp_applied=True, settings=settings)
        assert result is False
        assert any(
            "mutually exclusive" in rec.getMessage() for rec in caplog.records
        ), "Expected a warning about hp-refinement / AMR being mutually exclusive"

    def test_invalid_settings_still_validated(self):
        # Missing ErrorTolerance / MaxIterations must raise via _check_AMR_settings
        with pytest.raises(ValueError, match="ErrorTolerance"):
            self._call(hp_applied=False, settings={"Active": True})


class TestScalarImpedanceOnly:
    """Guard that the stimulation pipeline only deals in scalar impedance.

    The NxN admittance-matrix path was auto-triggered when
    ``n_active + n_floating > 2``, which routed CC + floating-impedance
    configurations through a matrix path that ``get_scale_factor`` and
    ``_store_solution_at_contacts`` could not consume. The full admittance
    matrix is now considered an analysis tool that lives outside
    ``run_full_analysis`` (see ``docs/impedance_analyzer_plan.md``).

    Invariants:
    1. ``run_full_analysis`` no longer exposes the ``multicontact_impedance``
       kwarg.
    2. ``get_scale_factor`` returns a plain scalar.
    3. ``compute_impedance`` raises ``NotImplementedError`` when called with
       anything other than exactly 2 active contacts.
    """

    def test_multicontact_kwarg_removed(self):
        import inspect

        from ossdbs.fem.volume_conductor.volume_conductor_model import (
            VolumeConductor,
        )

        sig = inspect.signature(VolumeConductor.run_full_analysis)
        assert "multicontact_impedance" not in sig.parameters, (
            "VolumeConductor.run_full_analysis must not expose "
            "``multicontact_impedance`` — the admittance-matrix mode is "
            "handled by a standalone analyzer."
        )

    def test_get_scale_factor_returns_scalar(self):
        """``get_scale_factor`` must return a plain scalar so
        ``_store_solution_at_contacts`` can assign into a single complex
        slot of ``_free_stimulation_variable``.
        """
        from ossdbs.fem.volume_conductor.volume_conductor_model import (
            VolumeConductor,
        )

        class _FakeContact:
            def __init__(self, current):
                self.current = current

        class _FakeContacts:
            def __init__(self, active):
                self.active = active

        class _Stub:
            current_controlled = True
            is_complex = True
            _impedances = np.array([10.0 + 0.0j, 12.0 + 1.0j], dtype=complex)
            contacts = _FakeContacts(
                [_FakeContact(1.0 + 0.0j), _FakeContact(-1.0 + 0.0j)]
            )

            impedances = VolumeConductor.impedances

        stub = _Stub()
        scale = VolumeConductor.get_scale_factor(stub, 0)
        assert np.ndim(scale) == 0, (
            f"get_scale_factor must return a plain scalar; got shape {np.shape(scale)}"
        )
        # Expected: Z * |I| = 10 * 1 = 10
        assert np.isclose(scale, 10.0 + 0.0j), f"Expected Z * |I| = 10+0j, got {scale}"

    def test_compute_impedance_raises_for_nontwo_active(self):
        """``compute_impedance`` only handles the 2-active scalar case.

        Multicontact configurations must go through the standalone
        admittance-matrix analyzer.
        """
        from ossdbs.fem.volume_conductor.volume_conductor_model import (
            VolumeConductor,
        )

        class _FakeContacts:
            def __init__(self, active):
                self.active = active

        class _Stub:
            contacts = _FakeContacts([object(), object(), object()])

        with pytest.raises(NotImplementedError, match="exactly 2 active"):
            VolumeConductor.compute_impedance(_Stub())


class _FakeContact:
    def __init__(self, name, current):
        self.name = name
        self.current = current
        self.voltage = 0.0


class _FakeContacts:
    def __init__(self, contacts):
        self._contacts = contacts
        self.active = contacts
        self.floating = []

    def __iter__(self):
        return iter(self._contacts)


class _SkipBookkeepingVolumeConductor(VolumeConductor):
    """Stub exercising only the octave-band skip-frequency bookkeeping."""

    def __init__(self, output_path):
        self._contacts = _FakeContacts(
            [_FakeContact("E0", 1.0), _FakeContact("E1", -1.0)]
        )
        self.is_complex = True
        self.output_path = output_path
        self._floating_potentials = None
        self.solved_frequencies = []

    def compute_solution(self, frequency):
        self.solved_frequencies.append(frequency)

    def update_space(self):
        pass

    def prepare_current_controlled_mode(self):
        pass

    def _save_report(self, timings):
        pass

    def compute_impedance(self):
        freq_idx = round(self.solved_frequencies[-1] / 100.0)
        return complex(100.0 + freq_idx, 0.0)

    def estimate_currents(self):
        freq_idx = round(self.solved_frequencies[-1] / 100.0)
        return {contact.name: 50.0 + freq_idx for contact in self.contacts}

    def _has_sigma_changed(self, computing_idx, frequency_indices, threshold):
        # Force exactly one skip, at the representative frequency 8 * 100 Hz,
        # so that computing_idx (4) and the previous representative
        # frequency's real array index (4) diverge from a naive
        # `computing_idx - 1` (3).
        return computing_idx != 4


class TestOctaveBandSkipBookkeeping:
    """Regression test for the octave-band "skip frequency" bookkeeping.

    ``run_full_analysis`` reuses the impedance/current of the previous
    *representative* frequency when the conductivity has not changed
    enough to warrant a new solve. Under octave-band approximation the
    loop counter (``computing_idx``, a position in the sparse
    ``frequency_indices`` array) and the actual frequency-array index
    (``freq_idx``) diverge, so the reused value must be looked up via
    ``frequency_indices[computing_idx - 1]`` rather than
    ``computing_idx - 1`` directly.
    """

    def _run(self, tmp_path):
        vc = _SkipBookkeepingVolumeConductor(str(tmp_path))
        signal = SimpleNamespace(
            frequencies=np.arange(20) * 100.0,
            amplitudes=np.ones(20),
            octave_band_approximation=True,
            current_controlled=True,
        )
        vc.run_full_analysis(signal, compute_impedance=True, estimate_currents=True)
        return vc

    def test_skipped_impedance_matches_previous_representative_frequency(
        self, tmp_path
    ):
        vc = self._run(tmp_path)
        # Band for freq_idx=8 (Hz 800) is [7, 8, 9, 10, 11]; it must copy
        # the impedance computed at freq_idx=4 (104), not the stale value
        # sitting at raw array position `computing_idx - 1 == 3` (102).
        assert np.allclose(vc.impedances[[7, 8, 9, 10, 11]], 104.0)
        assert not np.any(np.isnan(vc.impedances))

    def test_skipped_currents_match_previous_representative_frequency(self, tmp_path):
        vc = self._run(tmp_path)
        currents = vc._currents["E0"]
        assert np.allclose(currents[[7, 8, 9, 10, 11]], 54.0)
        assert not np.any(np.isnan(currents))
