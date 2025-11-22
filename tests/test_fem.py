import json
import os

import pytest

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
        os.getcwd(), "input_files/homogeneous.nii.gz"
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
