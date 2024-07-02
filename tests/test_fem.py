import pytest

import ossdbs
from ossdbs.fem.solver import CGSolver
from ossdbs.fem.volume_conductor.floating import VolumeConductorFloating


class TestSolver:
    def test_CGsolver(self):
        with pytest.raises(AttributeError):
            CGSolver(
                precond_par="Wrong Input",
                maxsteps="Wrong Input",
                precision="Wrong Input",
            )


class TestConductivity:
    def test_conductivityCF(self):
        with pytest.raises(TypeError):
            ossdbs.ConductivityCF(
                mri_image="Wrong Input",
                brain_region="Wrong Input",
                dielectric_model="Wrong Input",
                materials="Wrong Input",
            )


class TestVolumeConductorModel:
    def test_volume_conductor_floating(self):
        with pytest.raises(AttributeError):
            VolumeConductorFloating(
                geometry="Wrong Input",
                conductivity="Wrong Input",
                solver="Wrong Input",
                order="Wrong Input",
                meshing_parameters="Wrong Input",
            )
