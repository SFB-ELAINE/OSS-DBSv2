from src.geometry import Geometry
from src.mesh import Mesh
from src.volume_conductor_model import VolumeConductorQS
import ngsolve
import numpy as np
import netgen.occ as occ

class Electrode:

    def __init__(self) -> None:
        radius = 225e-6 / 2 * 1.1   # [m]
        direction = occ.Z
        height = 0.02
        shaft = occ.Cylinder(p=(-0.005, 0, radius),
                             d=direction,
                             r=radius,
                             h=height)
        shaft.bc('shaft')
        shaft.maxh = radius * 2 * 50
        contact = occ.Sphere(c=(-0.005, 0, 0), r=radius) - shaft
        contact.bc('contact')
        contact.maxh = 50 * radius / 20
        self.solid = occ.Glue([contact, shaft])


class Wire:

    def __init__(self) -> None:
        r_wire = 250e-6 / 2
        height = 0.02
        wire = occ.Cylinder((0.005, 0, 0), occ.Z, r=r_wire, h=height)
        wire.bc('wire')
        wire.maxh = 50 * r_wire
        self.solid = wire


class Brain:

    def __init__(self) -> None:
        sphere = occ.Sphere((0, 0, 0), 0.02)
        sphere.bc('brain')
        sphere.mat('saline')
        self.solid = sphere


class SimpleGeometry(Geometry):

    def __init__(self) -> None:
        electrode = Electrode()
        wire = Wire()
        brain = Brain()
        model = brain.solid - electrode.solid - wire.solid
        self.__geometry = occ.OCCGeometry(model)

    def generate_mesh(self):
        return self.__geometry.GenerateMesh()


def test_ngsolve_model():
    geometry = SimpleGeometry()
    boundaries = {"contact": 1.0, "wire": 0.0}
    conductivities = {"saline": 1278*1e-6/1e-2}

    with ngsolve.TaskManager():
        mesh = Mesh(geometry=geometry, order=2)
        conduct = [conductivities[mat] for mat in mesh.materials()]
        sigma = ngsolve.CoefficientFunction(coef=conduct)
    #     model = VolumeConductor(conductivity2=sigma)
    #     potential, error = model.evaluate_potential(mesh=mesh,
    #                                                 boundaries=boundaries)

    #     P = ngsolve.Integrate(ngsolve.grad(potential) *
    #                           ngsolve.Conj(sigma * ngsolve.grad(potential)),
    #                           mesh.ngsolvemesh())

    # np.testing.assert_allclose(1/P, 2227, atol=1)
