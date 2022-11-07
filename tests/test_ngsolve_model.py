from src.geometry import SimpleGeometry
from src.mesh import Mesh
from src.volume_conductor_model import VolumeConductor
import ngsolve


def test_ngsolve_model():
    geometry = SimpleGeometry()
    boundaries = {"contact": 1.0, "wire": 0.0}
    conductivities = {"saline": 1278*1e-6/1e-2}

    with ngsolve.TaskManager():
        mesh = Mesh(geometry=geometry,
                    order=2,
                    boundaries=boundaries)

        conduct = [conductivities[mat] for mat in mesh.materials()]
        sigma = ngsolve.CoefficientFunction(coef=conduct)
        model = VolumeConductor(conductivity=sigma)
        potential, error = model.evaluate_potential(mesh=mesh)

        P = ngsolve.Integrate(ngsolve.grad(potential) *
                              ngsolve.Conj(sigma * ngsolve.grad(potential)),
                              mesh.ngsolvemesh())

    assert round(1/P, 3) == 5676.935
