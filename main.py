from src.geometry import SimpleGeometry
from src.mesh import Mesh
from src.volume_conductor_model import VolumeConductor
from src.brain_imaging.magnetic_resonance_imaging import MagneticResonanceImage
from src.brain_imaging.diffusion_tensor_imaging import DiffusionTensorImage
import ngsolve


def main():

    # mri = MagneticResonanceImage(file_path='')
    # dti = DiffusionTensorImage(file_path='')
    # electrode = factory(electrode_information)

    geometry = SimpleGeometry()
    boundaries = {"contact": 1.0, "wire": 0.0}
    conductivities = {"saline": 1278*1e-6/1e-2}

    mesh = Mesh(geometry=geometry, order=2, boundaries=boundaries)
    model = VolumeConductor(mesh=mesh, conductivity=conductivities)

    with ngsolve.TaskManager():
        field, contact, P, potential = model.evaluate_potential()

    print('field: ', field)
    print('voltage_contact: ', contact)
    print('impedance: ', 1 / P)


if __name__ == '__main__':
    main()
