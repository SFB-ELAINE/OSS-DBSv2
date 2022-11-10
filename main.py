
from src.brain_model import BrainModel
from src.configuration import Configuration
from src.volume_conductor_model import VolumeConductor
from src.brainsubstance import Material
import ngsolve
import os


INPUT = {
    'Electrodes': [
        {
            'Name': 'MicroProbesCustomRodent',
            'Rotation': 0.0,
            'Direction': [0.0, 0.0, 1.0],
            'Translation': [5.0, 5.0, 5.0],
            'Contact_values': [1.0, ],
            'Body_value': 0.0,
        },
    ],
    'MagneticResonanceImage':
        {
            'Path': ''
        },
    'DiffusionTensorImage':
        {
            'Path': ''
        },
    'Output_directoy': 'result'
}


def main():

    configuration = Configuration(input=INPUT)
    electrode = configuration.electrodes()[0]
    mri = configuration.magnetic_resonance_image()
    dti = configuration.diffusion_tensor_image()

    with ngsolve.TaskManager():
        brain_model = BrainModel(mri=mri, electrode=electrode)
        print('MRI_DIMENSION', brain_model.bounding_box())
        mesh = brain_model.generate_mesh(order=2)

        csf_position = brain_model.material_distribution(Material.CSF)
        mesh.mark_elements_by_position(position=csf_position)
        mesh.refine()

        conductivity = brain_model.conductivity(frequency=0)
        conductivities = ngsolve.VoxelCoefficient(start=conductivity.start,
                                                  end=conductivity.end,
                                                  values=conductivity.data,
                                                  linear=False)
        diffusion = ngsolve.VoxelCoefficient(start=conductivity.start,
                                             end=conductivity.end,
                                             values=dti.diffusion(),
                                             linear=False)
      #  with diffusion the computation burden explodes 
      #  conductivities = diffusion * conductivities

        model = VolumeConductor(conductivity=conductivities)
        potential, error = model.evaluate_potential(mesh)

    P = ngsolve.Integrate(ngsolve.grad(potential) *
                          ngsolve.Conj(conductivities *
                                       ngsolve.grad(potential)),
                          mesh.ngsolvemesh())

    if P:
        print('impedance: ', 1 / P)
    else:
        print('impedance: ', 'inf')

    # directory = os.path.dirname(os.path.realpath(__file__))
    # output_path = configuration.output_path()
    # file_dir = os.path.join(directory, os.path.dirname(output_path))
    # if not os.path.exists(file_dir):
    #     os.mkdir(file_dir)

    # file_name = os.path.join(file_dir, os.path.basename(output_path))

    output = ngsolve.VTKOutput(ma=mesh.ngsolvemesh(),
                               coefs=[potential, -ngsolve.grad(potential)*1e3],
                               names=["potential", "field"],
                               filename=configuration.output_path(),
                               subdivision=0)
    output.Do()


if __name__ == '__main__':
    main()
