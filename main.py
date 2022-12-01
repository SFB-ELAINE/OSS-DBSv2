
from src.dielectric_distribution import DielectricDistribution
from src.brain_model import BrainModel
from src.configuration import Configuration
from src.volume_conductor_model import VolumeConductor
from src.brainsubstance import Material
from src.fastfouriertransform import FFT
import json
import ngsolve
import sys
import os
import numpy as np


def main(json_path: str) -> None:

    with open(json_path, 'r') as json_file:
        input = json.load(json_file)

    configuration = Configuration(input=input)
    electrodes = configuration.electrodes()
    mri = configuration.magnetic_resonance_image()
    boundary_values = configuration.boundary_values()
    signal = configuration.stimulation_signal()

    brain_model = BrainModel(mri=mri, electrodes=electrodes)
    csf_position = mri.material_distribution(Material.CSF)
    mesh = brain_model.generate_mesh(order=2)
    mesh.mark_elements_by_position(position=csf_position)
    mesh.refine()

    dielectric = DielectricDistribution(mri=mri)
    waves = FFT(signal).sine_waves()
    with ngsolve.TaskManager():
        frequency = waves[0].frequency
        conductivity = dielectric.complex_conductivity(frequency=frequency)
        model = VolumeConductor(conductivity=conductivity, mesh=mesh)
        potential, error = model.evaluate_potential(boundary_values)
        data = potential.vec.data * np.real(waves[0].amplitude) / 2
        potential_sum = ngsolve.GridFunction(space=mesh.sobolev_space())
        potential_sum.vec.data += data

        for wave in waves[1:1]:
            frequency = wave.frequency
            conductivity = dielectric.complex_conductivity(frequency)
            model = VolumeConductor(conductivity=conductivity, mesh=mesh)
            potential, error = model.evaluate_potential(boundary_values)
            potential_sum.vec.data += potential.vec.data * np.real(wave.amplitude)

    conductivities = ngsolve.VoxelCoefficient(start=conductivity.start,
                                              end=conductivity.end,
                                              values=conductivity.data,
                                              linear=False)

    P = ngsolve.Integrate(ngsolve.grad(potential_sum) *
                          ngsolve.Conj(conductivities *
                                       ngsolve.grad(potential_sum)),
                          mesh.ngsolvemesh())

    if P:
        print('impedance: ', 1 / P)
    else:
        print('impedance: ', 'inf')

    # # directory = os.path.dirname(os.path.realpath(__file__))
    # # output_path = configuration.output_path()
    # # file_dir = os.path.join(directory, os.path.dirname(output_path))
    # # if not os.path.exists(file_dir):
    # #     os.mkdir(file_dir)

    # # file_name = os.path.join(file_dir, os.path.basename(output_path))

    # output = ngsolve.VTKOutput(ma=mesh.ngsolvemesh(),
    #                            coefs=[potential, -ngsolve.grad(potential)*1e3],
    #                            names=["potential", "field"],
    #                            filename=configuration.output_path(),
    #                            subdivision=0)
    # output.Do()


if __name__ == '__main__':
    main(sys.argv[1])
