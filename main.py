
from src.brain_model import BrainModel
from src.input import Input
from src.strategy import QS_Strategy
from src.brainsubstance import Material
import ngsolve
import sys
import os


def main(json_path: str) -> None:

    input = Input(json_path=json_path)

    mri = input.mri()
    csf_position = mri.material_distribution(Material.CSF)
    brain_model = BrainModel(mri=input.mri())
    [brain_model.add_electrode(electrode) for electrode in input.electrodes()]

    mesh = brain_model.generate_mesh(order=2)
    mesh.mark_elements_by_position(position=csf_position)
    mesh.refine()

    strategy = QS_Strategy(mesh=mesh,
                           boundary_values=input.boundary_values(),
                           brain_model=brain_model,
                           signal=input.stimulation_signal())

    potential = strategy.potential()

    conductivity = brain_model.complex_conductivity(frequency=0)
    conductivities = ngsolve.VoxelCoefficient(start=conductivity.start,
                                              end=conductivity.end,
                                              values=conductivity.data,
                                              linear=False)
    P = ngsolve.Integrate(ngsolve.grad(potential) *
                          ngsolve.Conj(conductivities *
                                       ngsolve.grad(potential)),
                          mesh.ngsolvemesh())
    
    print('impedance:', 'inf' if not P else 1 / P)


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
    with ngsolve.TaskManager():
        main(sys.argv[1])
