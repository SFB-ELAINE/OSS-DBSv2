
from src.brain_model import BrainModel
from src.input import Input
from src.strategy import QS_Strategy, EQS_Strategy
from src.brainsubstance import Material
import ngsolve
import sys


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

    result = strategy.result()

    conductivity = brain_model.complex_conductivity(frequency=0)
    conductivities = ngsolve.VoxelCoefficient(start=conductivity.start,
                                              end=conductivity.end,
                                              values=conductivity.data,
                                              linear=False)
    P = ngsolve.Integrate(ngsolve.grad(result.potential) *
                          ngsolve.Conj(conductivities *
                                       ngsolve.grad(result.potential)),
                          mesh.ngsolvemesh())

    print('impedance:', 'inf' if not P else 1 / P)
    result.save(input.output_path())


if __name__ == '__main__':
    with ngsolve.TaskManager():
        main(sys.argv[1])
