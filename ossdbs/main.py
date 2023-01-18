
import ngsolve
import sys

from ossdbs.volume_conductor_model import VolumeConductor
from ossdbs.input import Input


def main() -> None:
    with ngsolve.TaskManager():
        ossdbs_fem(sys.argv[1])


def ossdbs_fem(json_path: str) -> None:

    input = Input(json_path=json_path)
    mesh = input.mesh()
    boundaries = list(input.boundary_values().keys())
    mesh.refine_by_boundaries(boundaries)

    conductivity = input.conductivity()
    volume_conductor = VolumeConductor(conductivity=conductivity, mesh=mesh)
    strategy = input.spectrum_mode()

    output = strategy.result(boundary_values=input.boundary_values(),
                             signal=input.stimulation_signal(),
                             volume_conductor=volume_conductor
                             )
    output.save(input.output_path())
    output.save_mesh()


if __name__ == '__main__':
    main()
