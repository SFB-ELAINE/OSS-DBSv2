
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
    solver = input.solver()
    contacts = input.contacts()
    conductivity = input.conductivity()

    mesh.refine_by_boundaries(contacts.active_contacts())
    volume_conductor = VolumeConductor(conductivity=conductivity,
                                       mesh=mesh,
                                       boundaries=contacts,
                                       solver=solver)
    strategy = input.spectrum_mode()
    output = strategy.result(signal=input.stimulation_signal(),
                             volume_conductor=volume_conductor
                             )
    output.save(input.output_path())
    output.save_mesh()


if __name__ == '__main__':
    main()
