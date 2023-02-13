
import ngsolve
import sys
from ossdbs.controller import Controller
from ossdbs.controller_impedance import ControllerImpedance


def main() -> None:

    if 'impedance' in sys.argv:
        with ngsolve.TaskManager():
            ossdbs_impedance(sys.argv[1])
        return

    with ngsolve.TaskManager():
       ossdbs_fem(sys.argv[1])


def ossdbs_fem(json_path: str) -> None:

    controller = Controller(json_path=json_path)
    mesh = controller.mesh()
    solver = controller.solver()
    contacts = controller.contacts()
    conductivity = controller.conductivity()
    strategy = controller.spectrum_mode()
    volume_conductor = controller.volume_conductor()
    mesh.refine_by_boundaries(contacts.active())
    points = mesh.included_points(controller.coordinates())
    included = mesh.is_included(controller.coordinates())
    volume_conductor = volume_conductor(conductivity=conductivity,
                                        mesh=mesh,
                                        contacts=contacts,
                                        solver=solver)
    result = strategy.compute(signal=controller.stimulation_signal(),
                              volume_conductor=volume_conductor,
                              points=points)
    result.save('result.h5')


    print(controller.categories())

    import numpy as np
    time_result = result.time_result()
    potential_t = np.zeros((len(included), 10000))
    current_density_t = np.zeros((len(included), 10000, 3))
    potential_t[included] = time_result.potential
    current_density_t[included] = time_result.current_density



def ossdbs_impedance(json_path: str) -> None:

    controller = ControllerImpedance(json_path=json_path)
    mesh = controller.mesh()
    solver = controller.solver()
    contacts = controller.contacts()
    conductivity = controller.conductivity()
    strategy = controller.spectrum_mode()
    volume_conductor = controller.volume_conductor()
    mesh.refine_by_boundaries(contacts.active())
    volume_conductor = volume_conductor(conductivity=conductivity,
                                        mesh=mesh,
                                        contacts=contacts,
                                        solver=solver)
    result = strategy.compute(signal=controller.stimulation_signal(),
                              volume_conductor=volume_conductor)
    result.save('result.csv')



if __name__ == '__main__':
    main()
