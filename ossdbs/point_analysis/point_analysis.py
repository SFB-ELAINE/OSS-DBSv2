
from dataclasses import dataclass
from ossdbs.point_analysis.controller import Controller
import numpy as np
import h5py


@dataclass
class Group:
    name: str
    potential: np.ndarray
    current_density: np.ndarray
    points: np.ndarray


def point_analysis(input: dict) -> None:

    controller = Controller(paramters=input)
    mesh = controller.mesh()
    solver = controller.solver()
    contacts = controller.contacts()
    conductivity = controller.conductivity()
    strategy = controller.spectrum_mode()
    volume_conductor = controller.volume_conductor()
    points = controller.coordinates()
    included = mesh.is_included(controller.coordinates())
    volume_conductor = volume_conductor(conductivity=conductivity,
                                        mesh=mesh,
                                        contacts=contacts,
                                        solver=solver)
    result = strategy.compute(signal=controller.stimulation_signal(),
                              volume_conductor=volume_conductor,
                              points=points[included])

    groups = new_func(controller, included, result)

    with h5py.File("test_result.hdf5", "w") as file:
        file.create_dataset('TimeSteps', data=result.time_steps)
        for group in groups:
            h5_group = file.create_group(group.name)
            h5_group.create_dataset('Points', data=group.points)
            h5_group.create_dataset('Potential', data=group.potential)
            h5_group.create_dataset('CurrentDensity',
                                    data=group.current_density)


def new_func(controller, included, result):
    potential_t = np.zeros(result.potential.shape)
    potential_t[included] = result.potential
    current_density_t = np.zeros(result.current_density.shape)
    current_density_t[included] = result.current_density
    points = controller.coordinates()

    groups = []
    start = 0
    for group in controller.categories():
        name = group[0]
        end = group[1][0] + start
        groups.append(Group(name=name,
                            potential=potential_t[start:end],
                            current_density=current_density_t[start:end],
                            points=points[start:end]))
                 
    return groups
