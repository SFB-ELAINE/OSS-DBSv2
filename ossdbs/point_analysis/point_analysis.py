

from ossdbs.point_analysis.controller import Controller


def point_analysis(input: dict) -> None:

    controller = Controller(paramters=input)
    mesh = controller.mesh()
    solver = controller.solver()
    contacts = controller.contacts()
    conductivity = controller.conductivity()
    strategy = controller.spectrum_mode()
    volume_conductor = controller.volume_conductor()
    points = controller.points()
    volume_conductor = volume_conductor(conductivity=conductivity,
                                        mesh=mesh,
                                        contacts=contacts,
                                        solver=solver)
    result = strategy.compute(signal=controller.stimulation_signal(),
                              volume_conductor=volume_conductor,
                              points=points)

    result.save_by_categories("test_result.hdf5", controller.categories())
