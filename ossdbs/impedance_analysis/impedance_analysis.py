from ossdbs.impedance_analysis.controller import Controller


def impedance_analysis(input: dict) -> None:

    controller = Controller(parameters=input)
    mesh = controller.mesh()
    solver = controller.solver()
    contacts = controller.contacts()
    conductivity = controller.conductivity()
    strategy = controller.spectrum_mode()
    volume_conductor = controller.volume_conductor()
    volume_conductor = volume_conductor(conductivity=conductivity,
                                        mesh=mesh,
                                        contacts=contacts,
                                        solver=solver)
    impedances = strategy.compute(signal=controller.stimulation_signal(),
                                  volume_conductor=volume_conductor)
    impedances.save('impedances.csv')
