
from src.brain_model import BrainModel
from src.volume_conductor_model import VolumeConductor
from src.brain_imaging import DefaultMagneticResonanceImage
from src.electrodes import ElectrodeCreator, ElectrodeParameters
from src.brainsubstance import Material
import ngsolve


INPUT = {
    'Electrodes': [
        {
            'Name': 'Rodden',
            'Rotation': 0.0,
            'Direction': [0.0, 0.0, 0.0],
            'Translation': [5.0, 5.0, 5.0],
            'Contact_values': [1.0, ]
        }
    ],
    'MagneticResonanceImage':
        {
            'Path': 'path'
        },
    'DiffusionTensorImage':
        {
            'Path': 'path'
        }
}


class Configuration:

    def __init__(self, input: dict) -> None:
        self.__input = input

    def electrodes(self) -> list:
        electrodes = []
        for index, electrode in enumerate(self.__input['Electrodes'], 1):
            par = ElectrodeParameters()
            par.name = electrode['Name']
            par.rotation = electrode['Rotation']
            par.direction = electrode['Direction']
            par.translation = tuple(electrode['Translation'])

            contact_values = electrode['Contact_values']
            par.contact_values = [("E{}C{}".format(index, i), v)
                                  for i, v in enumerate(contact_values, 1)]
            electrodes.append(ElectrodeCreator.create(parameters=par))
        return electrodes


def main():
    # Boundary names & values in geometry and Electrode Parameters
    boundaries = {"Contact": 1.0, "Body": 0.0}

    electrode = Configuration(input=INPUT).electrodes()[0]
    mri = DefaultMagneticResonanceImage(file_path='')

    with ngsolve.TaskManager():
        brain_model = BrainModel(mri=mri, electrode=electrode)
        mesh = brain_model.generate_mesh(order=2, boundary_values=boundaries)

        position = brain_model.material_distribution(Material.CSF)
        mesh.mark_elements_by_position(position=position)
        mesh.refine()

        start, end = brain_model.bounding_box()
        conductivity = brain_model.conductivity(frequency=0)
        conductivities = ngsolve.VoxelCoefficient(start=tuple(start),
                                                  end=tuple(end),
                                                  values=conductivity,
                                                  linear=False)

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


if __name__ == '__main__':
    main()
