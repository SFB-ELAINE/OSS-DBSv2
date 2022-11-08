
from src.brain_model import BrainModel
from src.volume_conductor_model import VolumeConductor
from src.brain_imaging import DefaultMagneticResonanceImage
from src.electrodes import ElectrodeCreator, ElectrodeParameters
from src.brainsubstance import Material
import ngsolve


def main():
    # Boundary names & values in geometry and Electrode Parameters
    boundaries = {"Contact": 1.0, "Body": 0.0}

    electrode_parameters = ElectrodeParameters(translation=(5, 5, 5))
    electrode = ElectrodeCreator.create(parameters=electrode_parameters)

    mri = DefaultMagneticResonanceImage(file_path='')

    with ngsolve.TaskManager():
        brain_model = BrainModel(mri=mri, electrode=electrode)
        mesh = brain_model.generate_mesh(order=2, boundary_values=boundaries)

        position = brain_model.material_distribution(Material.CSF)
        start, end = brain_model.bounding_box()

        mesh.mark_elements_by_position(position=position, start=start, end=end)
        mesh.refine()

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
