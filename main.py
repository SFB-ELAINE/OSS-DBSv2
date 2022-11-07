
from src.brain_model import BrainModel
from src.volume_conductor_model import VolumeConductor
from src.brain_imaging import DefaultMagneticResonanceImage
from src.electrodes import ElectrodeCreator, ElectrodeParameters
from src.brainsubstance import Material
import ngsolve


def main():

    boundaries = {"contact": 1.0, "wire": 0.0}
    conductivities = {"saline": 1278*1e-6/1e-2}

    electrode_parameters = ElectrodeParameters()
    electrode = ElectrodeCreator(parameters=electrode_parameters)

    mri = DefaultMagneticResonanceImage(file_path='')

    brain_model = BrainModel(mri=mri, electrode=electrode)
    mesh = brain_model.generate_mesh(order=2, boundaries=boundaries)

    position = brain_model.material_distribution(Material.CEREBRO_SPINAL_FLUID)
    start, end = brain_model.bounding_box()
    mesh.mark_elements_by_position(position=position, start=start, end=end)
    mesh.refine()

    conductivity = brain_model.conductivity(frequency=0)
    conductivities = ngsolve.VoxelCoefficient(start=start,
                                              end=end,
                                              values=conductivity,
                                              linear=False)

    # conduct = [conductivities[mat] for mat in mesh.materials()]
    # conductivities = ngsolve.CoefficientFunction(coef=conduct)

    model = VolumeConductor(conductivity=conductivities)

    with ngsolve.TaskManager():
        potential, error = model.evaluate_potential(mesh)

    P = ngsolve.Integrate(ngsolve.grad(potential) *
                          ngsolve.Conj(conductivities *
                                       ngsolve.grad(potential)),
                          mesh.ngsolvemesh())

    print('impedance: ', 1 / P)


if __name__ == '__main__':
    main()
