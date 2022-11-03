
from src.brain_model import BrainModel
from src.volume_conductor_model import VolumeConductor
from src.brain_imaging.magnetic_resonance_imaging \
    import DefaultNMagneticResonanceImage
from src.electrode_creator import ElectrodeCreator, ElectrodeParameters
import ngsolve


def main():

    boundaries = {"contact": 1.0, "wire": 0.0}
    conductivities = {"saline": 1278*1e-6/1e-2}

    electrode_parameters = ElectrodeParameters()
    electrode = ElectrodeCreator(parameters=electrode_parameters)

    mri = DefaultNMagneticResonanceImage(file_path='')

    brain_model = BrainModel(mri=mri, electrode=electrode)
    mesh = brain_model.generate_mesh(order=2, boundaries=boundaries)

    model = VolumeConductor(mesh=mesh, conductivity=conductivities)

    with ngsolve.TaskManager():
        field, contact, P, potential = model.evaluate_potential()

    print('field: ', field)
    print('voltage_contact: ', contact)
    print('impedance: ', 1 / P)


if __name__ == '__main__':
    main()
