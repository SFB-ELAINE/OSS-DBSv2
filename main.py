
from src.brain_model import BrainModel
from src.volume_conductor_model import VolumeConductor
from src.brain_imaging import DefaultMagneticResonanceImage
from src.electrodes import ElectrodeCreator, ElectrodeParameters
from src.brainsubstance import Material
from src.dielectric_model import DielectricModel1
import ngsolve
import numpy as np


def main():

    boundaries = {"contact": 1.0, "wire": 0.0}
    conductivities = {"saline": 1278*1e-6/1e-2}

    electrode_parameters = ElectrodeParameters()
    electrode = ElectrodeCreator(parameters=electrode_parameters)

    mri = DefaultMagneticResonanceImage(file_path='')

    brain_model = BrainModel(mri=mri, electrode=electrode)
    mesh = brain_model.generate_mesh(order=2, boundaries=boundaries)

    mesh.mark_elements_by_material(mri=mri,
                                   material=Material.CEREBRO_SPINAL_FLUID)
    mesh.refine()


    dielectric_csf = DielectricModel1.create_model(material=Material.CEREBRO_SPINAL_FLUID)
    dielectric_gm = DielectricModel1.create_model(material=Material.GRAY_MATTER)
    dielectric_wm = DielectricModel1.create_model(material=Material.WHITE_MATTER)

    conductivity = np.zeros(mri.data_map.shape)
    conductivity[mri.data_map == Material.CEREBRO_SPINAL_FLUID] = dielectric_csf.conductivity(frequency=0)
    conductivity[mri.data_map == Material.GRAY_MATTER] = dielectric_gm.conductivity(frequency=0)
    conductivity[mri.data_map == Material.WHITE_MATTER] = dielectric_gm.conductivity(frequency=0)
    conductivity[mri.data_map == Material.UNKNOWN] = dielectric_wm.conductivity(frequency=0)

    start, end = mri.bounding_box()
    conductivities = ngsolve.VoxelCoefficient(start=start,
                                              end=end,
                                              values=conductivity,
                                              linear=False)

    # conduct = [conductivities[mat] for mat in mesh.materials()]
    # conductivities = ngsolve.CoefficientFunction(coef=conduct)

    model = VolumeConductor(conductivity=conductivities)

    with ngsolve.TaskManager():
        field, contact, P, potential = model.evaluate_potential(mesh)

    print('field: ', field)
    print('voltage_contact: ', contact)
    print('impedance: ', 1 / P)


if __name__ == '__main__':
    main()
