#! /usr/bin/python3


import ngsolve
import sys

from ossdbs.volume_conductor_model import VolumeConductor

from ossdbs.conductivity import Conductivity
from ossdbs.input import Input
from ossdbs.brain_geometry import BrainGeometry


def main(json_path: str) -> None:

    input = Input(json_path=json_path)
    brain_model = BrainGeometry(region=input.region_of_interest())
    brain_model.set_electrodes(input.electrodes())
    mesh = brain_model.generate_mesh(input.mesh_order())
    boundaries = list(input.boundary_values().keys())
    mesh.refine_by_boundaries(boundaries)
    conductivity = Conductivity(input.mri())
    conductivity.set_complex(input.complex_mode())
    mesh.set_complex(input.complex_mode())
    volume_conductor = VolumeConductor(conductivity=conductivity, mesh=mesh)
    strategy = input.spectrum_mode()
    output = strategy.result(boundary_values=input.boundary_values(),
                             signal=input.stimulation_signal(),
                             volume_conductor=volume_conductor
                             )
    output.save(input.output_path())


if __name__ == '__main__':
    with ngsolve.TaskManager():
        main(sys.argv[1])
