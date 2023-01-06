#! /usr/bin/python3


import ngsolve
import sys

from ossdbs.volume_conductor_model import VolumeConductor
from ossdbs.input import Input
from ossdbs.brain_geometry import BrainGeometry


def main(json_path: str) -> None:

    input = Input(json_path=json_path)
    brain_geometry = BrainGeometry(region=input.region_of_interest())
    brain_geometry.set_electrodes(input.electrodes())

    mesh = brain_geometry.generate_mesh(input.mesh_order())
    boundaries = list(input.boundary_values().keys())
    mesh.refine_by_boundaries(boundaries)
    mesh.set_complex(input.complex_mode())

    conductivity = input.conductivity()
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
