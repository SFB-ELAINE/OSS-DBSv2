#! /usr/bin/python3


import ngsolve
import sys
try:
    import ossdbs
except ImportError:
    print("ossdbs has not been properly installed")

from ossdbs.conductivity import Conductivity
from ossdbs.input import Input
from ossdbs.strategy import AllFrequenciesStrategy
from ossdbs.brain_geometry import BrainGeometry


def main(json_path: str) -> None:

    input = Input(json_path=json_path)
    brain_model = BrainGeometry(region=input.region_of_interest())
    brain_model.set_electrodes(input.electrodes())
    mesh = brain_model.generate_mesh(input.mesh_order())
    boundaries = list(input.boundary_values().keys())
    mesh.refine_by_boundaries(boundaries)
    # mesh.refine_by_mri(input.mri())
    conductivity = Conductivity(input.mri())

    vc_type = input.volume_conductor_type()
    volume_conductor = vc_type(conductivity=conductivity, mesh=mesh)
    output = AllFrequenciesStrategy(boundary_values=input.boundary_values(),
                                    signal=input.stimulation_signal(),
                                    volume_conductor=volume_conductor
                                    ).result()
    output.save(input.output_path())


if __name__ == '__main__':
    with ngsolve.TaskManager():
        main(sys.argv[1])
