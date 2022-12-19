
from src.brain_model import BrainModel
from src.conductivity import Conductivity
from src.input import Input
from src.mesh_refinement import MeshRefinement
from src.strategy import AllFrequenciesStrategy, StrategyOctavevands
import ngsolve
import sys


def main(json_path: str) -> None:

    input = Input(json_path=json_path)
    brain_model = BrainModel(region=input.region_of_interest())
    brain_model.add_electrodes(input.electrodes())
    mesh = brain_model.generate_mesh(input.mesh_order())
    boundaries = list(input.boundary_values().keys())
    MeshRefinement(mesh).refine_by_boundaries(boundaries)
    # MeshRefinement(mesh).refine_by_mri(mri)
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
