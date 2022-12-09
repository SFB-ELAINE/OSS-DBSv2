
from src.brain_model import BrainModel
from src.input import Input
from src.strategy import QS_Strategy, EQS_Strategy
import ngsolve
import sys


def main(json_path: str) -> None:

    input = Input(json_path=json_path)
    brain_model = BrainModel(mri=input.mri())
    brain_model.add_electrodes(input.electrodes())

    result = QS_Strategy(boundary_values=input.boundary_values(),
                         brain_model=brain_model,
                         signal=input.stimulation_signal()
                         ).result()

    result.save(input.output_path())


if __name__ == '__main__':
    with ngsolve.TaskManager():
        main(sys.argv[1])
