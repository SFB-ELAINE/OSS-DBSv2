# Boston Scientific (Marlborough, Massachusetts, USA) vercise

from .boston_scientific_model import BostonScientificVerciseParameters
from .boston_scientific_model import BostonScientificVerciseModel
from .boston_scientific_model import BostonScientificVerciseDirectedModel


class BostonScientificVercise(BostonScientificVerciseModel):

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        parameters = BostonScientificVerciseParameters(tip_length=1.1,
                                                       contact_length=1.5,
                                                       contact_spacing=0.5,
                                                       lead_diameter=1.3,
                                                       total_length=100.0)
        super().__init__(parameters, rotation, direction, position)


class BostonScientificVerciseDirected(BostonScientificVerciseDirectedModel):

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        parameters = BostonScientificVerciseParameters(tip_length=1.5,
                                                       contact_length=1.5,
                                                       contact_spacing=0.5,
                                                       lead_diameter=1.3,
                                                       total_length=100.0)
        super().__init__(parameters, rotation, direction, position)
