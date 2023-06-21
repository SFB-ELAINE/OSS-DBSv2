from dataclasses import dataclass


@dataclass
class EncapsulationLayer:
    """Class to store information about encapsulation layer
    """

    name: str
    material: str = None
    dielectric_model: str = None
    dielectric_parameters: dict = None
    max_h: float = 1e6
