# Medtronic 3387

from .medtronic_model import MedtronicParameters, MedtronicModel


class Medtronic3387(MedtronicModel):

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        parameters = MedtronicParameters(tip_length=1.5,
                                         contact_length=1.5,
                                         contact_spacing=0.5,
                                         lead_diameter=1.27,
                                         total_length=100.0)
        super().__init__(parameters, rotation, direction, position)


class Medtronic3389(MedtronicModel):

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        parameters = MedtronicParameters(tip_length=1.5,
                                         contact_length=1.5,
                                         contact_spacing=1.5,
                                         lead_diameter=1.27,
                                         total_length=100.0)
        super().__init__(parameters, rotation, direction, position)


class Medtronic3391(MedtronicModel):

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        parameters = MedtronicParameters(tip_length=1.5,
                                         contact_length=3.0,
                                         contact_spacing=3.0,
                                         lead_diameter=1.27,
                                         total_length=100.0)
        super().__init__(parameters, rotation, direction, position)
