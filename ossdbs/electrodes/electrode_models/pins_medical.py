# PINS Medical L301

from .pins_medical_model import PINSMedicalParameters, PINSMedicalModel


class PINSMedicalL301(PINSMedicalModel):

    def __init__(self,
                rotation: float = 0.0,
                direction: tuple = (0, 0, 1),
                position: tuple = (0, 0, 0)) -> None:
        parameters = PINSMedicalParameters(tip_length=1.1,
                                          contact_length=1.5,
                                          contact_spacing=0.5,
                                          lead_diameter=1.3,
                                          total_length=100.0)
        super().__init__(parameters, rotation, direction, position)


class PINSMedicalL302(PINSMedicalModel):

    def __init__(self,
                rotation: float = 0.0,
                direction: tuple = (0, 0, 1),
                position: tuple = (0, 0, 0)) -> None:
        parameters = PINSMedicalParameters(tip_length=1.1,
                                          contact_length=1.5,
                                          contact_spacing=1.5,
                                          lead_diameter=1.3,
                                          total_length=100.0)
        super().__init__(parameters, rotation, direction, position)


class PINSMedicalL303(PINSMedicalModel):

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        parameters = PINSMedicalParameters(tip_length=1.1,
                                          contact_length=3.0,
                                          contact_spacing=3.0,
                                          lead_diameter=1.3,
                                          total_length=100.0)
        super().__init__(parameters, rotation, direction, position)
