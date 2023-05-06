from .abbott_stjude_model import AbbottStJudeParameters
from .abbott_stjude_model import AbbottStJudeActiveTipModel
from .abbott_stjude_model import AbbottStJudeDirectedModel


class AbbottStJudeActiveTip6146_6149(AbbottStJudeActiveTipModel):

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        parameters = AbbottStJudeParameters(tip_length=1.1,
                                            contact_length=1.5,
                                            contact_spacing=1.5,
                                            lead_diameter=1.3,
                                            total_length=100.0)
        super().__init__(parameters, rotation, direction, position)


class AbbottStJudeActiveTip6142_6145(AbbottStJudeActiveTipModel):

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        parameters = AbbottStJudeParameters(tip_length=2.6,
                                            contact_length=1.5,
                                            contact_spacing=0.5,
                                            lead_diameter=1.3,
                                            total_length=100.0)
        super().__init__(parameters, rotation, direction, position)


class AbbottStJudeDirected6172(AbbottStJudeDirectedModel):

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        parameters = AbbottStJudeParameters(tip_length=1.1,
                                            contact_length=1.5,
                                            contact_spacing=0.5,
                                            lead_diameter=1.3,
                                            total_length=100.0)
        super().__init__(parameters, rotation, direction, position)


class AbbottStJudeDirected6173(AbbottStJudeDirectedModel):

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        parameters = AbbottStJudeParameters(tip_length=1.1,
                                            contact_length=1.5,
                                            contact_spacing=1.5,
                                            lead_diameter=1.3,
                                            total_length=100.0)
        super().__init__(parameters, rotation, direction, position)
