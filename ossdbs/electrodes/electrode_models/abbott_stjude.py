from .abbott_stjude_model import AbbottStjudeParameters
from .abbott_stjude_model import AbbottStjudeActiveTipModel
from .abbott_stjude_model import AbbottStjudeDirectedModel


class AbbottStjudeActiveTip6146_6149(AbbottStjudeActiveTipModel):

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        parameters = AbbottStjudeParameters(tip_length=1.1,
                                           contact_length=1.5,
                                           contact_spacing=1.5,
                                           lead_diameter=1.3,
                                           total_length=100.0)
        super().__init__(parameters, rotation, direction, position)


class AbbottStjudeActiveTip6142_6145(AbbottStjudeActiveTipModel):

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        parameters = AbbottStjudeParameters(tip_length=2.6,
                                           contact_length=1.5,
                                           contact_spacing=0.5,
                                           lead_diameter=1.3,
                                           total_length=100.0)
        super().__init__(parameters, rotation, direction, position)


class AbbottStjudeDirected6172(AbbottStjudeDirectedModel):

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        parameters = AbbottStjudeParameters(tip_length=1.1,
                                           contact_length=1.5,
                                           contact_spacing=0.5,
                                           lead_diameter=1.3,
                                           total_length=100.0)
        super().__init__(parameters, rotation, direction, position)


class AbbottStjudeDirected6173(AbbottStjudeDirectedModel):

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        parameters = AbbottStjudeParameters(tip_length=1.1,
                                           contact_length=1.5,
                                           contact_spacing=1.5,
                                           lead_diameter=1.3,
                                           total_length=100.0)
        super().__init__(parameters, rotation, direction, position)
