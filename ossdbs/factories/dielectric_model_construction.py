
from ossdbs.dielectric_model import DielectricModel
from ossdbs.dielectric_model import ColeCole4Model
from ossdbs.dielectric_model import ColeCole4ModelCustom
from ossdbs.dielectric_model import ConstantModel
from ossdbs.dielectric_model import ConstantModelCustom


class DielectricModelFactory:

    def create(cls, dielectricum_model: str) -> DielectricModel:
        return {'ColeCole4': ColeCole4Model(),
                'ColeCole4Custom': ColeCole4ModelCustom(),
                'Constant': ConstantModel(),
                'ConstantCustom': ConstantModelCustom(),
                }[dielectricum_model]
