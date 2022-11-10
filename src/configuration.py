from src.brain_imaging.diffusion_tensor_imaging \
    import DefaultDiffusionTensorImage, DiffusionTensorImage
from src.brain_imaging.magnetic_resonance_imaging \
    import DefaultMagneticResonanceImage, MagneticResonanceImage
from src.electrodes.electrodes import ElectrodeCreator, ElectrodeParameters


class Configuration:

    def __init__(self, input: dict) -> None:
        self.__input = input

    def electrodes(self) -> list:
        electrodes = []
        for index, electrode in enumerate(self.__input['Electrodes'], 1):
            par = ElectrodeParameters()
            par.name = electrode['Name']
            par.rotation = electrode['Rotation']
            par.direction = electrode['Direction']
            par.translation = tuple(electrode['Translation'])

            values = electrode['Contact_values']
            par.contact_values = [("E{}C{}".format(index, i), value)
                                  for i, value in enumerate(values, 1)]
            electrodes.append(ElectrodeCreator.create(parameters=par))
        return electrodes

    def magnetic_resonance_image(self):
        path = self.__input['MagneticResonanceImage']['Path']
        if not path:
            return DefaultMagneticResonanceImage()
        return MagneticResonanceImage(file_path=path)

    def diffusion_tensor_image(self):
        path = self.__input['MagneticResonanceImage']['Path']
        if not path:
            return DefaultDiffusionTensorImage()
        return DiffusionTensorImage(file_path=path)

    def output_path(self):
        return self.__input['Output_directoy']
