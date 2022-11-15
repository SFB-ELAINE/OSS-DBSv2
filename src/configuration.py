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
            names = {'Contact_{}'.format(i): "E{}C{}".format(index, i)
                     for i, _ in enumerate(electrode['Contacts']['Active'], 1)}
            names['Body'] = 'E{}B'.format(index)
            electrode = ElectrodeCreator.create(parameters=par)
            electrode.rename_boundaries(names)
            electrodes.append(electrode)
        return electrodes

    def boundary_values(self):
        boundary_values = {}
        for index, electrode in enumerate(self.__input['Electrodes'], 1):
            if electrode['Body']['Active']:
                value = electrode['Body']['Value']
                boundary_values.update({'E{}B'.format(index): value})

            contacts = {'E{}C{}'.format(index, i): value
                        for i, value
                        in enumerate(electrode['Contacts']['Value'], 1)
                        if electrode['Contacts']['Active'][i-1]}
            boundary_values.update(contacts)

        if self.__input['BrainSurface']['Active']:
            value = self.__input['BrainSurface_value']['Value']
            boundary_values.update({'Brain': value})

        return boundary_values

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
