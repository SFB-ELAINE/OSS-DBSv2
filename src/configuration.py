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
                     for i in range(1, len(electrode['Contact_values']) + 1)}
            names['Body'] = 'E{}B'.format(index)
            electrode = ElectrodeCreator.create(parameters=par)
            electrode.rename_boundaries(names)
            electrodes.append(electrode)
        return electrodes

    def boundary_values(self):
        boundary_values = {}
        for index, electrode in enumerate(self.__input['Electrodes'], 1):
            body_value = electrode['Body_value']
            body = {}
            if body_value is not None:
                body.update({'E{}B'.format(index): body_value})
            contact_values = electrode['Contact_values']
            contacts = {'E{}C{}'.format(index, i): value
                        for i, value in enumerate(contact_values, 1)
                        if value is not None}
            boundary_values.update(body | contacts)
        if self.__input['BrainSurface_value'] is not None:
            boundary_values.update({'Brain':
                                    self.__input['BrainSurface_value']})
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
