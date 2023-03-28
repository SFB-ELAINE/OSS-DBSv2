
import json
from typing import List
from ossdbs.electrodes.contacts import Contact, Contacts
from ossdbs.electrodes.electrode_models import ElectrodeModel
from ossdbs.electrodes.electrodes import Electrodes
from ossdbs.factories.custom_electrode_factory import CustomElectrodeFactory
from ossdbs.factories.electrode_factory import ElectrodeFactory


class ElectrodesFactory:

    @classmethod
    def create(cls, electrode_paramters: dict) -> Electrodes:

        electrodes = []
        contacts = []

        for index, input_par in enumerate(electrode_paramters):
            electrode = cls.__create_electrode(index, input_par)
            electrode.set_contact_names(cls.__boundaries(index))
            electrodes.append(electrode)
            contacts.extend([cls.__create_contact(index, contact_par)
                             for contact_par in input_par['Contacts']])

        contacts = Contacts([contact for contact in contacts
                             if contact.active or contact.floating])

        return Electrodes(electrodes=electrodes, contacts=contacts)

    @staticmethod
    def __create_contact(index: int, contact_par: dict) -> Contact:
        name = 'E{}C{}'.format(index, contact_par['Contact_ID'] - 1)
        active = contact_par['Active']
        floating = contact_par['Floating'] and not active
        current = contact_par['Current[A]']
        voltage = contact_par['Voltage[V]']
        real = contact_par['SurfaceImpedance[Ωm]']['real']
        imag = contact_par['SurfaceImpedance[Ωm]']['imag']
        surface_impedance = real + 1j * imag
        return Contact(name=name,
                       active=active,
                       floating=floating,
                       current=current,
                       voltage=voltage,
                       surface_impedance=surface_impedance)

    def __create_electrode(cls, electrode_parameters: dict) -> ElectrodeModel:
        dir = electrode_parameters['Direction']
        pos = electrode_parameters['TipPosition']
        direction = (dir['x[mm]'], dir['y[mm]'], dir['z[mm]'])
        position = (pos['x[mm]'], pos['y[mm]'], pos['z[mm]'])
        rotation = electrode_parameters['Rotation[Degrees]']

        if 'Custom' in electrode_parameters['Name']:
            model_parameters = cls.read_json(electrode_parameters)
            return CustomElectrodeFactory(name=electrode_parameters['Name'],
                                          model_parameters=model_parameters,
                                          direction=direction,
                                          position=position,
                                          rotation=rotation)

        return ElectrodeFactory.create(name=electrode_parameters['Name'],
                                       direction=direction,
                                       position=position,
                                       rotation=rotation)

    @staticmethod
    def read_json(electrode_parameters):
        path = electrode_parameters['PathToCustomParameters']
        with open(path, 'r') as json_file:
            model_parameters = json.load(json_file)
        return model_parameters

    @staticmethod
    def __boundaries(index: int) -> dict:
        return {'Contact_1': "E{}C0".format(index),
                'Contact_2': "E{}C1".format(index),
                'Contact_3': "E{}C2".format(index),
                'Contact_4': "E{}C3".format(index),
                'Contact_5': "E{}C4".format(index),
                'Contact_6': "E{}C5".format(index),
                'Contact_7': "E{}C6".format(index),
                'Contact_8': "E{}C7".format(index),
                'Body': 'E{}B'.format(index)}
