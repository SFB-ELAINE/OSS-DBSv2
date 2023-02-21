
from typing import List
from ossdbs.electrode_contacts import ElectrodeContact
from ossdbs.electrode_creation import ElectrodeFactory


class Electrodes:

    def __init__(self, electrode_paramters: dict) -> None:

        self.__electrodes = []
        self.__contacts = []

        for index, input_par in enumerate(electrode_paramters):
            electrode = self.__create_electrode(input_par)
            boundary_names = self.__boundaries(index)
            electrode.rename_boundaries(boundary_names)
            self.__electrodes.append(electrode)
            self.__contacts.extend(self.__create_contacts(index, electrode))

    def __create_contacts(self, index, electrode):
        for contact_par in electrode['Contacts']:
            name = 'E{}C{}'.format(index, contact_par['Contact_ID'] - 1)
            active = contact_par['Active']
            floating = contact_par['Floating'] and not active
            current = contact_par['Current[A]']
            voltage = contact_par['Voltage[V]']
            real = contact_par['SurfaceImpedance[Ωm]']['real']
            imag = contact_par['SurfaceImpedance[Ωm]']['imag']
            surface_impedance = real + 1j * imag
            return ElectrodeContact(name=name,
                                    active=active,
                                    floating=floating,
                                    current=current,
                                    voltage=voltage,
                                    surface_impedance=surface_impedance)

    def __create_electrode(self, input_par):
        dir = input_par['Direction']
        pos = input_par['TipPosition']
        direction = (dir['x[mm]'], dir['y[mm]'], dir['z[mm]'])
        position = (pos['x[mm]'], pos['y[mm]'], pos['z[mm]'])
        rotation = input_par['Rotation[Degrees]']
        return ElectrodeFactory.create(name=input_par['Name'],
                                       direction=direction,
                                       position=position,
                                       rotation=rotation)

    def __boundaries(self, index: int) -> dict:
        return {'Contact_1': "E{}C0".format(index),
                'Contact_2': "E{}C1".format(index),
                'Contact_3': "E{}C2".format(index),
                'Contact_4': "E{}C3".format(index),
                'Contact_5': "E{}C4".format(index),
                'Contact_6': "E{}C5".format(index),
                'Contact_7': "E{}C6".format(index),
                'Contact_8': "E{}C7".format(index),
                'Body': 'E{}B'.format(index)}

    def current_values(self) -> dict:
        return {contact.name: contact.current
                for contact in self.__contacts if contact.active}

    def voltage_values(self) -> dict:
        return {contact.name: contact.voltage
                for contact in self.__contacts if contact.active}

    def floating_impedance_values(self) -> dict:
        return {contact.name: contact.surface_impedance
                for contact in self.__contacts
                if contact.floating and not contact.active}

    def capsule_geometries(self) -> List:
        return [electrode.capsule_geometry(thickness=0.1, max_h=0.1)
                for electrode in self.__electrodes]

    def electrode_geometries(self) -> List:
        return [electrode.geometry() for electrode in self.__electrodes]
