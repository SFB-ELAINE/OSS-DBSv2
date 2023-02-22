
from ossdbs.electrodes import Electrode
from ossdbs.electrodes import AbbottStjudeActiveTip6142_6145
from ossdbs.electrodes import AbbottStjudeActiveTip6146_6149
from ossdbs.electrodes import AbbottStjudeDirected6172
from ossdbs.electrodes import AbbottStjudeDirected6173
from ossdbs.electrodes import BostonScientificVercise
from ossdbs.electrodes import BostonScientificVerciseDirected
from ossdbs.electrodes import Medtronic3387
from ossdbs.electrodes import Medtronic3389
from ossdbs.electrodes import Medtronic3391
from ossdbs.electrodes import PINSMedicalL301
from ossdbs.electrodes import PINSMedicalL302
from ossdbs.electrodes import PINSMedicalL303
from ossdbs.electrodes import MicroProbesCustomRodent
from ossdbs.electrodes import MicroProbesSNEX_100

from ossdbs.electrodes import AbbottStjudeActiveCustom
from ossdbs.electrodes import AbbottStjudeDirectedCustom
from ossdbs.electrodes import BostonScientificVerciseCustom
from ossdbs.electrodes import BostonScientificVerciseDirectedCustom
from ossdbs.electrodes import MedtronicCustom
from ossdbs.electrodes import MicroProbesCustomRodentCustom
from ossdbs.electrodes import MicroProbesSNEX_100Custom
from ossdbs.electrodes import PINSMedicalCustom

from dataclasses import dataclass
from typing import List


class ElectrodeFactory:
    """Creates a list of Electrode objects."""

    ELECTRODES = {'AbbottStjudeActiveTip6142_6145':
                  AbbottStjudeActiveTip6142_6145,
                  'AbbottStjudeActiveTip6146_6149':
                  AbbottStjudeActiveTip6146_6149,
                  'AbbottStjudeDirected6172':
                  AbbottStjudeDirected6172,
                  'AbbottStjudeDirected6173':
                  AbbottStjudeDirected6173,
                  'BostonScientificVercise':
                  BostonScientificVercise,
                  'BostonScientificVerciseDirected':
                  BostonScientificVerciseDirected,
                  'Medtronic3387':
                  Medtronic3387,
                  'Medtronic3389':
                  Medtronic3389,
                  'Medtronic3391':
                  Medtronic3391,
                  'MicroProbesSNEX_100':
                  MicroProbesSNEX_100,
                  'PINSMedicalL301':
                  PINSMedicalL301,
                  'PINSMedicalL302':
                  PINSMedicalL302,
                  'PINSMedicalL303':
                  PINSMedicalL303,
                  'MicroProbesCustomRodent':
                  MicroProbesCustomRodent,
                  'AbbottStjudeActiveCustom':
                  AbbottStjudeActiveCustom,
                  'AbbottStjudeDirectedCustom':
                  AbbottStjudeDirectedCustom,
                  'BostonScientificVerciseCustom':
                  BostonScientificVerciseCustom,
                  'BostonScientificVerciseDirectedCustom':
                  BostonScientificVerciseDirectedCustom,
                  'MedtronicCustom':
                  MedtronicCustom,
                  'MicroProbesCustomRodentCustom':
                  MicroProbesCustomRodentCustom,
                  'MicroProbesSNEX_100Custom':
                  MicroProbesSNEX_100Custom,
                  'PINSMedicalCustom':
                  PINSMedicalCustom
                  }

    @classmethod
    def create(cls,
               name: str,
               direction: tuple,
               position: tuple,
               rotation: float) -> Electrode:
        """create a list of Electrode objects.

        Parameters
        ----------
        parameters : ElectrodeParameters

        Returns
        -------
        Electrode
            Electrode objects.
        """

        electrode_type = cls.ELECTRODES[name]
        return electrode_type(direction=direction,
                              position=position,
                              rotation=rotation)


@dataclass
class ElectrodeContact:
    name: str = ''
    active: bool = False
    current: float = 0.0
    floating: bool = False
    voltage: float = 0.0
    surface_impedance: float = 0.0


class Electrodes:

    def __init__(self, electrode_paramters: dict) -> None:

        self.__electrodes = []
        self.__contacts = []

        for index, input_par in enumerate(electrode_paramters):
            electrode = self.__create_electrode(input_par)
            boundary_names = self.__boundaries(index)
            electrode.rename_boundaries(boundary_names)
            electrode_contacts = self.__create_contacts(index, input_par)
            self.__electrodes.append(electrode)
            self.__contacts.extend(electrode_contacts)

    def __create_contacts(self, index, electrode):
        contacts = []
        for contact_par in electrode['Contacts']:
            name = 'E{}C{}'.format(index, contact_par['Contact_ID'] - 1)
            active = contact_par['Active']
            floating = contact_par['Floating'] and not active
            current = contact_par['Current[A]']
            voltage = contact_par['Voltage[V]']
            real = contact_par['SurfaceImpedance[Ωm]']['real']
            imag = contact_par['SurfaceImpedance[Ωm]']['imag']
            surface_impedance = real + 1j * imag
            contacts.append(ElectrodeContact(name=name,
                                             active=active,
                                             floating=floating,
                                             current=current,
                                             voltage=voltage,
                                             surface_impedance=\
                                             surface_impedance))
        return contacts

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

    def contacts(self) -> List[ElectrodeContact]:
        return [contact for contact in self.__contacts
                if contact.active or contact.floating]

    def electrodes(self) -> List[Electrode]:
        return self.__electrodes

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
