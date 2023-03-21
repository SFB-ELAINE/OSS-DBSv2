
from typing import List
from ossdbs.electrodes.contacts import Contact, Contacts
from ossdbs.electrodes.electrode_models import ElectrodeModel
from ossdbs.electrodes.electrode_models import AbbottStjudeActiveTip6142_6145
from ossdbs.electrodes.electrode_models import AbbottStjudeActiveTip6146_6149
from ossdbs.electrodes.electrode_models import AbbottStjudeDirected6172
from ossdbs.electrodes.electrode_models import AbbottStjudeDirected6173
from ossdbs.electrodes.electrode_models import BostonScientificVercise
from ossdbs.electrodes.electrode_models import BostonScientificVerciseDirected
from ossdbs.electrodes.electrode_models import Medtronic3387
from ossdbs.electrodes.electrode_models import Medtronic3389
from ossdbs.electrodes.electrode_models import Medtronic3391
from ossdbs.electrodes.electrode_models import PINSMedicalL301
from ossdbs.electrodes.electrode_models import PINSMedicalL302
from ossdbs.electrodes.electrode_models import PINSMedicalL303
from ossdbs.electrodes.electrode_models import MicroProbesCustomRodent
from ossdbs.electrodes.electrode_models import MicroProbesSNEX_100
from ossdbs.electrodes.electrode_models import AbbottStjudeActiveCustom
from ossdbs.electrodes.electrode_models import AbbottStjudeDirectedCustom
from ossdbs.electrodes.electrode_models import BostonScientificVerciseCustom
from ossdbs.electrodes.electrode_models import BostonScientificVerciseDirectedCustom
from ossdbs.electrodes.electrode_models import MedtronicCustom
from ossdbs.electrodes.electrode_models import MicroProbesCustomRodentCustom
from ossdbs.electrodes.electrode_models import MicroProbesSNEX_100Custom
from ossdbs.electrodes.electrode_models import PINSMedicalCustom

from ossdbs.electrodes.electrodes import Electrodes


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
               rotation: float
               ) -> ElectrodeModel:
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


class ElectrodesFactory:

    @classmethod
    def create(cls, electrode_paramters: dict) -> Electrodes:

        electrodes = []
        contacts = []

        for index, input_par in enumerate(electrode_paramters):
            electrode = cls.__create_electrode(input_par)
            electrode.set_contact_names(cls.__boundaries(index))
            electrodes.append(electrode)
            contacts.extend(cls.__create_contacts(index, input_par))

        contacts = Contacts([contact for contact in contacts
                             if contact.active or contact.floating])

        return Electrodes(electrodes=electrodes, contacts=contacts)

    @classmethod
    def __create_contacts(cls,
                          index: int,
                          electrode_parameters: dict
                          ) -> List[Contact]:
        return [cls.create_contact(index, contact_par)
                for contact_par in electrode_parameters['Contacts']]

    @staticmethod
    def create_contact(index: int, contact_par: dict) -> Contact:
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

    @staticmethod
    def __create_electrode(electrode_parameters: dict) -> ElectrodeModel:
        dir = electrode_parameters['Direction']
        pos = electrode_parameters['TipPosition']
        direction = (dir['x[mm]'], dir['y[mm]'], dir['z[mm]'])
        position = (pos['x[mm]'], pos['y[mm]'], pos['z[mm]'])
        rotation = electrode_parameters['Rotation[Degrees]']
        return ElectrodeFactory.create(name=electrode_parameters['Name'],
                                       direction=direction,
                                       position=position,
                                       rotation=rotation)

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
