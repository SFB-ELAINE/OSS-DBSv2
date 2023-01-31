# PINS Medical L303
from ossdbs.electrodes.electrode import Electrode
import netgen
import numpy as np
import os
import json


class PINSMedicalCustom(Electrode):
    """PINS Medical L302 electrode.

    Attributes
    ----------
    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    translation : tuple
        Translation vector (x,y,z) of electrode.
    """

    # dimensions [m]
    TIP_LENGTH = 1.1e-3
    CONTACT_LENGTH = 3.0e-3
    CONTACT_SPACING = 3.0e-3
    LEAD_DIAMETER = 1.3e-3
    TOTAL_LENGHTH = 100.0e-3

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        self.__position = position
        norm = np.linalg.norm(direction)
        self.__direction = tuple(direction / norm) if norm else (0, 0, 1)
        self.__load_parameters()
        self.__boundaries = {'Body': 'Body',
                             'Contact_1': 'Contact_1',
                             'Contact_2': 'Contact_2',
                             'Contact_3': 'Contact_3',
                             'Contact_4': 'Contact_4'}

    def rename_boundaries(self, boundaries: dict) -> None:
        self.__boundaries.update(boundaries)

    def generate_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of electrode.

        Returns
        -------
        geometry : netgen.libngpy._NgOCC.TopoDS_Shape
        """
        contacts = self.__contacts()
        body = self.__body() - contacts
        electrode = netgen.occ.Glue([body, contacts])
        return electrode.Move(self.__position)

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self.LEAD_DIAMETER * 0.5
        center = tuple(np.array(self.__direction) * radius)
        tip = netgen.occ.Sphere(c=center, r=radius)
        lead = netgen.occ.Cylinder(p=center,
                                   d=self.__direction,
                                   r=radius,
                                   h=self.TOTAL_LENGHTH - self.TIP_LENGTH)
        body = tip + lead
        body.bc(self.__boundaries['Body'])
        return body

    def __contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        contact = netgen.occ.Cylinder(p=(0, 0, 0),
                                      d=self.__direction,
                                      r=self.LEAD_DIAMETER * 0.5,
                                      h=self.CONTACT_LENGTH)

        length = (self.CONTACT_LENGTH + self.CONTACT_SPACING)
        n_contacts = 4
        distrances = np.arange(n_contacts) * length + self.TIP_LENGTH
        contacts = [contact.Move(tuple(np.array(self.__direction) * distance))
                    for distance in distrances]

        for index, contact in enumerate(contacts, 1):
            contact.bc(self.__boundaries['Contact_{}'.format(index)])

        return netgen.occ.Glue(contacts)

    def __load_parameters(self) -> None:
        dir_name = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(dir_name, 'pins_medical_custom.json')

        with open(path, 'r') as json_file:
            parameters = json.load(json_file)

        self.CONTACT_LENGTH = parameters['ContactLength[m]']
        self.TIP_LENGTH = parameters['TipLength[m]']
        self.CONTACT_SPACING = parameters['ContactSpacingAxial[m]']
        self.LEAD_DIAMETER = parameters['LeadDiameter[m]']
        self.TOTAL_LENGHTH = parameters['LeadLength[m]']
