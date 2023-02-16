# PINS Medical L302
from ossdbs.electrodes.electrode import Electrode
import netgen
import numpy as np


class PINSMedicalL302(Electrode):
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
    CONTACT_LENGTH = 1.5e-3
    CONTACT_SPACING = 1.5e-3
    LEAD_DIAMETER = 1.3e-3
    TOTAL_LENGHTH = 100.0e-3

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        self.__position = position
        norm = np.linalg.norm(direction)
        self.__direction = tuple(direction / norm) if norm else (0, 0, 1)
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

        distance_1 = self.TIP_LENGTH
        distance_2 = distance_1 + self.CONTACT_LENGTH + self.CONTACT_SPACING
        distance_3 = distance_2 + self.CONTACT_LENGTH + self.CONTACT_SPACING
        distance_4 = distance_3 + self.CONTACT_LENGTH + self.CONTACT_SPACING

        vector_1 = tuple(np.array(self.__direction) * distance_1)
        vector_2 = tuple(np.array(self.__direction) * distance_2)
        vector_3 = tuple(np.array(self.__direction) * distance_3)
        vector_4 = tuple(np.array(self.__direction) * distance_4)

        contacts = [contact.Move(vector_1),
                    contact.Move(vector_2),
                    contact.Move(vector_3),
                    contact.Move(vector_4)]

        for index, contact in enumerate(contacts, 1):
            contact.bc(self.__boundaries['Contact_{}'.format(index)])

            for edge in contact.edges:
                edge.name = self.__boundaries['Contact_{}'.format(index)]
        return netgen.occ.Glue(contacts)

    def encapsulating_geometry(self, thickness: float) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self.LEAD_DIAMETER * 0.5
        center = tuple(np.array(self.__direction) * radius)
        lead = netgen.occ.Cylinder(p=center,
                                   d=self.__direction,
                                   r=radius + thickness,
                                   h=self.TOTAL_LENGHTH - self.TIP_LENGTH)
        tip = netgen.occ.Sphere(c=center, r=radius + thickness,)
        capsule = tip + lead
        capsule.bc('Capsule')
        capsule.mat('Capsule')
        capsule.maxh = 0.0001
        return capsule.Move(self.__position)
