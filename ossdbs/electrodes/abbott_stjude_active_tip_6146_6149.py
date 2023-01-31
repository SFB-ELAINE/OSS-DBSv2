# Abbott/St Jude Active Tip 6146-6149
from ossdbs.electrodes.electrode import Electrode
import netgen
import numpy as np


class AbbottStjudeActiveTip6146_6149(Electrode):
    """Abbott/St Jude Active Tip 6146-6149 electrode.

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
        point = tuple(np.array(self.__direction) * radius)
        body = netgen.occ.Cylinder(p=point,
                                   d=self.__direction,
                                   r=radius,
                                   h=self.TOTAL_LENGHTH - self.TIP_LENGTH)
        body.bc(self.__boundaries['Body'])
        return body

    def __contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self.LEAD_DIAMETER * 0.5
        center = tuple(np.array(self.__direction) * radius)
        tip = netgen.occ.Sphere(c=center, r=radius)
        contact_1 = netgen.occ.Cylinder(p=center,
                                        d=self.__direction,
                                        r=radius,
                                        h=self.TIP_LENGTH + radius)
        active_tip = tip + contact_1

        contact = netgen.occ.Cylinder(p=(0, 0, 0),
                                      d=self.__direction,
                                      r=radius,
                                      h=self.CONTACT_LENGTH)
        length = (self.CONTACT_LENGTH + self.CONTACT_SPACING)
        offset = self.TIP_LENGTH + self.CONTACT_SPACING
        n_contacts = 4
        distrances = np.arange(1, n_contacts) * length + offset
        contacts = [active_tip] + \
                   [contact.Move(tuple(np.array(self.__direction) * distance))
                    for distance in distrances]

        for index, contact in enumerate(contacts, 1):
            contact.bc(self.__boundaries['Contact_{}'.format(index)])

        return netgen.occ.Fuse(contacts)
