# Abbott/St Jude Active Tip 6142-6145
from ossdbs.electrodes.electrode import Electrode
import netgen
import numpy as np


class AbbottStjudeActiveTip6142_6145(Electrode):
    """Abbott/St Jude Active Tip 6142-6149 electrode.

    Attributes
    ----------
    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    translation : tuple
        Translation vector (x,y,z) of electrode.
    """

    # dimensions [mm]
    TIP_LENGTH = 1.1
    CONTACT_LENGTH = 1.5
    CONTACT_SPACING = 0.5
    LEAD_DIAMETER = 1.3
    TOTAL_LENGHTH = 100.0
    N_CONTACTS = 4

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 translation: tuple = (0, 0, 0)) -> None:
        self.__translation = translation
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
        return electrode.Move(self.__translation)

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
        distrances = np.arange(1, self.N_CONTACTS) * length + offset
        contacts = [active_tip] + \
                   [contact.Move(tuple(np.array(self.__direction) * distance))
                    for distance in distrances]

        for index, contact in enumerate(contacts, 1):
            contact.bc(self.__boundaries['Contact_{}'.format(index)])

        return netgen.occ.Fuse(contacts)
