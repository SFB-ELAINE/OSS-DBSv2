# Abbott/St Jude Active Tip 6142-6145
from src.electrodes.abstract_electrode import AbstractElectrode
import netgen
import numpy as np


class AbbottStjudeActiveTip6142_6145(AbstractElectrode):
    """Abbott/St Jude Active Tip 6142-6149 electrode.

    Attributes
    ----------
    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    translation : tuple
        Translation vector (x,y,z) of electrode.

    Methods
    -------
    generate_geometry()
        Generate geometry of electrode.
    """

    # dimensions [mm]
    TIP_LENGTH = 1.1
    CONTACT_LENGTH = 1.5
    CONTACT_SPACING = 0.5
    LEAD_DIAMETER = 1.3
    TOTAL_LENGHTH = 100.0
    TUBE_THICKNESS = 0.0
    TUBE_FREE_LENGTH = 100.0
    N_CONTACTS = 4

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 translation: tuple = (0, 0, 0)) -> None:
        self.__translation = translation
        norm = np.linalg.norm(direction)
        self.__direction = tuple(direction / norm) if norm else (0, 0, 1)

    def generate_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of electrode.

        Returns
        -------
        geometry : netgen.libngpy._NgOCC.TopoDS_Shape
        """
        return self.__construct_geometry()

    def __construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        tube = self.__tube()

        contacts = self.__contacts(diameter=self.LEAD_DIAMETER)

        body = self.__body(diameter=self.LEAD_DIAMETER) - contacts
        electrode = netgen.occ.Glue([body, contacts]) - tube

        diameter = self.LEAD_DIAMETER + 2 * self.TUBE_THICKNESS
        contacts_tube = self.__contacts(diameter=diameter)
        body_tube = self.__body(diameter=diameter) - contacts_tube
        electrode_tube = netgen.occ.Glue([body_tube, contacts_tube]) * tube

        electrode = netgen.occ.Glue([electrode, electrode_tube])
        moved_electrode = electrode.Move(self.__translation)

        return moved_electrode

    def __body(self, diameter: float) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        point = tuple(np.array(self.__direction) * self.LEAD_DIAMETER * 0.5)
        height = self.TOTAL_LENGHTH - self.TIP_LENGTH
        body = netgen.occ.Cylinder(p=point,
                                   d=self.__direction,
                                   r=diameter * 0.5,
                                   h=height)
        body.bc("Body")
        return body

    def __contacts(self, diameter: float)\
            -> netgen.libngpy._NgOCC.TopoDS_Shape:

        point = tuple(np.array(self.__direction) * self.TIP_LENGTH)
        space = netgen.occ.HalfSpace(p=point, n=self.__direction)
        center = tuple(np.array(self.__direction) * self.LEAD_DIAMETER * 0.5)
        tip = netgen.occ.Sphere(c=center, r=diameter * 0.5) * space
        height = max(self.TIP_LENGTH - self.LEAD_DIAMETER * 0.5, 0)

        if height > 0:
            contact_1 = netgen.occ.Cylinder(p=center,
                                            d=self.__direction,
                                            r=diameter * 0.5,
                                            h=height)
            active_tip = netgen.occ.Fuse([tip, contact_1])
        else:
            active_tip = tip

        active_tip.bc('Contact_1')

        contact = netgen.occ.Cylinder(p=(0, 0, 0),
                                      d=self.__direction,
                                      r=diameter * 0.5,
                                      h=self.CONTACT_LENGTH)
        length = (self.CONTACT_LENGTH + self.CONTACT_SPACING)
        offset = self.TIP_LENGTH + self.CONTACT_SPACING
        contacts = [contact.Move(tuple(np.array(self.__direction) * distance))
                    for distance
                    in np.arange(self.N_CONTACTS - 1) * length + offset]

        for index, contact in enumerate(contacts, 2):
            contact.bc("Contact_{}".format(index))

        return netgen.occ.Fuse([active_tip] + contacts)

    def __tube(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self.LEAD_DIAMETER * 0.5 + self.TUBE_THICKNESS
        height = self.TOTAL_LENGHTH - self.TUBE_FREE_LENGTH
        point = tuple(np.array(self.__direction) * self.TUBE_FREE_LENGTH)
        tube = netgen.occ.Cylinder(p=point,
                                   d=self.__direction,
                                   r=radius,
                                   h=height)
        tube.bc('Body')

        offset = self.TIP_LENGTH + self.CONTACT_SPACING
        distance = self.CONTACT_LENGTH + self.CONTACT_SPACING
        distances = offset + np.arange(self.N_CONTACTS - 1) * distance
        lower_limit = np.append(0, distances)
        upper_limit = np.append(self.TIP_LENGTH,
                                distances + self.CONTACT_LENGTH)
        intersection = np.logical_and(lower_limit < self.TUBE_FREE_LENGTH,
                                      upper_limit > self.TUBE_FREE_LENGTH)
        if np.any(intersection):
            tube.bc('Contact_{}'.format(np.argmax(intersection) + 1))

        return tube
