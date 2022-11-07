# Abbott/St Jude Directed 6173
from src.electrodes.abstract_electrode import AbstractElectrode
import netgen
import numpy as np


class AbbottStjudeDirected6173(AbstractElectrode):
    """Abbott/St Jude Directed 6173 electrode.

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
    CONTACT_SPACING = 1.5
    LEAD_DIAMETER = 1.3
    TOTAL_LENGTH = 20.0
    TUBE_THICKNESS = 0.0
    TUBE_FREE_LENGTH = 20.0
    CONTACT_SPACING_RADIAL = 0.25

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 translation: tuple = (0, 0, 0)) -> None:
        self.__translation = translation
        norm = np.linalg.norm(direction)
        self.__direction = tuple(direction / norm) if norm else (0, 0, 1)

    def generate_geometry(self) -> netgen.libngpy._meshing.Mesh:
        """Generate geometry of electrode.

        Returns
        -------
        geometry : netgen.libngpy._NgOCC.TopoDS_Shape
        """
        return self.__construct_geometry()

    def __construct_geometry(self) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:
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
        center = tuple(np.array(self.__direction) * self.LEAD_DIAMETER * 0.5)
        tip = netgen.occ.Sphere(c=center, r=diameter * 0.5)
        height = self.TOTAL_LENGTH - self.TIP_LENGTH
        lead = netgen.occ.Cylinder(p=center,
                                   d=self.__direction,
                                   r=diameter * 0.5,
                                   h=height)
        body = tip + lead
        body.bc("Body")
        return body

    def __contacts(self, diameter: float) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:

        point = tuple(np.array(self.__direction) * self.TIP_LENGTH)
        contact = netgen.occ.Cylinder(p=(0, 0, 0),
                                      d=self.__direction,
                                      r=diameter * 0.5,
                                      h=self.CONTACT_LENGTH)

        distance_1 = (self.TIP_LENGTH)

        contact_1 = contact.Move(tuple(np.array(self.__direction) *
                                       distance_1))

        distance_8 = (self.TIP_LENGTH +
                      3 * self.CONTACT_SPACING +
                      3 * self.CONTACT_LENGTH)

        contact_8 = contact.Move(tuple(np.array(self.__direction) *
                                       distance_8))

        axis = netgen.occ.Axis((0, 0, 0), self.__direction)
        distance = self.TIP_LENGTH + self.CONTACT_SPACING + self.CONTACT_LENGTH
        point = tuple(np.array(self.__direction) * distance)
        contact = self.__contact(diameter, self.CONTACT_LENGTH).Move(point)
        distance = self.CONTACT_LENGTH + self.CONTACT_SPACING
        point = tuple(np.array(self.__direction) * distance)
        contacts = [contact_1,
                    contact,
                    contact.Rotate(axis, 120),
                    contact.Rotate(axis, 240),
                    contact.Move(point),
                    contact.Rotate(axis, 120).Move(point),
                    contact.Rotate(axis, 240).Move(point),
                    contact_8
                    ]

        for index, contact in enumerate(contacts, 1):
            contact.bc("Contact_{}".format(index))

        return netgen.occ.Fuse(contacts)

    def __contact(self, diameter: float, height: float) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = diameter * 0.5
        body = netgen.occ.Cylinder(p=(0, 0, 0),
                                   d=self.__direction,
                                   r=radius,
                                   h=height)

        direction2 = self.__secound_direction()
        new_direction = tuple(np.cross(direction2, self.__direction))
        eraser = netgen.occ.HalfSpace(p=(0, 0, 0), n=new_direction)
        delta = self.CONTACT_SPACING_RADIAL / self.LEAD_DIAMETER * 180 / np.pi
        angle = 30 + delta
        axis = netgen.occ.Axis((0, 0, 0), self.__direction)
        return body - eraser.Rotate(axis, angle) - eraser.Rotate(axis, -angle)

    def __secound_direction(self):
        x, y, z = self.__direction

        if not x and not y:
            return (0, 1, 0)

        if not x and not z:
            return (0, 0, 1)

        if not y and not z:
            return (0, 1, 0)

        return (x, y, not z)

    def __tube(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:

        radius = self.LEAD_DIAMETER * 0.5 + self.TUBE_THICKNESS
        height = self.TOTAL_LENGTH - self.TUBE_FREE_LENGTH
        point = tuple(np.array(self.__direction) * self.TUBE_FREE_LENGTH)
        tube = netgen.occ.Cylinder(p=point,
                                   d=self.__direction,
                                   r=radius,
                                   h=height)

        offset = self.TIP_LENGTH
        distance = self.CONTACT_LENGTH + self.CONTACT_SPACING
        lower_limit = offset + np.arange(4) * distance
        upper_limit = lower_limit + self.CONTACT_LENGTH
        intersection = np.logical_and(lower_limit < self.TUBE_FREE_LENGTH,
                                      upper_limit > self.TUBE_FREE_LENGTH)
        index = np.argmax(intersection)

        if not np.any(intersection):
            index = -1

        if index == 0:
            tube.bc('Contact_1')

        elif index == 1:
            tube = self.tube_with_multiple_contacts(first_contact=2)

        elif index == 2:
            tube = self.tube_with_multiple_contacts(first_contact=5)

        elif index == 3:
            tube.bc('Contact_8')
        else:
            tube.bc('Body')

        return tube

    def tube_with_multiple_contacts(self, first_contact: int):
        axis = netgen.occ.Axis((0, 0, 0), self.__direction)
        diameter = self.LEAD_DIAMETER + 2 * self.TUBE_THICKNESS
        height = self.TOTAL_LENGTH - self.TUBE_FREE_LENGTH
        point = tuple(np.array(self.__direction) * self.TUBE_FREE_LENGTH)
        contact = self.__contact(diameter, height).Move(point)
        contacts = [contact,
                    contact.Rotate(axis, 120),
                    contact.Rotate(axis, 240),
                    ]

        for index, contact in enumerate(contacts, first_contact):
            contact.bc("Contact_{}".format(index))

        tube = netgen.occ.Glue(contacts)
        radius = self.LEAD_DIAMETER * 0.5 + self.TUBE_THICKNESS
        height = self.TOTAL_LENGTH - self.TUBE_FREE_LENGTH
        point = tuple(np.array(self.__direction) * self.TUBE_FREE_LENGTH)
        body = netgen.occ.Cylinder(p=point,
                                   d=self.__direction,
                                   r=radius,
                                   h=height)

        body.bc('Body')
        tube = netgen.occ.Glue(contacts)
        body = body - tube
        return netgen.occ.Glue([tube, body])
