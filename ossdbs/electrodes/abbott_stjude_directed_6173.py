# Abbott/St Jude Directed 6173
from ossdbs.electrodes.abstract_electrode import Electrode
import numpy as np
import netgen


class AbbottStjudeDirected6173(Electrode):
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
    TIP_LENGTH = 1.5
    CONTACT_LENGTH = 1.5
    CONTACT_SPACING = 1.5
    LEAD_DIAMETER = 1.29
    TOTAL_LENGTH = 40.0 # Updated
    CONTACT_SPACING_RADIAL = 0.337721210261 # UPDATED to 30°

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 translation: tuple = (0, 0, 0)) -> None:
        self.__translation = translation
        self.__rotation = rotation
        norm = np.linalg.norm(direction)
        self.__direction = tuple(direction / norm) if norm else (0, 0, 1)
        self.__boundaries = {'Body': 'Body',
                             'Contact_1': 'Contact_1',
                             'Contact_2': 'Contact_2',
                             'Contact_3': 'Contact_3',
                             'Contact_4': 'Contact_4',
                             'Contact_5': 'Contact_5',
                             'Contact_6': 'Contact_6',
                             'Contact_7': 'Contact_7',
                             'Contact_8': 'Contact_8'}

    def rename_boundaries(self, boundaries: dict) -> None:
        self.__boundaries.update(boundaries)

    def generate_geometry(self) -> netgen.libngpy._meshing.Mesh:
        """Generate geometry of electrode.

        Returns
        -------
        geometry : netgen.libngpy._NgOCC.TopoDS_Shape
        """
        contacts = self.__contacts()
        body = self.__body() - contacts
        electrode = netgen.occ.Glue([body, contacts])
        axis = netgen.occ.Axis(p=(0, 0, 0), d=self.__direction)
        rotated_electrode = electrode.Rotate(axis=axis, ang=self.__rotation)
        return rotated_electrode.Move(self.__translation)

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self.LEAD_DIAMETER * 0.5
        center = tuple(np.array(self.__direction) * radius)
        tip = netgen.occ.Sphere(c=center, r=radius)
        lead = netgen.occ.Cylinder(p=center,
                                   d=self.__direction,
                                   r=radius,
                                   h=self.TOTAL_LENGTH - self.TIP_LENGTH)
        body = tip + lead
        body.bc(self.__boundaries['Body'])
        return body

    def __contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:

        axis = netgen.occ.Axis(p=(0, 0, 0), d=self.__direction)

        distance = self.CONTACT_LENGTH + self.CONTACT_SPACING
        distance_2 = self.TIP_LENGTH + distance
        distance_3 = self.TIP_LENGTH + 2 * distance
        distance_4 = self.TIP_LENGTH + 3 * distance

        vector_1 = tuple(np.array(self.__direction) * self.TIP_LENGTH)
        vector_2 = tuple(np.array(self.__direction) * distance_2)
        vector_3 = tuple(np.array(self.__direction) * distance_3)
        vector_4 = tuple(np.array(self.__direction) * distance_4)

        contact = netgen.occ.Cylinder(p=(0, 0, 0),
                                      d=self.__direction,
                                      r=self.LEAD_DIAMETER * 0.5,
                                      h=self.CONTACT_LENGTH)

        contact_directed = self.__contact_directed()

        contacts = [contact.Move(vector_1),
                    contact_directed.Move(vector_2),
                    contact_directed.Rotate(axis, 120).Move(vector_2),
                    contact_directed.Rotate(axis, 240).Move(vector_2),
                    contact_directed.Move(vector_3),
                    contact_directed.Rotate(axis, 120).Move(vector_3),
                    contact_directed.Rotate(axis, 240).Move(vector_3),
                    contact.Move(vector_4)
                    ]

        for index, contact in enumerate(contacts, 1):
            contact.bc(self.__boundaries['Contact_{}'.format(index)])

        return netgen.occ.Fuse(contacts)

    def __contact_directed(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        body = netgen.occ.Cylinder(p=(0, 0, 0),
                                   d=self.__direction,
                                   r=self.LEAD_DIAMETER * 0.5,
                                   h=self.CONTACT_LENGTH)
        new_direction = tuple(np.cross(self.__direction_2(), self.__direction))
        eraser = netgen.occ.HalfSpace(p=(0, 0, 0), n=new_direction)
        delta = self.CONTACT_SPACING_RADIAL / self.LEAD_DIAMETER * 180 / np.pi
        angle = 30 + delta
        axis = netgen.occ.Axis((0, 0, 0), self.__direction)
        return body - eraser.Rotate(axis, angle) - eraser.Rotate(axis, -angle)

    def __direction_2(self):
        x, y, z = self.__direction

        if not x and not y:
            return (0, 1, 0)

        if not x and not z:
            return (0, 0, 1)

        if not y and not z:
            return (0, 1, 0)

        return (x, y, not z)
