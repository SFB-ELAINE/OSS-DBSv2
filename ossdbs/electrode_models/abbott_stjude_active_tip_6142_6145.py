# Abbott/St Jude Active Tip 6142-6145
from ossdbs.electrode_models.electrode import ElectrodeModel
import netgen.occ as occ
import netgen
import numpy as np


class AbbottStjudeActiveTip6142_6145(ElectrodeModel):
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
    TIP_LENGTH = 2.6
    CONTACT_LENGTH = 1.5
    CONTACT_SPACING = 0.5
    LEAD_DIAMETER = 1.3
    TOTAL_LENGHTH = 100.0

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

    def geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of electrode.

        Returns
        -------
        geometry : netgen.libngpy._NgOCC.TopoDS_Shape
        """
        contacts = self.__contacts()
        electrode = occ.Glue([self.__body() - contacts, contacts])
        return electrode.Move(v=self.__position)

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self.LEAD_DIAMETER * 0.5
        point = tuple(np.array(self.__direction) * radius)
        height = self.TOTAL_LENGHTH - self.TIP_LENGTH
        body = occ.Cylinder(p=point, d=self.__direction, r=radius, h=height)
        body.bc(self.__boundaries['Body'])
        return body

    def __contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self.LEAD_DIAMETER * 0.5
        direction = self.__direction

        center = tuple(np.array(direction) * radius)
        contact_tip = occ.Sphere(c=center, r=radius)
        h_pt2 = self.TIP_LENGTH - radius
        contact_pt2 = occ.Cylinder(p=center, d=direction, r=radius, h=h_pt2)
        contact_1 = contact_tip + contact_pt2

        height = self.CONTACT_LENGTH
        contact = occ.Cylinder(p=(0, 0, 0), d=direction, r=radius, h=height)

        distance_1 = self.TIP_LENGTH + self.CONTACT_SPACING
        distance_2 = distance_1 + self.CONTACT_LENGTH + self.CONTACT_SPACING
        distance_3 = distance_2 + self.CONTACT_LENGTH + self.CONTACT_SPACING

        vector_1 = tuple(np.array(direction) * distance_1)
        vector_2 = tuple(np.array(direction) * distance_2)
        vector_3 = tuple(np.array(direction) * distance_3)

        contact_2 = contact.Move(v=vector_1)
        contact_3 = contact.Move(v=vector_2)
        contact_4 = contact.Move(v=vector_3)

        contacts = [contact_1, contact_2, contact_3, contact_4]

        for index, contact in enumerate(contacts, 1):
            name = self.__boundaries['Contact_{}'.format(index)]
            contact.bc(name)
            for edge in contact.edges:
                edge.name = name

        return occ.Glue(contacts)

    def capsule_geometry(self, thickness: float, max_h: float = 0.1) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:
        center = tuple(np.array(self.__direction) * self.LEAD_DIAMETER * 0.5)
        radius = self.LEAD_DIAMETER * 0.5 + thickness
        height = self.TOTAL_LENGHTH - self.TIP_LENGTH
        tip = occ.Sphere(c=center, r=radius)
        lead = occ.Cylinder(p=center, d=self.__direction, r=radius, h=height)
        capsule = tip + lead
        capsule.bc('Capsule')
        capsule.mat('Capsule')
        capsule.maxh = max_h
        return capsule.Move(v=self.__position)
