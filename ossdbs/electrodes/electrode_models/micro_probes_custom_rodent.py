
from .electrode import ElectrodeModel
import netgen
import netgen.occ as occ
import numpy as np


class MicroProbesCustomRodent(ElectrodeModel):
    """MicroProbes Custom Rodent electrode.

    Attributes
    ----------
    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    position : tuple
        Position vector (x,y,z) of electrode tip.
    """

    # dimensions [m]
    CONTACT_LENGTH = 0.1125
    LEAD_DIAMETER = 0.225
    TOTAL_LENGHTH = 13.3
    TUBE_THICKNESS = .01

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        self.__position = tuple(position)
        norm = np.linalg.norm(direction)
        self.__direction = tuple(direction / norm) if norm else (0, 0, 1)
        self.__boundaries = {'Body': 'Body', 'Contact_1': 'Contact_1'}

    def capsule_geometry(self, thickness: float) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of capsule layer around electrode.

        Parameters
        ----------
        thickness : float
            Thickness of encapsulating layer.

        Returns
        -------
        netgen.libngpy._NgOCC.TopoDS_Shape
        """
        center = tuple(np.array(self.__direction) * self.LEAD_DIAMETER * 0.5)
        radius = self.LEAD_DIAMETER * 0.5 + thickness
        height = self.TOTAL_LENGHTH - self.LEAD_DIAMETER * 0.5
        lead = occ.Cylinder(p=center, d=self.__direction, r=radius, h=height)
        tip = occ.Sphere(c=center, r=radius + thickness,)
        capsule = tip + lead
        capsule.bc('Capsule')
        capsule.mat('Capsule')
        return capsule.Move(v=self.__position) - self.geometry()

    def geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of electrode.

        Returns
        -------
        netgen.libngpy._NgOCC.TopoDS_Shape
        """
        contacts = self.__contacts()
        electrode = netgen.occ.Glue([self.__body() - contacts, contacts])
        return electrode.Move(v=self.__position)

    def set_contact_names(self, boundaries: dict) -> None:
        """Set the names of electrode contacts.

        Parameters
        ----------
        contact_names : dict
            {'Body': 'body_name',
             'Contact_1': 'contact_name',
             'Contact_2': ...}
        """
        self.__boundaries.update(boundaries)

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self.LEAD_DIAMETER * 0.5
        direction = self.__direction
        tube_radius = radius + self.TUBE_THICKNESS
        center = tuple(np.array(self.__direction) * radius)
        tip = occ.Sphere(c=center, r=tube_radius)
        height = self.TOTAL_LENGHTH - radius
        lead = occ.Cylinder(p=center, d=direction, r=tube_radius, h=height)
        point = tuple(np.array(direction) * self.CONTACT_LENGTH)
        space = occ.HalfSpace(p=point, n=direction)
        body = tip + lead - space
        body.bc(self.__boundaries['Body'])
        return body

    def __contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self.LEAD_DIAMETER * 0.5
        point = tuple(np.array(self.__direction) * self.CONTACT_LENGTH)
        space = netgen.occ.HalfSpace(p=point, n=self.__direction)
        center = tuple(np.array(self.__direction) * radius)
        tip = occ.Sphere(c=center, r=radius) * space
        height = self.CONTACT_LENGTH - radius
        lead = occ.Cylinder(p=center, d=self.__direction, r=radius, h=height)

        if self.CONTACT_LENGTH <= radius:
            contact = tip
        else:
            contact = tip + lead

        contact.bc(self.__boundaries['Contact_1'])
        for edge in contact.edges:
            edge.name = self.__boundaries['Contact_1']
        return contact
