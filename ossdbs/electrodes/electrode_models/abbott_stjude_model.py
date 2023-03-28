# Abbott/St Jude Active Tip 6142-6145
from dataclasses import dataclass
from .electrode_model import ElectrodeModel
import netgen.occ as occ
import netgen
import numpy as np


@dataclass
class AbbottStjudeParameters():
    # dimensions [mm]
    tip_length: float
    contact_length: float
    contact_spacing: float
    lead_diameter: float
    total_length: float


class AbbottStjudeActiveTipModel(ElectrodeModel):
    """Abbott/St Jude Active Tip electrode.

    Attributes
    ----------
    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    position : tuple
        Position vector (x,y,z) of electrode tip.
    """

    def __init__(self,
                 parameters: AbbottStjudeParameters,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)
                 ) -> None:
        self._position = position
        norm = np.linalg.norm(direction)
        self._direction = tuple(direction / norm) if norm else (0, 0, 1)
        self._boundaries = {'Body': 'Body',
                             'Contact_1': 'Contact_1',
                             'Contact_2': 'Contact_2',
                             'Contact_3': 'Contact_3',
                             'Contact_4': 'Contact_4'}
        self._tip_length = parameters.tip_length
        self._contact_length = parameters.contact_length
        self._contact_spacing = parameters.contact_spacing
        self._lead_diameter = parameters.lead_diameter
        self._total_length = parameters.total_length

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
        center = tuple(np.array(self._direction) * self._lead_diameter * 0.5)
        radius = self._lead_diameter * 0.5 + thickness
        height = self._total_length - self._tip_length
        tip = occ.Sphere(c=center, r=radius)
        lead = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)
        capsule = tip + lead
        capsule.bc('Capsule')
        capsule.mat('Capsule')
        return capsule.Move(v=self._position) - self.geometry()

    def geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of electrode.

        Returns
        -------
        netgen.libngpy._NgOCC.TopoDS_Shape
        """
        contacts = self.__contacts()
        electrode = occ.Glue([self.__body() - contacts, contacts])
        return electrode.Move(v=self._position)

    def set_contact_names(self, boundaries: dict) -> None:
        """Set the names of electrode contacts.

        Parameters
        ----------
        contact_names : dict
            {'Body': 'body_name',
             'Contact_1': 'contact_name',
             'Contact_2': ...}
        """
        self._boundaries.update(boundaries)

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self._lead_diameter * 0.5
        point = tuple(np.array(self._direction) * radius)
        height = self._total_length - self._tip_length
        body = occ.Cylinder(p=point, d=self._direction, r=radius, h=height)
        body.bc(self._boundaries['Body'])
        return body

    def __contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self._lead_diameter * 0.5
        direction = self._direction

        center = tuple(np.array(direction) * radius)
        contact_tip = occ.Sphere(c=center, r=radius)
        h_pt2 = self._tip_length - radius
        contact_pt2 = occ.Cylinder(p=center, d=direction, r=radius, h=h_pt2)
        contact_1 = contact_tip + contact_pt2

        height = self._contact_length
        contact = occ.Cylinder(p=(0, 0, 0), d=direction, r=radius, h=height)

        distance_1 = self._tip_length + self._contact_spacing
        distance_2 = distance_1 + self._contact_length + self._contact_spacing
        distance_3 = distance_2 + self._contact_length + self._contact_spacing

        vector_1 = tuple(np.array(direction) * distance_1)
        vector_2 = tuple(np.array(direction) * distance_2)
        vector_3 = tuple(np.array(direction) * distance_3)

        contact_2 = contact.Move(v=vector_1)
        contact_3 = contact.Move(v=vector_2)
        contact_4 = contact.Move(v=vector_3)

        contacts = [contact_1, contact_2, contact_3, contact_4]

        for index, contact in enumerate(contacts, 1):
            name = self._boundaries['Contact_{}'.format(index)]
            contact.bc(name)
            for edge in contact.edges:
                edge.name = name

        return occ.Glue(contacts)


class AbbottStjudeDirectedModel(ElectrodeModel):
    """Abbott/St Jude Directed electrode.

    Attributes
    ----------
    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    position : tuple
        Position vector (x,y,z) of electrode tip.
    """

    def __init__(self,
                 parameters: AbbottStjudeParameters,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)
                 ) -> None:
        self._position = position
        self._rotation = rotation
        norm = np.linalg.norm(direction)
        self._direction = tuple(direction / norm) if norm else (0, 0, 1)
        self._boundaries = {'Body': 'Body',
                             'Contact_1': 'Contact_1',
                             'Contact_2': 'Contact_2',
                             'Contact_3': 'Contact_3',
                             'Contact_4': 'Contact_4',
                             'Contact_5': 'Contact_5',
                             'Contact_6': 'Contact_6',
                             'Contact_7': 'Contact_7',
                             'Contact_8': 'Contact_8'}
        self._tip_length = parameters.tip_length
        self._contact_length = parameters.contact_length
        self._contact_spacing = parameters.contact_spacing
        self._lead_diameter = parameters.lead_diameter
        self._total_length = parameters.total_length

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
        center = tuple(np.array(self._direction) * self._lead_diameter * 0.5)
        radius = self._lead_diameter * 0.5 + thickness
        height = self._total_length - self._tip_length
        tip = occ.Sphere(c=center, r=radius)
        lead = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)
        capsule = tip + lead
        capsule.bc('Capsule')
        capsule.mat('Capsule')
        return capsule.Move(v=self._position) - self.geometry()

    def geometry(self) -> netgen.libngpy._meshing.Mesh:
        """Generate geometry of electrode.

        Returns
        -------
        netgen.libngpy._NgOCC.TopoDS_Shape
        """
        contacts = self.__contacts()
        body = self.__body() - contacts
        electrode = occ.Glue([body, contacts])
        axis = occ.Axis(p=(0, 0, 0), d=self._direction)
        rotated_electrode = electrode.Rotate(axis=axis, ang=self._rotation)
        return rotated_electrode.Move(v=self._position)

    def set_contact_names(self, boundaries: dict) -> None:
        """Set the names of electrode contacts.

        Parameters
        ----------
        contact_names : dict
            {'Body': 'body_name',
             'Contact_1': 'contact_name',
             'Contact_2': ...}
        """
        self._boundaries.update(boundaries)

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self._lead_diameter * 0.5
        center = tuple(np.array(self._direction) * radius)
        tip = occ.Sphere(c=center, r=radius)
        height = self._total_length - self._tip_length
        lead = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)
        body = tip + lead
        body.bc(self._boundaries['Body'])
        return body

    def __contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:

        distance_1 = self._tip_length
        distance_2 = distance_1 + self._contact_length + self._contact_spacing
        distance_3 = distance_2 + self._contact_length + self._contact_spacing
        distance_4 = distance_3 + self._contact_length + self._contact_spacing

        vector_1 = tuple(np.array(self._direction) * distance_1)
        vector_2 = tuple(np.array(self._direction) * distance_2)
        vector_3 = tuple(np.array(self._direction) * distance_3)
        vector_4 = tuple(np.array(self._direction) * distance_4)

        point = (0, 0, 0)
        radius = self._lead_diameter * 0.5
        height = self._contact_length
        contact = occ.Cylinder(p=point, d=self._direction, r=radius, h=height)
        axis = occ.Axis(p=point, d=self._direction)

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
            name = self._boundaries['Contact_{}'.format(index)]
            contact.bc(name)
            for edge in contact.edges:
                edge.name = name

        return occ.Fuse(contacts)

    def __contact_directed(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        point = (0, 0, 0)
        radius = self._lead_diameter * 0.5
        height = self._contact_length
        body = occ.Cylinder(p=point, d=self._direction, r=radius, h=height)
        new_direction = tuple(np.cross(self.__direction_2(), self._direction))
        eraser = occ.HalfSpace(p=point, n=new_direction)
        delta = 15
        angle = 30 + delta
        axis = occ.Axis(p=point, d=self._direction)
        return body - eraser.Rotate(axis, angle) - eraser.Rotate(axis, -angle)

    def __direction_2(self):
        x, y, z = self._direction

        if not x and not y:
            return (0, 1, 0)

        if not x and not z:
            return (0, 0, 1)

        if not y and not z:
            return (0, 1, 0)

        return (x, y, not z)
