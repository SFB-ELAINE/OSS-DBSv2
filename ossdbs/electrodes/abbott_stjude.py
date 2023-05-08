# Abbott/St Jude Active Tip 6142-6145
from dataclasses import dataclass
from .electrode_model_template import ElectrodeModel
import netgen.occ as occ
import netgen
import numpy as np


@dataclass
class AbbottStJudeParameters():
    # dimensions [mm]
    tip_length: float
    contact_length: float
    contact_spacing: float
    lead_diameter: float
    total_length: float


class AbbottStJudeActiveTipModel(ElectrodeModel):
    """Abbott/St Jude Active Tip electrode.

    Attributes
    ----------
    parameters : AbbottStJudeParameters
        Parameters for Abbott Stjude Geometry.

    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    position : tuple
        Position vector (x,y,z) of electrode tip.
    """

    _n_contacts = 4

    def _construct_encapsulation_geometry(self, thickness: float) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:
        center = tuple(np.array(self._direction) * self._parameters.lead_diameter * 0.5)
        radius = self._parameters.lead_diameter * 0.5 + thickness
        height = self._parameters.total_length - self._parameters.tip_length
        tip = occ.Sphere(c=center, r=radius)
        lead = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)
        encapsulation = tip + lead
        encapsulation.bc('EncapsulationLayerSurface')
        encapsulation.mat('EncapsulationLayer')
        return encapsulation.Move(v=self._position) - self.geometry

    def _construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        contacts = self.__contacts()
        electrode = occ.Glue([self.__body() - contacts, contacts])
        return electrode.Move(v=self._position)

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self._parameters.lead_diameter * 0.5
        point = tuple(np.array(self._direction) * radius)
        height = self._parameters.total_length - self._parameters.tip_length
        body = occ.Cylinder(p=point, d=self._direction, r=radius, h=height)
        body.bc(self._boundaries['Body'])
        return body

    def __contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self._parameters.lead_diameter * 0.5
        direction = self._direction

        center = tuple(np.array(direction) * radius)
        contact_tip = occ.Sphere(c=center, r=radius)
        h_pt2 = self._parameters.tip_length - radius
        contact_pt2 = occ.Cylinder(p=center, d=direction, r=radius, h=h_pt2)
        contact_1 = contact_tip + contact_pt2

        height = self._parameters.contact_length
        contact = occ.Cylinder(p=(0, 0, 0), d=direction, r=radius, h=height)

        # TODO shorten
        distance_1 = self._parameters.tip_length + self._parameters.contact_spacing
        distance_2 = distance_1 + self._parameters.contact_length + self._parameters.contact_spacing
        distance_3 = distance_2 + self._parameters.contact_length + self._parameters.contact_spacing

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


class AbbottStJudeDirectedModel(ElectrodeModel):
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

    _n_contacts = 8

    def _construct_encapsulation_geometry(self, thickness: float) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:
        center = tuple(np.array(self._direction) * self._parameters.lead_diameter * 0.5)
        radius = self._parameters.lead_diameter * 0.5 + thickness
        height = self._parameters.total_length - self._parameters.tip_length
        tip = occ.Sphere(c=center, r=radius)
        lead = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)
        encapsulation = tip + lead
        encapsulation.mat('EncapsulationLayer')
        return encapsulation.Move(v=self._position) - self.geometry

    def _construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        contacts = self.__contacts()
        # TODO check
        body = self.__body() - contacts
        electrode = occ.Glue([body, contacts])
        axis = occ.Axis(p=(0, 0, 0), d=self._direction)
        rotated_electrode = electrode.Rotate(axis=axis, ang=self._rotation)
        return rotated_electrode.Move(v=self._position)

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self._parameters.lead_diameter * 0.5
        center = tuple(np.array(self._direction) * radius)
        tip = occ.Sphere(c=center, r=radius)
        height = self._parameters.total_length - self._parameters.tip_length
        lead = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)
        body = tip + lead
        body.bc(self._boundaries['Body'])
        return body

    def __contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:

        # TODO shorten
        distance_1 = self._parameters.tip_length
        distance_2 = distance_1 + self._parameters.contact_length + self._parameters.contact_spacing
        distance_3 = distance_2 + self._parameters.contact_length + self._parameters.contact_spacing
        distance_4 = distance_3 + self._parameters.contact_length + self._parameters.contact_spacing

        vector_1 = tuple(np.array(self._direction) * distance_1)
        vector_2 = tuple(np.array(self._direction) * distance_2)
        vector_3 = tuple(np.array(self._direction) * distance_3)
        vector_4 = tuple(np.array(self._direction) * distance_4)

        point = (0, 0, 0)
        radius = self._parameters.lead_diameter * 0.5
        height = self._parameters.contact_length
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
        radius = self._parameters.lead_diameter * 0.5
        height = self._parameters.contact_length
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
