# Medtronic 3387
from dataclasses import dataclass, asdict
from .electrode_model_template import ElectrodeModel
import netgen
import netgen.occ as occ
import numpy as np


@dataclass
class MedtronicParameters():
    # dimensions [mm]
    tip_length: float
    contact_length: float
    contact_spacing: float
    lead_diameter: float
    total_length: float


class MedtronicModel(ElectrodeModel):
    """Medtronic DBS electrode.

    Attributes
    ----------
    parameters : MedtronicParameters
        Parameters for the Medtronic geometry.

    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    position : tuple
        Position vector (x,y,z) of electrode tip.
    """

    _n_contacts = 4

    def parameter_check(self):
        # Check to ensure that all parameters are at least 0
        for param in (asdict(self._parameters).values()):
            if (param < 0):
                raise ValueError("Parameter values cannot be less than zero")

    def _construct_encapsulation_geometry(self, thickness: float) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of encapsulation layer around electrode.

        Parameters
        ----------
        thickness : float
            Thickness of encapsulation layer.

        Returns
        -------
        netgen.libngpy._NgOCC.TopoDS_Shape
        """
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
        # TODO check
        electrode = contacts
        electrode = netgen.occ.Glue([self.__body() - contacts, contacts])
        return electrode.Move(v=self._position)

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
        point = (0, 0, 0)
        radius = self._parameters.lead_diameter * 0.5
        height = self._parameters.contact_length
        contact = occ.Cylinder(p=point, d=self._direction, r=radius, h=height)
        distance = self._parameters.tip_length
        contacts = []
        for count in range(self._n_contacts):
            name = self._boundaries['Contact_{}'.format(count + 1)]
            contact.bc(name)
            min_edge_z_val = float("inf")
            max_edge_z_val = float("-inf")
            for edge in contact.edges:
                if (edge.center.z < min_edge_z_val):
                    min_edge_z_val = edge.center.z
                    min_edge = edge   
                if (edge.center.z > max_edge_z_val):
                    max_edge_z_val = edge.center.z
                    max_edge = edge    
            # Only name edge with the min and max z values (edge between the non-contact and contact surface)
            min_edge.name = name
            max_edge.name = name
            vector = tuple(np.array(self._direction) * distance)
            contacts.append(contact.Move(vector))
            distance += self._parameters.contact_length + self._parameters.contact_spacing

        return netgen.occ.Glue(contacts)


class MedtronicSenSightModel(ElectrodeModel):
    """Medtronic SenSight electrode.

    Attributes
    ----------
    parameters : MedtronicParameters
        Parameters for the Medtronic geometry.

    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    position : tuple
        Position vector (x,y,z) of electrode tip.
    """
    _n_contacts = 8

    def parameter_check(self):
        # Check to ensure that all parameters are at least 0
        for param in (asdict(self._parameters).values()):
            if (param < 0):
                raise ValueError("Parameter values cannot be less than zero")

    def _construct_encapsulation_geometry(self, thickness: float) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of encapsulation layer around electrode.

        Parameters
        ----------
        thickness : float
            Thickness of encapsulation layer.

        Returns
        -------
        netgen.libngpy._NgOCC.TopoDS_Shape
        """
        radius = self._parameters.lead_diameter * 0.5 + thickness
        center = tuple(np.array(self._direction) * self._parameters.lead_diameter * 0.5)
        height = self._parameters.total_length - self._parameters.tip_length
        tip = netgen.occ.Sphere(c=center, r=radius)
        lead = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)
        encapsulation = tip + lead
        encapsulation.bc('EncapsulationLayerSurface')
        encapsulation.mat('EncapsulationLayer')
        return encapsulation.Move(v=self._position) - self.geometry

    def _construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        contacts = self.__contacts()
        # TODO check
        electrode = occ.Glue([self.__body() - contacts, contacts])
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
        vectors = []
        distance = self._parameters.tip_length
        for index in range(0, 4):
            vectors.append(tuple(np.array(self._direction) * distance))
            distance += self._parameters.contact_length + self._parameters.contact_spacing

        point = (0, 0, 0)
        radius = self._parameters.lead_diameter * 0.5
        height = self._parameters.contact_length
        contact = occ.Cylinder(p=point, d=self._direction, r=radius, h=height)
        axis = occ.Axis(p=point, d=self._direction)

        contact_directed = self.__contact_directed()
        contacts = [contact.Move(vectors[0]),
                    contact_directed.Move(vectors[1]),
                    contact_directed.Rotate(axis, 120).Move(vectors[1]),
                    contact_directed.Rotate(axis, 240).Move(vectors[1]),
                    contact_directed.Move(vectors[2]),
                    contact_directed.Rotate(axis, 120).Move(vectors[2]),
                    contact_directed.Rotate(axis, 240).Move(vectors[2]),
                    contact.Move(vectors[3])
                    ]

        for index, contact in enumerate(contacts, 1):
            name = self._boundaries['Contact_{}'.format(index)]
            contact.bc(name)
            # Label max z value and min z value for contact_1 and contact_8
            if (name == 'Contact_1' or name == 'Contact_8'):
                min_edge_z_val = float("inf")
                for edge in contact.edges:
                    if (edge.center.z < min_edge_z_val):
                        min_edge_z_val = edge.center.z
                        min_edge = edge
                min_edge.name = name
                max_edge_z_val = float("-inf")
                for edge in contact.edges:
                    if (edge.center.z > max_edge_z_val):
                        max_edge_z_val = edge.center.z
                        max_edge = edge
                max_edge.name = name
            else:
                # Label all the named contacts appropriately
                for edge in contact.edges:
                    if (edge.name is not None):
                        edge.name = name
        return netgen.occ.Fuse(contacts)

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

        contact = body - eraser.Rotate(axis, angle) - eraser.Rotate(axis, -angle)
        # Centering contact to label edges
        contact = contact.Rotate(axis, angle)
        # Find  max z, min z, max x, and max y values and label min x and min y edge
        max_z_val = max_y_val = max_x_val = float("-inf")
        min_z_val = float("inf")
        for edge in contact.edges:
            if (edge.center.z > max_z_val):
                max_z_val = edge.center.z       
            if (edge.center.z < min_z_val):
                min_z_val = edge.center.z  
            if (edge.center.x > max_x_val):
                max_x_val = edge.center.x
                max_x_edge = edge
            if (edge.center.y > max_y_val):
                max_y_val = edge.center.y
                max_y_edge = edge
        max_x_edge.name = "max x"
        max_y_edge.name = "max y"
        # Label only the outer edges of the contact with min z and max z values
        for edge in contact.edges:
            if (np.isclose(edge.center.z, max_z_val)
                and not (np.isclose(edge.center.x, radius/2)
                or np.isclose(edge.center.y, radius/2))):
                edge.name = "max z"

            elif (np.isclose(edge.center.z, min_z_val)
                  and not (np.isclose(edge.center.x, radius/2)
                  or np.isclose(edge.center.y, radius/2))):
                edge.name = "min z"
        # Reseting position so that 0 deg lies in the middle of contact
        contact = contact.Rotate(axis, angle)
        # TODO: check that the starting axis of the contacts are correct according to the documentation

        return contact

    def __direction_2(self):
        x, y, z = self._direction

        if not x and not y:
            return (0, 1, 0)

        if not x and not z:
            return (0, 0, 1)

        if not y and not z:
            return (0, 1, 0)

        return (x, y, not z)
    