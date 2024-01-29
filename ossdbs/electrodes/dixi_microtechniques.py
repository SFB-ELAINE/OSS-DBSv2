# Dixi Microtechniques SEEG
from dataclasses import asdict, dataclass

import netgen
import netgen.occ as occ
import numpy as np

from .electrode_model_template import ElectrodeModel


@dataclass
class DixiSEEGParameters:
    # dimensions [mm]
    tip_length: float
    contact_length: float
    contact_spacing: float
    lead_diameter: float
    total_length: float
    n_contacts: int

    def get_center_first_contact(self) -> float:
        """Returns distance between electrode tip and center of first contact."""
        return 0.5 * self.contact_length

    def get_distance_l1_l4(self) -> float:
        """Returns distance between first level contact and fourth level contact."""
        return 3 * (self.contact_length + self.contact_spacing)


class DixiSEEGModel(ElectrodeModel):
    """Dixi Microtechniques SEEG electrode.

    Attributes
    ----------
    parameters : DixiSEEGParameters
        Parameters for Dixi Microtechniques SEEG electrode geometry.

    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    position : tuple
        Position vector (x,y,z) of electrode tip.
    """

    _n_contacts = 18  # TODO set to actual value from each dataclass
    # def __init__(self, parameters: DixiSEEGParameters, rotation: float, direction: tuple, position: tuple):
    #    self._n_contacts = parameters.n_contacts

    def parameter_check(self):
        # Check to ensure that all parameters are at least 0
        for param in asdict(self._parameters).values():
            if param < 0:
                raise ValueError("Parameter values cannot be less than zero")

    def _construct_encapsulation_geometry(
        self, thickness: float
    ) -> netgen.libngpy._NgOCC.TopoDS_Shape:
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
        encapsulation.bc("EncapsulationLayerSurface")
        encapsulation.mat("EncapsulationLayer")
        return encapsulation.Move(v=self._position) - self.geometry

    def _construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        contacts = self.__contacts()
        # TODO check
        electrode = netgen.occ.Glue([self.__body() - contacts, contacts])
        return electrode.Move(v=self._position)

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self._parameters.lead_diameter * 0.5
        center = tuple(np.array(self._direction) * radius)
        tip = occ.Sphere(c=center, r=radius)
        height = self._parameters.total_length - self._parameters.tip_length
        lead = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)
        body = tip + lead
        body.bc(self._boundaries["Body"])
        return body

    def __contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self._parameters.lead_diameter * 0.5
        direction = self._direction

        center = tuple(np.array(direction) * radius)
        # define half space at tip_center to use to construct a hemsiphere as part of the contact tip
        half_space = netgen.occ.HalfSpace(p=center, n=direction)
        contact_tip = occ.Sphere(c=center, r=radius) * half_space
        h_pt2 = self._parameters.contact_length - radius
        contact_pt2 = occ.Cylinder(p=center, d=direction, r=radius, h=h_pt2)
        # defining first contact
        contact = contact_tip + contact_pt2
        height = self._parameters.contact_length
        contact_cyl = occ.Cylinder(p=(0, 0, 0), d=direction, r=radius, h=height)

        contacts = []
        for count in range(self._parameters.n_contacts):
            name = self._boundaries[f"Contact_{count + 1}"]
            contact.bc(name)
            min_edge_z_val = float("inf")
            max_edge_z_val = float("-inf")
            for edge in contact.edges:
                if edge.center.z < min_edge_z_val:
                    min_edge_z_val = edge.center.z
                    min_edge = edge
                if edge.center.z > max_edge_z_val:
                    max_edge_z_val = edge.center.z
                    max_edge = edge
                    # Only name edge with the max z value for contact_1
                max_edge.name = name
            # first contact is different from others
            if count == 0:
                distance = (
                    self._parameters.contact_length + self._parameters.contact_spacing
                )
                contacts.append(contact)
                contact = contact_cyl
            else:
                min_edge.name = name
                vector = tuple(np.array(self._direction) * distance)
                contacts.append(contact.Move(vector))
                distance += (
                    self._parameters.contact_length + self._parameters.contact_spacing
                )

        return occ.Glue(contacts)
