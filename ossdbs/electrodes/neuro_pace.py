# Boston Scientific (Marlborough, Massachusetts, USA) vercise
from dataclasses import asdict, dataclass

import netgen
import netgen.occ as occ
import numpy as np

from .electrode_model_template import ElectrodeModel


@dataclass
class NeuroPaceParameters:
    # dimensions [mm]
    tip_length: float
    contact_length: float
    contact_spacing: float
    lead_diameter: float
    total_length: float

    def get_center_first_contact(self) -> float:
        """Returns distance between electrode tip and center of first contact."""
        return 0.5 * self.tip_length

    def get_distance_l1_l4(self) -> float:
        """Returns distance between first level contact and fourth level contacts."""
        return 3 * (self.contact_length + self.contact_spacing)


class NeuroPaceModel(ElectrodeModel):
    """Neuro Pace DL-344 electrode .

    Attributes
    ----------
    parameters : BostonScientificVerciseParameters
        Parameters for the Boston Scientific Vercise geometry.

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
        geometry : netgen.libngpy._NgOCC.TopoDS_Shape
        """
        radius = self._parameters.lead_diameter * 0.5 + thickness
        height = self._parameters.total_length - self._parameters.tip_length
        center = tuple(np.array(self._direction) * self._parameters.lead_diameter * 0.5)
        tip = occ.Sphere(c=center, r=radius)
        lead = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)
        encapsulation = tip + lead
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
        point = (0, 0, 0)
        radius = self._parameters.lead_diameter * 0.5
        height = self._parameters.contact_length
        contact = occ.Cylinder(p=point, d=self._direction, r=radius, h=height)
        contacts = []
        distance = self._parameters.tip_length
        for count in range(self._n_contacts):
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
            # Only name edge with the min and max z values (represents the edge between the non-contact and contact surface)
            min_edge.name = name
            max_edge.name = name
            vector = tuple(np.array(self._direction) * distance)
            contacts.append(contact.Move(vector))
            distance += (
                self._parameters.contact_length + self._parameters.contact_spacing
            )
        return netgen.occ.Glue(contacts)