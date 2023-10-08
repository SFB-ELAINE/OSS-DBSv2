from dataclasses import asdict, dataclass

import netgen
import netgen.occ as occ
import numpy as np

from .electrode_model_template import ElectrodeModel


@dataclass
class MicroElectrodeParameters:
    # dimensions [mm]
    tip_length: float
    tip_diameter: float
    contact_length: float
    lead_diameter: float
    total_length: float

    def get_center_first_contact(self) -> float:
        """Returns distance between electrode tip and center of first contact."""
        return 0.5 * self.tip_length

    def get_distance_l1_l4(self) -> float:
        """Returns distance between first level contact and fourth level contact."""
        return -1.0


class MicroElectrodeModel(ElectrodeModel):
    """MicroElectrode.

    Attributes
    ----------
    parameters : MicroProbesRodentElectrodeParameters
        Parameters for MicroProbes Rodent Electrode geometry.

    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    position : tuple
        Position vector (x,y,z) of electrode tip.
    """

    _n_contacts = 1

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
        center = tuple(np.array(self._direction) * self._parameters.tip_length * 0.5)
        radius = self._parameters.tip_length + thickness
        height = self._parameters.total_length - self._parameters.tip_length * 0.5
        tip = occ.Sphere(c=center, r=radius)
        lead = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)
        encapsulation = tip + lead
        encapsulation.bc("EncapsulationLayerSurface")
        encapsulation.mat("EncapsulationLayer")
        return encapsulation.Move(v=self._position) - self.geometry

    def _construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        contact = self.__contacts()
        electrode = netgen.occ.Glue([self.__body(), contact])
        return electrode.Move(v=self._position)

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius_lead = self._parameters.lead_diameter * 0.5
        center = tuple(np.array(self._direction) * self._parameters.tip_length)
        height_lead = self._parameters.total_length - self._parameters.tip_length
        lead = occ.Cylinder(p=center, d=self._direction, r=radius_lead, h=height_lead)
        lead.bc(self._boundaries["Body"])
        return lead

    def __contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        point = (0, 0, 0)
        radius = self._parameters.tip_diameter * 0.5
        height = self._parameters.contact_length
        contact = occ.Cylinder(p=point, d=self._direction, r=radius, h=height)
        contact.bc(self._boundaries["Contact_1"])
        min_edge_z_val = float("inf")
        max_edge_z_val = float("-inf")
        for edge in contact.edges:
            if edge.center.z < min_edge_z_val:
                min_edge_z_val = edge.center.z
                min_edge = edge
            if edge.center.z > max_edge_z_val:
                max_edge_z_val = edge.center.z
                max_edge = edge
        min_edge.name = "Contact_1"
        contact = contact.MakeFillet(contact.edges["Contact_1"], 0.05)
        max_edge.name = "Contact_1"
        return contact
