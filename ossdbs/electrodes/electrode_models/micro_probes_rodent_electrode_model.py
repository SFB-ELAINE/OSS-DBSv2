
from dataclasses import dataclass
from .electrode_model import ElectrodeModel
import netgen
import netgen.occ as occ
import numpy as np


@dataclass
class MicroProbesRodentElectrodeParameters():
    # dimensions [mm]
    tube_thickness: float
    contact_length: float
    lead_diameter: float
    total_length: float


class MicroProbesRodentElectrodeModel(ElectrodeModel):
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

    def __init__(self,
                 parameters: MicroProbesRodentElectrodeParameters,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        self._position = tuple(position)
        norm = np.linalg.norm(direction)
        self._direction = tuple(direction / norm) if norm else (0, 0, 1)
        self._boundaries = {'Body': 'Body', 'Contact_1': 'Contact_1'}
        self._contact_length = parameters.contact_length
        self._lead_diameter = parameters.lead_diameter
        self._total_length = parameters.total_length
        self._tube_thickness = parameters.tube_thickness

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
        height = self._total_length - self._lead_diameter * 0.5
        lead = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)
        tip = occ.Sphere(c=center, r=radius + thickness,)
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
        electrode = netgen.occ.Glue([self.__body() - contacts, contacts])
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
        direction = self._direction
        tube_radius = radius + self._tube_thickness
        center = tuple(np.array(self._direction) * radius)
        tip = occ.Sphere(c=center, r=tube_radius)
        height = self._total_length - radius
        lead = occ.Cylinder(p=center, d=direction, r=tube_radius, h=height)
        point = tuple(np.array(direction) * self._contact_length)
        space = occ.HalfSpace(p=point, n=direction)
        body = tip + lead - space
        body.bc(self._boundaries['Body'])
        return body

    def __contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self._lead_diameter * 0.5
        point = tuple(np.array(self._direction) * self._contact_length)
        space = netgen.occ.HalfSpace(p=point, n=self._direction)
        center = tuple(np.array(self._direction) * radius)
        tip = occ.Sphere(c=center, r=radius) * space
        height = self._contact_length - radius
        lead = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)

        if self._contact_length <= radius:
            contact = tip
        else:
            contact = tip + lead

        contact.bc(self._boundaries['Contact_1'])
        for edge in contact.edges:
            edge.name = self._boundaries['Contact_1']
        return contact
