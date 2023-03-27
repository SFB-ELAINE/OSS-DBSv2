# PINS Medical L301
from abc import ABC
from dataclasses import dataclass
from .electrode_model import ElectrodeModel
import netgen
import netgen.occ as occ
import numpy as np


@dataclass
class PINSMedicalParameters(ABC):
    # dimensions [mm]
    tip_length: float
    contact_length: float
    contact_spacing: float
    lead_diameter: float
    total_length: float


class PINSMedicalModel(ElectrodeModel):
    """PINS Medical L301 electrode.

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
                 parameters: PINSMedicalParameters,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
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
        center = tuple(np.array(self._direction) * radius)
        tip = occ.Sphere(c=center, r=radius)
        height = self._total_length - self._tip_length
        lead = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)
        body = tip + lead
        body.bc(self._boundaries['Body'])
        return body

    def __contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        point = (0, 0, 0)
        radius = self._lead_diameter * 0.5
        height = self._contact_length
        contact = occ.Cylinder(p=point, d=self._direction, r=radius, h=height)

        distance_1 = self._tip_length
        distance_2 = distance_1 + self._contact_length + self._contact_spacing
        distance_3 = distance_2 + self._contact_length + self._contact_spacing
        distance_4 = distance_3 + self._contact_length + self._contact_spacing

        vector_1 = tuple(np.array(self._direction) * distance_1)
        vector_2 = tuple(np.array(self._direction) * distance_2)
        vector_3 = tuple(np.array(self._direction) * distance_3)
        vector_4 = tuple(np.array(self._direction) * distance_4)

        contacts = [contact.Move(v=vector_1),
                    contact.Move(v=vector_2),
                    contact.Move(v=vector_3),
                    contact.Move(v=vector_4)]

        for index, contact in enumerate(contacts, 1):
            name = self._boundaries['Contact_{}'.format(index)]
            contact.bc(name)
            for edge in contact.edges:
                edge.name = name

        return netgen.occ.Glue(contacts)
