# Copyright 2023, 2024 Konstantin Butenko, Shruthi Chakravarthy
# Copyright 2023, 2024 Jan Philipp Payonk, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

# SceneRay 1242
from dataclasses import dataclass

import netgen
import netgen.occ as occ
import numpy as np

from .electrode_model_template import ElectrodeModel
from .utilities import get_highest_edge, get_lowest_edge


@dataclass
class SceneRay1242Parameters:
    """Electrode geometry parameters."""

    # dimensions [mm]
    tip_length: float
    contact_length: float
    first_contact_spacing: float
    contact_spacing: float
    lead_diameter: float
    total_length: float

    def get_center_first_contact(self) -> float:
        """Returns distance between electrode tip and center of first contact."""
        return self.tip_length + 0.5 * self.contact_length

    def get_distance_l1_l4(self) -> float:
        """Returns distance between first level contact and fourth level contact."""
        return (
            2 * (self.contact_length + self.contact_spacing)
            + self.first_contact_spacing
            + self.contact_length
        )


class SceneRay1242Model(ElectrodeModel):
    """SceneRay DBS electrode.

    Attributes
    ----------
    parameters : SceneRayParameters
        Parameters for the SceneRay geometry.

    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    position : tuple
        Position vector (x,y,z) of electrode tip.
    """

    _n_contacts = 4

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
        contacts = self._contacts()
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
        body.bc(self._boundaries["Body"])
        return body

    def _contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        origin = (0, 0, 0)
        direction = (0, 0, 1)
        radius = self._parameters.lead_diameter * 0.5
        height = self._parameters.contact_length
        contact = occ.Cylinder(p=origin, d=direction, r=radius, h=height)
        distance = self._parameters.tip_length
        contacts = []
        for count in range(self._n_contacts):
            name = self._boundaries[f"Contact_{count + 1}"]
            contact.bc(name)
            min_edge = get_lowest_edge(contact)
            max_edge = get_highest_edge(contact)
            # Only name edge with the min and max z values
            # (edge between the non-contact and contact surface)
            min_edge.name = name
            max_edge.name = name
            vector = tuple(np.array(direction) * distance)
            contacts.append(contact.Move(vector))

            if count == 0:
                # SceneRay1242 has unique spacing between the distal pair of contacts
                distance += (
                    self._parameters.contact_length
                    + self._parameters.first_contact_spacing
                )
            else:
                distance += (
                    self._parameters.contact_length + self._parameters.contact_spacing
                )

        if np.allclose(self._direction, direction):
            return netgen.occ.Fuse(contacts)
        # rotate electrode to match orientation
        # e.g. from z-axis to y-axis
        rotation = tuple(
            np.cross(direction, self._direction)
            / np.linalg.norm(np.cross(direction, self._direction))
        )
        angle = np.degrees(np.arccos(self._direction[2]))
        return netgen.occ.Fuse(contacts).Rotate(occ.Axis(p=origin, d=rotation), angle)
