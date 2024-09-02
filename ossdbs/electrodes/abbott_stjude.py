# Copyright 2023, 2024 Konstantin Butenko, Shruthi Chakravarthy
# Copyright 2023, 2024 Jan Philipp Payonk, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

# Abbott/St Jude Active Tip 6142-6145
import logging
from dataclasses import dataclass

import netgen
import netgen.occ as occ
import numpy as np

from .electrode_model_template import ElectrodeModel
from .utilities import get_electrode_spin_angle, get_highest_edge, get_lowest_edge

_logger = logging.getLogger(__name__)


@dataclass
class AbbottStJudeActiveTipParameters:
    """Electrode geometry parameters."""

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
        """Returns distance between first level contact and fourth level contact."""
        return 3 * (self.contact_length + self.contact_spacing)


class AbbottStJudeActiveTipModel(ElectrodeModel):
    """Abbott/St Jude Active Tip electrode.

    Attributes
    ----------
    parameters : AbbottStJudeActiveTipParameters
        Parameters for Abbott Stjude Geometry.

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
        electrode = occ.Glue([self.__body() - contacts, contacts])
        return electrode.Move(v=self._position)

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self._parameters.lead_diameter * 0.5
        center = tuple(np.array(self._direction) * radius)
        height = self._parameters.total_length - self._parameters.tip_length
        body = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)
        body.bc(self._boundaries["Body"])
        return body

    def _contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self._parameters.lead_diameter * 0.5

        origin = (0, 0, 0)
        direction = (0, 0, 1)

        center = tuple(np.array(direction) * radius)
        # define half space at tip_center
        # to construct a hemisphere as part of the contact tip
        half_space = netgen.occ.HalfSpace(p=center, n=direction)
        contact_tip = occ.Sphere(c=center, r=radius) * half_space
        h_pt2 = self._parameters.tip_length - radius
        contact_pt2 = occ.Cylinder(p=center, d=direction, r=radius, h=h_pt2)
        # defining first contact
        contact = contact_tip + contact_pt2
        height = self._parameters.contact_length
        contact_cyl = occ.Cylinder(p=(0, 0, 0), d=direction, r=radius, h=height)

        contacts = []
        for count in range(self._n_contacts):
            name = self._boundaries[f"Contact_{count + 1}"]
            contact.bc(name)

            max_edge = get_highest_edge(contact)
            max_edge.name = name

            # first contact is different from others
            if count == 0:
                distance = (
                    self._parameters.tip_length + self._parameters.contact_spacing
                )
                contacts.append(contact)
                contact = contact_cyl
            else:
                min_edge = get_lowest_edge(contact)
                min_edge.name = name
                vector = tuple(np.array(direction) * distance)
                contacts.append(contact.Move(vector))
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


@dataclass
class AbbottStJudeParameters:
    """Electrode geometry parameters."""

    # dimensions [mm]
    tip_length: float
    contact_length: float
    contact_spacing: float
    lead_diameter: float
    total_length: float

    def get_center_first_contact(self) -> float:
        """Returns distance between electrode tip and center of first contact."""
        return self.tip_length + 0.5 * self.contact_length

    def get_distance_l1_l4(self) -> float:
        """Returns distance between first level contact and fourth level contact."""
        return 3 * (self.contact_length + self.contact_spacing)


class AbbottStJudeDirectedModel(ElectrodeModel):
    """Abbott/St Jude Directed electrode.

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

    _n_contacts = 8

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
        center = tuple(np.array(self._direction) * self._parameters.lead_diameter * 0.5)
        radius = self._parameters.lead_diameter * 0.5 + thickness
        height = self._parameters.total_length - self._parameters.tip_length
        tip = occ.Sphere(c=center, r=radius)
        lead = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)
        encapsulation = tip + lead
        encapsulation.mat("EncapsulationLayer")
        return encapsulation.Move(v=self._position) - self.geometry

    def _construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        contacts = self._contacts()
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
        body.bc(self._boundaries["Body"])
        return body

    def _contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        vectors = []
        direction = (0, 0, 1)
        distance = self._parameters.tip_length
        for _ in range(0, 4):
            vectors.append(tuple(np.array(direction) * distance))
            distance += (
                self._parameters.contact_length + self._parameters.contact_spacing
            )

        origin = (0, 0, 0)
        radius = self._parameters.lead_diameter * 0.5
        height = self._parameters.contact_length
        contact = occ.Cylinder(p=origin, d=direction, r=radius, h=height)
        axis = occ.Axis(p=origin, d=direction)

        contact_directed = self._contact_directed()
        contacts = [
            contact.Move(vectors[0]),
            contact_directed.Move(vectors[1]),
            contact_directed.Rotate(axis, 240).Move(vectors[1]),
            contact_directed.Rotate(axis, 120).Move(vectors[1]),
            contact_directed.Move(vectors[2]),
            contact_directed.Rotate(axis, 240).Move(vectors[2]),
            contact_directed.Rotate(axis, 120).Move(vectors[2]),
            contact.Move(vectors[3]),
        ]

        for index, contact in enumerate(contacts, 1):
            name = self._boundaries[f"Contact_{index}"]
            contact.bc(name)
            # Label max z value and min z value for contact_1 and contact_8
            if name == "Contact_1" or name == "Contact_8":
                min_edge = get_lowest_edge(contact)
                min_edge.name = name

                max_edge = get_highest_edge(contact)
                max_edge.name = name
            else:
                # Label all the named contacts appropriately
                for edge in contact.edges:
                    if edge.name == "Rename":
                        edge.name = name

        if np.allclose(self._direction, direction):
            return netgen.occ.Fuse(contacts)
        else:
            rotation = tuple(
                np.cross(direction, self._direction)
                / np.linalg.norm(np.cross(direction, self._direction))
            )
            angle = np.degrees(np.arccos(self._direction[2]))
            rotated_geo = netgen.occ.Fuse(contacts).Rotate(
                occ.Axis(p=origin, d=rotation), angle
            )
            rotation_angle = get_electrode_spin_angle(rotation, angle, self._direction)
            if np.isclose(rotation_angle, 0):
                return rotated_geo
            return rotated_geo.Rotate(
                occ.Axis(p=(0, 0, 0), d=self._direction), rotation_angle
            )

    def _contact_directed(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        origin = (0, 0, 0)
        direction = (0, 0, 1)
        radius = self._parameters.lead_diameter * 0.5
        height = self._parameters.contact_length
        body = occ.Cylinder(p=origin, d=direction, r=radius, h=height)
        # tilted y-vector marker is in YZ-plane and orthogonal to _direction
        new_direction = (0, 1, 0)
        eraser = occ.HalfSpace(p=origin, n=new_direction)
        angle = 45
        axis = occ.Axis(p=origin, d=direction)

        contact = body - eraser.Rotate(axis, angle) - eraser.Rotate(axis, -angle)

        # Label all outer edges
        for edge in contact.edges:
            edge_center = np.array([edge.center.x, edge.center.y, edge.center.z])

            # Skip center edge
            if np.allclose(np.cross(edge_center, direction), 0):
                continue

            new_center = np.dot(edge_center, direction) * np.array(direction)

            # Mark only outer edges
            if not np.isclose(np.linalg.norm(edge_center - new_center), radius / 2):
                edge.name = "Rename"

        return contact
