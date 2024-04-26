# Copyright 2023, 2024 Konstantin Butenko, Shruthi Chakravarthy
# Copyright 2023, 2024 Jan Philipp Payonk, Johannes Reding
# Copyright 2023, 2024 Tom Reincke, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

# Boston Scientific (Marlborough, Massachusetts, USA) vercise
from dataclasses import dataclass

import netgen
import netgen.occ as occ
import numpy as np

from .electrode_model_template import ElectrodeModel
from .utilities import get_highest_edge, get_lowest_edge


@dataclass
class BostonScientificVerciseDirectedParameters:
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
        """Returns distance between first level contact and fourth level contacts."""
        return 3 * (self.contact_length + self.contact_spacing)


class BostonScientificVerciseDirectedModel(ElectrodeModel):
    """Boston Scientific vercise directional lead electrode.

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
        netgen.libngpy._NgOCC.TopoDS_Shape
        """
        radius = self._parameters.lead_diameter * 0.5 + thickness
        center = tuple(np.array(self._direction) * self._parameters.lead_diameter * 0.5)
        height = self._parameters.total_length - self._parameters.tip_length
        tip = netgen.occ.Sphere(c=center, r=radius)
        lead = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)
        encapsulation = tip + lead
        encapsulation.bc("EncapsulationLayerSurface")
        encapsulation.mat("EncapsulationLayer")
        return encapsulation.Move(v=self._position) - self.geometry

    def _construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        contacts = self._contacts()
        # TODO check
        electrode = occ.Glue([self.__body() - contacts, contacts])
        axis = occ.Axis(p=(0, 0, 0), d=self._direction)
        rotated_electrode = electrode.Rotate(axis=axis, ang=self._rotation)
        return rotated_electrode.Move(v=self._position)

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self._parameters.lead_diameter * 0.5
        center = tuple(np.array(self._direction) * radius)
        height = self._parameters.total_length - self._parameters.tip_length
        body = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)
        body.bc(self._boundaries["Body"])
        return body

    def _contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self._parameters.lead_diameter * 0.5
        direction = self._direction
        center = tuple(np.array(direction) * radius)
        # define half space at tip_center
        # to construct a hemisphere as part of the contact tip
        half_space = netgen.occ.HalfSpace(p=center, n=direction)
        contact_tip = occ.Sphere(c=center, r=radius) * half_space
        h_pt2 = self._parameters.tip_length - radius
        contact_pt2 = occ.Cylinder(p=center, d=direction, r=radius, h=h_pt2)
        contact_1 = contact_tip + contact_pt2

        vectors = []
        distance = self._parameters.tip_length + self._parameters.contact_spacing
        for _ in range(0, 3):
            vectors.append(tuple(np.array(self._direction) * distance))
            distance += (
                self._parameters.contact_length + self._parameters.contact_spacing
            )

        point = (0, 0, 0)
        height = self._parameters.contact_length
        axis = occ.Axis(p=point, d=self._direction)
        contact_8 = occ.Cylinder(p=point, d=self._direction, r=radius, h=height)
        contact_directed = self._contact_directed()

        contacts = [
            contact_1,
            contact_directed.Move(v=vectors[0]),
            contact_directed.Rotate(axis, 120).Move(v=vectors[0]),
            contact_directed.Rotate(axis, 240).Move(v=vectors[0]),
            contact_directed.Move(v=vectors[1]),
            contact_directed.Rotate(axis, 120).Move(v=vectors[1]),
            contact_directed.Rotate(axis, 240).Move(v=vectors[1]),
            contact_8.Move(v=vectors[2]),
        ]

        for index, contact in enumerate(contacts, 1):
            name = self._boundaries[f"Contact_{index}"]
            contact.bc(name)
            # Label max z value and min z value for contact_8
            if name == "Contact_8":
                min_edge = get_lowest_edge(contact)
                min_edge.name = name
            # Only label contact edge with maximum z value for contact_1
            if name == "Contact_1" or name == "Contact_8":
                max_edge = get_highest_edge(contact)
                max_edge.name = name
            else:
                # Label all the named contacts appropriately
                for edge in contact.edges:
                    if edge.name is not None:
                        edge.name = name
        return netgen.occ.Fuse(contacts)

    def _contact_directed(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        point = (0, 0, 0)
        radius = self._parameters.lead_diameter * 0.5
        height = self._parameters.contact_length
        body = occ.Cylinder(p=point, d=self._direction, r=radius, h=height)
        # tilted y-vector marker is in YZ-plane and orthogonal to _direction
        new_direction = (0, self._direction[2], -self._direction[1])
        eraser = occ.HalfSpace(p=point, n=new_direction)
        delta = 15
        angle = 30 + delta
        axis = occ.Axis(p=point, d=self._direction)

        contact = body - eraser.Rotate(axis, angle) - eraser.Rotate(axis, -angle)
        # Centering contact to label edges
        contact = contact.Rotate(axis, angle)
        # TODO refactor / wrap in function
        # Find  max z, min z, max x, and max y values and label min x and min y edge
        max_z_val = max_y_val = max_x_val = float("-inf")
        min_z_val = float("inf")
        for edge in contact.edges:
            if edge.center.z > max_z_val:
                max_z_val = edge.center.z
            if edge.center.z < min_z_val:
                min_z_val = edge.center.z
            if edge.center.x > max_x_val:
                max_x_val = edge.center.x
                max_x_edge = edge
            if edge.center.y > max_y_val:
                max_y_val = edge.center.y
                max_y_edge = edge
        max_x_edge.name = "max x"
        max_y_edge.name = "max y"
        # Label only the outer edges of the contact with min z and max z values
        for edge in contact.edges:
            if np.isclose(edge.center.z, max_z_val) and not (
                np.isclose(edge.center.x, radius / 2)
                or np.isclose(edge.center.y, radius / 2)
            ):
                edge.name = "max z"
            elif np.isclose(edge.center.z, min_z_val) and not (
                np.isclose(edge.center.x, radius / 2)
                or np.isclose(edge.center.y, radius / 2)
            ):
                edge.name = "min z"

        # TODO check that the starting axis of the contacts
        # are correct according to the documentation
        contact = contact.Rotate(axis, -angle)

        return contact


@dataclass
class BostonScientificVerciseParameters:
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
        """Returns distance between first level contact and fourth level contacts."""
        return 7 * (self.contact_length + self.contact_spacing)


class BostonScientificVerciseModel(ElectrodeModel):
    """Boston Scientific vercise standard lead electrode.

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
        contacts = self._contacts()
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

    def _contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        point = (0, 0, 0)
        radius = self._parameters.lead_diameter * 0.5
        height = self._parameters.contact_length
        contact = occ.Cylinder(p=point, d=self._direction, r=radius, h=height)
        contacts = []
        distance = self._parameters.tip_length
        for count in range(self._n_contacts):
            name = self._boundaries[f"Contact_{count + 1}"]
            contact.bc(name)
            min_edge = get_lowest_edge(contact)
            max_edge = get_highest_edge(contact)
            # Only name edge with the min and max z values
            # (represents the edge between the non-contact and contact surface)
            min_edge.name = name
            max_edge.name = name
            vector = tuple(np.array(self._direction) * distance)
            contacts.append(contact.Move(vector))
            distance += (
                self._parameters.contact_length + self._parameters.contact_spacing
            )
        return netgen.occ.Glue(contacts)
