# Copyright 2023, 2024 Jan Philipp Payonk, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

# Dixi Microtechniques SEEG
from dataclasses import asdict, dataclass

import netgen
import netgen.occ as occ
import numpy as np

from .electrode_model_template import ElectrodeModel
from .utilities import get_highest_edge, get_lowest_edge


@dataclass
class DixiSEEGParameters:
    """Parameters for Dixi electrodes."""

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

    _n_contacts = None  # set to actual value for each electrode

    def __init__(
        self,
        parameters: DixiSEEGParameters,
        rotation: float = 0,
        direction: tuple = (0, 0, 1),
        position: tuple = (0, 0, 0),
    ) -> None:
        # set number of contacts
        self._n_contacts = parameters.n_contacts
        # set all other parameters
        super().__init__(parameters, rotation, direction, position)

    def parameter_check(self):
        """Check electrode parameters."""
        # ensure that all parameters are at least 0
        for param in asdict(self._parameters).values():
            if param < 0:
                raise ValueError("Parameter values cannot be less than zero")
        # check that number of contacts has been set correctly
        if not isinstance(self.n_contacts, int):
            raise ValueError(
                "The number of contacts has to be supplied as an integer value."
            )

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
        radius = self._parameters.lead_diameter * 0.5
        origin = (0, 0, 0)
        direction = (0, 0, 1)

        center = tuple(np.array(direction) * radius)
        # define half space at tip_center
        # to construct a hemisphere as part of the contact tip
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
            max_edge = get_highest_edge(contact)
            max_edge.name = name

            # first contact is different from others
            if count == 0:
                distance = (
                    self._parameters.contact_length + self._parameters.contact_spacing
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
