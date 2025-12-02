# Copyright 2025 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

from dataclasses import dataclass

import netgen
import netgen.occ as occ
import numpy as np

from .electrode_model_template import ElectrodeModel


@dataclass
class NeuroNexusParameters:
    """Electrode geometry parameters."""

    # Shank geometry
    shank_thickness: float
    max_width: float  # Width at the top of the electrode
    min_width: float  # Width at the first contact
    total_length: float
    angle_tip: float
    tip_length: float  # distance tip to center of first contact

    # Contact parameters
    contact_diameter: float
    contact_spacing: float

    @property
    def lead_diameter(self) -> float:
        """Instead of lead diameter we use the widest part."""
        return max(self.shank_thickness, self.max_width)

    def get_center_first_contact(self) -> float:
        """Returns distance between electrode tip and center of first contact."""
        return self.tip_length

    def get_distance_l1_l4(self) -> float:
        """Returns distance between first level contact and fourth level contact."""
        return -1.0


class NeuroNexusElectrodeModel(ElectrodeModel):
    """NeuroNexus electrode.

    Attributes
    ----------
    parameters : NeuroNexusParameters
        Parameters for NeuroNexus electrode geometry.

    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    position : tuple
        Position vector (x,y,z) of electrode tip.
    """

    _n_contacts = 8

    def parameter_check(self):
        """Check geometry parameters."""
        if not np.isclose(self._parameters.angle_tip, 30.0):
            raise NotImplementedError("So far, the tip angle is hard-coded")
        if not np.isclose(self._parameters.tip_length, 50.0e-3):
            raise NotImplementedError("So far, the tip length is hard-coded")

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
        raise NotImplementedError()

    def _construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Construct electrode geometry."""
        # Define the position of the first contact
        first_contact_position = self._parameters.get_center_first_contact()
        last_contact_position = (
            first_contact_position
            + (self._n_contacts - 1) * self._parameters.contact_spacing
        )

        # Create a custom profile for the tapered electrode
        profile = occ.WorkPlane()

        # Start at the tip (origin)
        profile.MoveTo(0, 0)

        # Create the sharp tip that reaches the first contact in 50 microns
        tip_half_width = self._parameters.min_width / 2  # Half-width at first contact
        profile.LineTo(tip_half_width, first_contact_position)

        # Calculate the taper from first contact to last contact
        y1 = first_contact_position
        y2 = last_contact_position
        w1 = self._parameters.min_width / 2  # Half-width at first contact
        w2 = self._parameters.max_width / 2  # Half-width at last contact

        # Add the tapered section from first to last contact
        profile.LineTo(w2, y2)

        # Add the top section - straight from last contact to end
        profile.LineTo(w2, self._parameters.total_length)

        # Add the right side, going back down
        profile.LineTo(-w2, self._parameters.total_length)

        # Mirror the profile for the other side
        profile.LineTo(-w2, y2)
        profile.LineTo(-w1, y1)
        profile.LineTo(0, 0)  # Return to origin to close the shape

        # Create face from the profile
        electrode_face = profile.Face()

        # Extrude to create 3D body
        body_3D = electrode_face.Extrude(self._parameters.shank_thickness)
        # Name passive part
        body_3D.mat("Body")
        body_3D.bc("Body")

        # Create electrode contacts
        contacts = []

        # Calculate contact positions
        contact_positions = []
        for i in range(self._n_contacts):
            # Distribute contacts evenly along the electrode
            y_pos = first_contact_position + i * self._parameters.contact_spacing
            contact_positions.append(y_pos)

        # Create circular contacts
        for i, y_pos in enumerate(contact_positions):
            # Create a circle at the position
            # align contacts with the impl. coordinate
            contact = occ.WorkPlane(occ.Axes((0, y_pos, 0.0), n=occ.Z))
            contact.Circle(self._parameters.contact_diameter / 2)
            contact_face = contact.Face()
            # Name the contact
            name = f"Contact_{i + 1}"
            # Add material and boundary condition
            contact_face.mat(name)
            contact_face.bc(name)
            for edge in contact_face.edges:
                edge.name = name
            # Add to list of contacts
            contacts.append(contact_face)

        electrode = occ.Glue([body_3D, *contacts])
        # Rotate electrode so that it points in z and faces x
        origin = (0, 0, 0)
        x_direction = (1, 0, 0)
        z_direction = (0, 0, 1)
        angle = 90
        electrode = electrode.Rotate(occ.Axis(p=origin, d=x_direction), angle)
        electrode = electrode.Rotate(occ.Axis(p=origin, d=z_direction), -angle)

        # to prevent crashes
        occ_electrode = occ.OCCGeometry(electrode)
        occ_electrode.Heal()
        electrode = occ_electrode.shape
        # some face names get lost, rename
        contact_areas = 0
        for face in electrode.faces:
            if face.name is None:
                face.name = "Body"
            elif "Contact" in face.name:
                contact_areas += face.mass
        # check if contact areas are correct
        if not np.isclose(
            contact_areas,
            self._n_contacts * np.pi * (self._parameters.contact_diameter * 0.5) ** 2,
        ):
            raise ValueError("NeuroNexus electrodes was not correctly built.")

        return electrode.Move(v=self._position)
