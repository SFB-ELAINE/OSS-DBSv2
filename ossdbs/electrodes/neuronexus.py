# Copyright 2025 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

from dataclasses import asdict, dataclass

import netgen
import netgen.occ as occ
import numpy as np

from .electrode_model_template import ElectrodeModel
from .utilities import get_highest_edge, get_lowest_edge


@dataclass
class NeuroNexusParameters:
    """Electrode geometry parameters."""

    shank_thickness = 50  # um
    max_width = 125  # um  # Width at the top of the electrode
    min_width = 33  # um   # Width at the first contact
    total_length = 1500  # um
    angle_tip = 30  # degrees

    # Contact parameters
    contact_diameter = 30  # um
    contact_spacing = 50  # um


    def get_center_first_contact(self) -> float:
        """Returns distance between electrode tip and center of first contact."""

    def get_distance_l1_l4(self) -> float:
        """Returns distance between first level contact and fourth level contact."""
        return -1.0


class NeuroNexusElectrodeModel(ElectrodeModel):
    """NeuroNexus electrode.

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

    # TODO
    # def parameter_check(self):
    #     """Check geometry parameters."""

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
        raise NotImplementedError()

    def _construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        # Define the position of the first contact
        first_contact_position = 50  # This is the tip length from point to first contact center
        last_contact_position = first_contact_position + (number_contacts - 1) * contact_spacing

        # Create a custom profile for the tapered electrode
        profile = occ.WorkPlane()

        # Start at the tip (origin)
        profile.MoveTo(0, 0)

        # Create the sharp tip that reaches the first contact in 50 microns
        tip_half_width = min_width / 2  # Half-width at first contact
        profile.LineTo(tip_half_width, first_contact_position)

        # Calculate the taper from first contact to last contact
        y1 = first_contact_position
        y2 = last_contact_position
        w1 = min_width / 2  # Half-width at first contact
        w2 = max_width / 2  # Half-width at last contact

        # Add the tapered section from first to last contact
        profile.LineTo(w2, y2)

        # Add the top section - straight from last contact to end
        profile.LineTo(w2, total_length)

        # Add the right side, going back down
        profile.LineTo(-w2, total_length)

        # Mirror the profile for the other side
        profile.LineTo(-w2, y2)
        profile.LineTo(-w1, y1)
        profile.LineTo(0, 0)  # Return to origin to close the shape

        # Create face from the profile
        electrode_face = profile.Face()

        # Extrude to create 3D body
        body_3D = electrode_face.Extrude(shank_thickness)

        # Create electrode contacts
        electrode_names = []
        contacts = []

        # Calculate contact positions
        contact_positions = []
        for i in range(number_contacts):
            # Distribute contacts evenly along the electrode
            y_pos = first_contact_position + i * contact_spacing
            contact_positions.append(y_pos)


        # Create circular contacts
        for i, y_pos in enumerate(contact_positions):
            # Create a circle at the position
            contact = occ.WorkPlane(occ.Axes((0, y_pos, shank_thickness), n=occ.Z))
            contact.Circle(contact_diameter / 2)
            contact_face = contact.Face()
            
            # Name the contact
            name = f"Contact_{i+1}"
            electrode_names.append(name)
            
            # Add material and boundary condition
            contact_face.mat(name)
            contact_face.bc(name)
            
            # Add to list of contacts
            contacts.append(contact_face)

        electrode = occ.Glue([body_3D] + contacts)
        return electrode.Move(v=self._position)
