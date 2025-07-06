# Copyright 2023, 2024 Konstantin Butenko, Shruthi Chakravarthy
# Copyright 2023, 2024 Jan Philipp Payonk, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

from dataclasses import asdict, dataclass

import netgen
import netgen.occ as occ
import numpy as np

from .electrode_model_template import ElectrodeModel
from .utilities import get_highest_edge, get_lowest_edge


@dataclass
class MicroProbesRodentElectrodeParameters:
    """Electrode geometry parameters."""

    # dimensions [mm]
    # exposed: The length of the exposed wire between tip and lead (if any)
    exposed_wire: float
    contact_radius: float
    lead_radius: float
    total_length: float
    wire_radius: float

    def get_center_first_contact(self) -> float:
        """Returns distance between electrode tip and center of first contact."""
        return 0.5 * self.contact_radius

    def get_distance_l1_l4(self) -> float:
        """Returns distance between first level contact and fourth level contact."""
        return -1.0

    @property
    def lead_diameter(self) -> float:
        """Lead diameter."""
        return 2.0 * self.lead_radius


class MicroProbesRodentElectrodeModel(ElectrodeModel):
    """MicroProbes Custom Rodent electrode.

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

    @property
    def wire_exists(self):
        """Check if parts of wire are exposed."""
        return self._parameters.exposed_wire != 0

    def parameter_check(self):
        """Check geometry parameters."""
        for param_name, param_value in asdict(self._parameters).items():
            if param_name != "exposed_wire" and param_value < 0:
                raise ValueError(f"Parameter {param_name} cannot be less than zero.")
            elif param_name == "exposed_wire":
                contact_radius = getattr(self._parameters, "contact_radius", None)
                if contact_radius is not None and param_value < -contact_radius:
                    raise ValueError(
                        f"Parameter {param_name} cannot be less than the negative of "
                        "the contact radius."
                    )

        # check that electrode is long enough
        if (
            self._parameters.total_length
            < self._parameters.contact_radius + self._parameters.exposed_wire
        ):
            raise ValueError(
                """Total length cannot be less
                   than the length of exposed wire and contact radius."""
            )
        # check that wire is thick enough
        if self._parameters.exposed_wire > 0:
            if np.isclose(self._parameters.wire_radius, 0):
                raise ValueError(
                    """If exposed wire length is greater than zero,
                    must specify wire radius to be greater than zero."""
                )
        # wire cannot be wider than contact
        if self._parameters.wire_radius > self._parameters.contact_radius:
            raise ValueError("Wire radius cannot be bigger than contact radius")

    _n_contacts = 1

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
        encap_tip_radius = self._parameters.contact_radius + thickness
        encap_tip_center = tuple(
            np.array(self._direction) * self._parameters.contact_radius
        )
        encap_tip = occ.Sphere(c=encap_tip_center, r=encap_tip_radius)

        # define half space at tip_center to construct a hemsiphere as the contact tip
        half_space = netgen.occ.HalfSpace(p=encap_tip_center, n=self._direction)
        encap_tip = encap_tip * half_space

        encap_lead_ht = self._parameters.total_length - (
            self._parameters.exposed_wire + self._parameters.contact_radius
        )
        encap_lead = occ.Cylinder(
            p=encap_tip_center,
            d=self._direction,
            r=encap_tip_radius,
            h=encap_lead_ht,
        )

        """
        Note: the following code was used to fillet the encapsulation layer

        # Find tip edge with the max z value for fillet
        max_edge_z_val = float("-inf")
        for edge in encap_tip.edges:
            if edge.center.z > max_edge_z_val:
                max_edge_z_val = edge.center.z
                fillet_tipE = edge
        fillet_tipE.name = "fillet_tipE"

        encap_lead_radius = self._parameters.lead_radius + thickness
        encap_lead_ht = self._parameters.total_length - (
            self._parameters.exposed_wire + self._parameters.contact_radius
        )
        encap_lead_start_pt = tuple(
            np.array(self._direction)
            * (self._parameters.exposed_wire + self._parameters.contact_radius)
        )
        encap_lead = occ.Cylinder(
            p=encap_lead_start_pt,
            d=self._direction,
            r=encap_lead_radius,
            h=encap_lead_ht,
        )

        # Find lead edge with the min z value for fillet
        min_edge_z_val = float("inf")
        for edge in encap_lead.edges:
            if edge.center.z < min_edge_z_val:
                min_edge_z_val = edge.center.z
                fillet_leadE = edge
        fillet_leadE.name = "fillet_leadE"

        if self.wire_exists:
            encap_wire_radius = self._parameters.wire_radius + thickness
            encap_wire_start_pt = encap_tip_center
            encap_wire_ht = self._parameters.exposed_wire
            encap_wire = occ.Cylinder(
                p=encap_wire_start_pt,
                d=self._direction,
                r=encap_wire_radius,
                h=encap_wire_ht,
            )
            encapsulation = encap_wire + encap_tip + encap_lead

            # Find wire edges with min and max z value for fillet
            max_edge_z_val = float("-inf")
            min_edge_z_val = float("inf")
            for edge in encap_wire.edges:
                if edge.center.z < min_edge_z_val:
                    min_edge_z_val = edge.center.z
                    fillet_wireE1 = edge

                if edge.center.z > max_edge_z_val:
                    max_edge_z_val = edge.center.z
                    fillet_wireE2 = edge
            fillet_wireE1.name = "fillet_wireE1"
            fillet_wireE2.name = "fillet_wireE2"

            # Only run MakeFillet if sharp edges are present
            # Command is very sensitive to input parameters,
            # may have to implement a check here
            if encap_wire_radius != encap_tip_radius:
                encapsulation = encapsulation.MakeFillet(
                    [fillet_wireE1], encap_wire_radius / 24
                )
                encapsulation = encapsulation.MakeFillet(
                    [fillet_tipE], encap_tip_radius / 24
                )

            if encap_wire_radius != encap_lead_radius:
                encapsulation = encapsulation.MakeFillet(
                    [fillet_wireE2], encap_wire_radius / 24
                )
                encapsulation = encapsulation.MakeFillet(
                    [fillet_leadE], encap_lead_radius / 24
                )
        else:
            encapsulation = encap_tip + encap_lead
            # if (encap_tip_radius != encap_lead_radius):
            #    encapsulation = encapsulation.MakeFillet([fillet_leadE],
                                                          encap_lead_radius / 50)
            # TODO: Issues with the following command
            # encapsulation = encapsulation.MakeFillet([fillet_tipE], 0.00001)
        """

        encapsulation = encap_tip + encap_lead
        encapsulation.bc("EncapsulationLayerSurface")
        encapsulation.mat("EncapsulationLayer")
        return encapsulation.Move(v=self._position) - self.geometry

    def _construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        contacts = self._contacts()
        electrode = netgen.occ.Glue([self.__body() - contacts, contacts])
        return electrode.Move(v=self._position)

    # Body is defined here to only include the lead
    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        direction = self._direction
        lead_radius = self._parameters.lead_radius
        lead_height = (
            self._parameters.total_length
            - max(self._parameters.exposed_wire, 0)
            - self._parameters.contact_radius
        )

        # If wire doesn't exist, start point will be the same as the tip center
        lead_start_pt = tuple(
            np.array(direction)
            * (self._parameters.exposed_wire + self._parameters.contact_radius)
        )
        body = occ.Cylinder(p=lead_start_pt, d=direction, r=lead_radius, h=lead_height)
        body.bc(self._boundaries["Body"])
        return body

    def _contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        origin = (0, 0, 0)
        direction = (0, 0, 1)
        contact_radius = self._parameters.contact_radius
        tip_center = tuple(np.array(direction) * self._parameters.contact_radius)
        tip = occ.Sphere(c=tip_center, r=contact_radius)
        # If exposed wire exists,
        # we include the wire and tip as part of the contact object
        if self.wire_exists:
            if self._parameters.exposed_wire > 0:
                # Standard case, exposed wire
                lead_start_pt = tuple(
                    np.array(direction)
                    * (self._parameters.exposed_wire + self._parameters.contact_radius)
                )
                half_space = netgen.occ.HalfSpace(
                    p=occ.gp_Pnt(*lead_start_pt), n=occ.gp_Vec(*direction)
                )
                wire_height = self._parameters.exposed_wire
                wire_start_pt = tip_center
                wire = occ.Cylinder(
                    p=wire_start_pt,
                    d=direction,
                    r=self._parameters.wire_radius,
                    h=wire_height,
                )
                contact = (tip * half_space) + wire
            else:
                # Negative exposed_wire, meaning part of the tip is covered
                covering_height = abs(self._parameters.exposed_wire)
                cover_start_pt = tuple(
                    np.array(direction)
                    * (self._parameters.contact_radius - covering_height)
                )
                # Convert to gp_Pnt and gp_Vec
                cover_start_pt_pnt = occ.gp_Pnt(*cover_start_pt)
                normal_vec = occ.gp_Vec(*np.array(direction))
                half_space = netgen.occ.HalfSpace(p=cover_start_pt_pnt, n=normal_vec)
                covered_tip = tip * half_space
                contact = covered_tip
        else:
            # No exposed wire, simple contact
            half_space = netgen.occ.HalfSpace(
                p=occ.gp_Pnt(*tip_center), n=occ.gp_Vec(*direction)
            )
            contact = tip * half_space

        contact.bc(self._boundaries["Contact_1"])

        for edge in contact.edges:
            edge.name = "Contact_1"

        if np.allclose(self._direction, direction):
            return contact

        # Rotate electrode to match orientation if required
        rotation = tuple(
            np.cross(direction, self._direction)
            / np.linalg.norm(np.cross(direction, self._direction))
        )
        angle = np.degrees(np.arccos(self._direction[2]))
        return contact.Rotate(occ.Axis(p=origin, d=rotation), angle)


@dataclass
class MicroProbesSNEX100Parameters:
    """Electrode geometry parameters."""

    # dimensions [mm]
    core_electrode_length: float
    core_electrode_diameter: float
    core_tubing_length: float
    core_tubing_diameter: float
    outer_electrode_length: float
    outer_electrode_diameter: float
    outer_tubing_diameter: float
    total_length: float

    def get_center_first_contact(self) -> float:
        """Returns distance between electrode tip and center of first contact."""
        return 0.5 * self.core_electrode_length

    def get_distance_l1_l4(self) -> float:
        """Returns distance between first level contact and fourth level contact."""
        return -1.0

    @property
    def lead_diameter(self) -> float:
        """Lead diameter, used outermost point of SNEX."""
        return self.outer_tubing_diameter


class MicroProbesSNEX100Model(ElectrodeModel):
    """MicroProbes SNEX 100 Concentric Bipolar electrode.

    Attributes
    ----------
    parameters : MicroProbesSNEX100Parameters
        Parameters for MicroProbesSNEX100 geometry.

    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    translation : tuple
        Translation vector (x,y,z) of electrode.
    """

    def parameter_check(self):
        """Check geometry parameters."""
        # Check to ensure that all parameters are at least 0
        for param in asdict(self._parameters).values():
            if param < 0:
                raise ValueError("Parameter values cannot be less than zero")
        if (
            self._parameters.total_length
            < self._parameters.core_electrode_length
            + self._parameters.core_tubing_length
            + self._parameters.outer_electrode_length
        ):
            raise ValueError(
                """Total length cannot be less than the sum of
                   the lengths of the core electrode, tubing, and outer electrode."""
            )
        if (
            self._parameters.core_tubing_diameter
            < self._parameters.core_electrode_diameter
        ):
            raise ValueError(
                "Core tubing diameter cannot be less than core electrode diameter"
            )
        if (
            self._parameters.outer_electrode_diameter
            < self._parameters.core_tubing_diameter
        ):
            raise ValueError(
                "Outer electrode diameter cannot be less than core tubing diameter"
            )
        if (
            self._parameters.outer_tubing_diameter
            < self._parameters.outer_electrode_diameter
        ):
            raise ValueError(
                "Outer tubing diameter cannot be less than outer electrode diameter"
            )

    _n_contacts = 2

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
        # Constructong core electrode
        direction = self._direction
        distance_1 = self._parameters.core_electrode_diameter * 0.5
        point_1 = tuple(np.array(direction) * distance_1)
        radius_1 = self._parameters.core_electrode_diameter * 0.5 + thickness
        part_0 = occ.Sphere(c=point_1, r=radius_1)
        height_1 = self._parameters.core_electrode_length - distance_1
        part_1 = occ.Cylinder(p=point_1, d=direction, r=radius_1, h=height_1)

        # Find max Z value for for edge between core electrode and core tubing
        # TODO check if this is a good idea to find this edge
        # TODO maybe define new object instead of passing added shapes
        max_CoreE = get_highest_edge(part_1 + part_0)
        # TODO why is naming the edge important?
        # max_CoreE.name = "fillet"

        # Constructing core tubing
        distance_2 = self._parameters.core_electrode_length
        point_2 = tuple(np.array(direction) * distance_2)
        radius_2 = self._parameters.core_tubing_diameter * 0.5 + thickness
        height_2 = self._parameters.core_tubing_length
        part_2 = occ.Cylinder(p=point_2, d=direction, r=radius_2, h=height_2)

        # Find min Z value for outer electrode rim
        # and max Z value for edge between outer tubing and outer electrode
        min_CoreTubeE = get_lowest_edge(part_2)
        max_CoreTubeE = get_highest_edge(part_2)
        # min_CoreTubeE.name = "fillet_edge"
        # max_CoreTubeE.name = "fillet_edge"

        # Constructing Outer Electrode
        distance_3 = distance_2 + self._parameters.core_tubing_length
        point_3 = tuple(np.array(direction) * distance_3)
        radius_3 = self._parameters.outer_electrode_diameter * 0.5 + thickness
        height_3 = self._parameters.outer_electrode_length
        part_3 = occ.Cylinder(p=point_3, d=direction, r=radius_3, h=height_3)

        min_OuterE = get_lowest_edge(part_3)
        max_OuterE = get_highest_edge(part_3)

        # min_OuterE.name = "fillet_edge"
        # max_OuterE.name = "fillet_edge"

        # Constructing Outer tubing
        distance_4 = distance_3 + self._parameters.outer_electrode_length
        point_4 = tuple(np.array(direction) * distance_4)
        radius_4 = self._parameters.outer_tubing_diameter * 0.5 + thickness
        height_4 = self._parameters.total_length - distance_4
        part_4 = occ.Cylinder(p=point_4, d=direction, r=radius_4, h=height_4)

        # Find min Z value for for edge between outer tubing rim
        min_OuterTubeE = get_lowest_edge(part_4)
        # min_OuterTubeE.name = "fillet_edge"

        encapsulation = part_0 + part_1 + part_2 + part_3 + part_4
        # Run MakeFillet on edges - command is very sensitive to input parameters
        # TODO check radius values
        # outer tubing
        encapsulation = encapsulation.MakeFillet([min_OuterTubeE], radius_4 / 12)
        # outer electrode
        encapsulation = encapsulation.MakeFillet([max_OuterE], radius_3 / 12)
        encapsulation = encapsulation.MakeFillet([min_OuterE], radius_3 / 8)

        # core tubing
        encapsulation = encapsulation.MakeFillet([max_CoreTubeE], radius_2 / 4)
        encapsulation = encapsulation.MakeFillet([min_CoreTubeE], radius_3 / 16)

        # core electrode
        encapsulation = encapsulation.MakeFillet(max_CoreE.edges, radius_1 / 12)
        encapsulation.bc("EncapsulationLayerSurface")
        encapsulation.mat("EncapsulationLayer")

        return encapsulation.Move(v=self._position) - self.geometry

    def _construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        electrode = occ.Glue([self.__body(), self._contacts()])
        return electrode.Move(v=self._position)

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        # Defining the core tubing
        # using the start point of the cylinder as the tip center
        direction = self._direction
        distance_1 = self._parameters.core_electrode_length
        point_1 = tuple(np.array(self._direction) * distance_1)
        radius_1 = self._parameters.core_tubing_diameter * 0.5
        height_1 = self._parameters.core_tubing_length
        body_pt1 = occ.Cylinder(p=point_1, d=direction, r=radius_1, h=height_1)
        # Defining the edge between the core tubing and the outer electrode (contact_2)
        max_edge_z = get_highest_edge(body_pt1)
        max_edge_z.name = self._boundaries["Contact_2"]

        # Defining the outer tubing
        distance_2 = (
            self._parameters.core_electrode_length
            + self._parameters.core_tubing_length
            + self._parameters.outer_electrode_length
        )
        point_2 = tuple(np.array(direction) * distance_2)
        radius_2 = self._parameters.outer_tubing_diameter * 0.5
        height_2 = self._parameters.total_length - distance_2
        body_pt2 = occ.Cylinder(p=point_2, d=direction, r=radius_2, h=height_2)
        body = occ.Fuse([body_pt1, body_pt2])
        body.bc(self._boundaries["Body"])
        return body

    def _contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        origin = (0, 0, 0)
        direction = (0, 0, 1)
        radius_1 = self._parameters.core_electrode_diameter * 0.5
        center = tuple(np.array(direction) * radius_1)
        # define half space at tip_center
        # to construct a hemsiphere as part of the contact tip
        half_space = netgen.occ.HalfSpace(p=center, n=direction)
        contact_tip = occ.Sphere(c=center, r=radius_1) * half_space
        height = self._parameters.core_electrode_length - radius_1
        contact = occ.Cylinder(p=center, d=direction, r=radius_1, h=height)
        contact_1 = contact_tip + contact

        distance = (
            self._parameters.core_electrode_length + self._parameters.core_tubing_length
        )
        point = tuple(np.array(direction) * distance)
        radius_2 = self._parameters.outer_electrode_diameter * 0.5
        height_2 = self._parameters.outer_electrode_length
        contact_2 = occ.Cylinder(p=point, d=direction, r=radius_2, h=height_2)

        contact_1.bc(self._boundaries["Contact_1"])
        contact_2.bc(self._boundaries["Contact_2"])
        # Find edge with max z value for contact_1
        max_edge_z = get_highest_edge(contact_1)
        # Only name edge with the maximum z value for contact_1
        # (represents the edge between the non-contact and contact surface)
        max_edge_z.name = self._boundaries["Contact_1"]

        # Find edge with max z value for contact_2
        max_edge_z = get_highest_edge(contact_2)
        max_edge_z.name = self._boundaries["Contact_2"]

        if np.allclose(self._direction, direction):
            return netgen.occ.Fuse([contact_1, contact_2])
        # rotate electrode to match orientation
        # e.g. from z-axis to y-axis
        rotation = tuple(
            np.cross(direction, self._direction)
            / np.linalg.norm(np.cross(direction, self._direction))
        )
        angle = np.degrees(np.arccos(self._direction[2]))
        return netgen.occ.Fuse([contact_1, contact_2]).Rotate(
            occ.Axis(p=origin, d=rotation), angle
        )

    def get_max_mesh_size_contacts(self, ratio: float) -> float:
        """Use electrode's contact size to estimate maximal mesh size.

        Parameters
        ----------
        ratio: float
            Ratio between characteristic contact size and maximal mesh size.

        """
        return self._parameters.core_electrode_diameter / ratio
