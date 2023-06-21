
from dataclasses import dataclass
from .electrode_model_template import ElectrodeModel
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

    def _construct_encapsulation_geometry(self, thickness: float) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:
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
        height = self._parameters.total_length - self._parameters.lead_diameter * 0.5
        lead = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)
        tip = occ.Sphere(c=center, r=radius)
        encapsulation = tip + lead
        encapsulation.bc('EncapsulationLayerSurface')
        encapsulation.mat('EncapsulationLayer')
        return encapsulation.Move(v=self._position) - self.geometry

    def _construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        contacts = self.__contacts()
        # TODO check
        electrode = netgen.occ.Glue([self.__body() - contacts, contacts])
        return electrode.Move(v=self._position)

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self._parameters.lead_diameter * 0.5
        direction = self._direction
        tube_radius = radius + self._parameters.tube_thickness
        center = tuple(np.array(self._direction) * radius)
        tip = occ.Sphere(c=center, r=tube_radius)
        height = self._parameters.total_length - radius
        lead = occ.Cylinder(p=center, d=direction, r=tube_radius, h=height)
        point = tuple(np.array(direction) * self._parameters.contact_length)
        space = occ.HalfSpace(p=point, n=direction)
        body = tip + lead - space
        body.bc(self._boundaries['Body'])
        return body

    def __contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self._parameters.lead_diameter * 0.5
        point = tuple(np.array(self._direction) * self._parameters.contact_length)
        space = netgen.occ.HalfSpace(p=point, n=self._direction)
        center = tuple(np.array(self._direction) * radius)
        tip = occ.Sphere(c=center, r=radius) * space
        height = self._parameters.contact_length - radius
        lead = occ.Cylinder(p=center, d=self._direction, r=radius, h=height)

        # TODO check
        if self._parameters.contact_length <= radius:
            contact = tip
        else:
            contact = tip + lead

        contact.bc(self._boundaries['Contact_1'])
        for edge in contact.edges:
            edge.name = self._boundaries['Contact_1']
        return contact

@dataclass
class MicroProbesSNEX100Parameters():
    # dimensions [mm]
    core_electrode_length: float
    core_electrode_diameter: float
    core_tubing_length: float
    core_tubing_diameter: float
    outer_electrode_length: float
    outer_electrode_diameter: float
    outer_tubing_diameter: float
    total_length: float


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

    _n_contacts = 2

    def _construct_encapsulation_geometry(self, thickness: float) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of encapsulation layer around electrode.

        Parameters
        ----------
        thickness : float
            Thickness of encapsulation layer.

        Returns
        -------
        netgen.libngpy._NgOCC.TopoDS_Shape
        """
        direction = self._direction
        distance_1 = self._parameters.core_electrode_diameter * 0.5
        point_1 = tuple(np.array(self._direction) * distance_1)
        radius_1 = self._parameters.core_electrode_diameter * 0.5 + thickness
        part_0 = occ.Sphere(c=point_1, r=radius_1)
        height_1 = self._parameters.core_electrode_length - distance_1
        part_1 = occ.Cylinder(p=point_1, d=direction, r=radius_1, h=height_1)

        distance_2 = self._parameters.core_electrode_length - thickness
        point_2 = tuple(np.array(self._direction) * distance_2)
        radius_2 = self._parameters.core_tubing_diameter * 0.5 + thickness
        height_2 = self._parameters.core_tubing_length
        part_2 = occ.Cylinder(p=point_2, d=direction, r=radius_2, h=height_2)

        distance_3 = distance_2 + self._parameters.core_tubing_length
        point_3 = tuple(np.array(self._direction) * distance_3)
        radius_3 = self._parameters.outer_electrode_diameter * 0.5 + thickness
        height_3 = self._parameters.outer_electrode_length
        part_3 = occ.Cylinder(p=point_3, d=direction, r=radius_3, h=height_3)

        distance_4 = distance_3 + self._parameters.outer_electrode_length
        point_4 = tuple(np.array(self._direction) * distance_4)
        radius_4 = self._parameters.outer_tubing_diameter * 0.5 + thickness
        height_4 = self._parameters.total_length - distance_4
        part_4 = occ.Cylinder(p=point_4, d=direction, r=radius_4, h=height_4)

        encapsulation = part_0 + part_1 + part_2 + part_3 + part_4
        encapsulation.mat('EncapsulationLayer')
        return encapsulation.Move(v=self._position) - self.geometry

    def _construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        electrode = occ.Glue([self.__body(), self.__contacts()])
        return electrode.Move(v=self._position)

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        direction = self._direction
        distance_1 = self._parameters.core_electrode_length
        point_1 = tuple(np.array(self._direction) * distance_1)
        radius_1 = self._parameters.core_tubing_diameter * 0.5
        height_1 = self._parameters.core_tubing_length
        body_pt1 = occ.Cylinder(p=point_1, d=direction, r=radius_1, h=height_1)

        distance_2 = (self._parameters.core_electrode_length
                      + self._parameters.core_tubing_length
                      + self._parameters.outer_electrode_length)
        point_2 = tuple(np.array(direction) * distance_2)
        radius_2 = self._parameters.outer_tubing_diameter * 0.5
        height_2 = self._parameters.total_length - distance_2
        body_pt2 = occ.Cylinder(p=point_2, d=direction, r=radius_2, h=height_2)
        body = occ.Fuse([body_pt1, body_pt2])
        body.bc(self._boundaries['Body'])
        return body

    def __contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        direction = self._direction
        radius_1 = self._parameters.core_electrode_diameter * 0.5
        center = tuple(np.array(self._direction) * radius_1)
        contatct_tip = occ.Sphere(c=center, r=radius_1)
        height = self._parameters.core_electrode_length - radius_1
        contatct = occ.Cylinder(p=center, d=direction, r=radius_1, h=height)
        contact_1 = contatct_tip + contatct

        distance = self._parameters.core_electrode_length + self._parameters.core_tubing_length
        point = tuple(np.array(self._direction) * distance)
        radius_2 = self._parameters.outer_electrode_diameter * 0.5
        height_2 = self._parameters.outer_electrode_length
        contact_2 = occ.Cylinder(p=point, d=direction, r=radius_2, h=height_2)

        contact_1.bc(self._boundaries['Contact_1'])
        contact_2.bc(self._boundaries['Contact_2'])

        for edge in contact_1.edges:
            edge.name = self._boundaries['Contact_1']

        for edge in contact_2.edges:
            edge.name = self._boundaries['Contact_2']

        return netgen.occ.Glue([contact_1, contact_2])

    def get_max_mesh_size_contacts(self, ratio: float) -> float:
        """Use electrode's contact size to estimate maximal mesh size.

        Parameters
        ----------

        ratio: float
            Ratio between characteristic contact size and maximal mesh size.

        """

        return self._parameters.core_electrode_diameter / ratio
