# MicroProbes SNEX 100 Concentric Bipolar Electrodes
from dataclasses import dataclass
from .electrode_model_template import ElectrodeModel
import netgen
import netgen.occ as occ
import numpy as np


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

    def __init__(self,
                 parameters: MicroProbesSNEX100Parameters,
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
                            'Contact_4': 'Contact_4',
                            'Contact_5': 'Contact_5',
                            'Contact_6': 'Contact_6',
                            'Contact_7': 'Contact_7',
                            'Contact_8': 'Contact_8'}
        self._core_electrode_length = parameters.core_electrode_length
        self._core_electrode_diameter = parameters.core_electrode_diameter
        self._core_tubing_length = parameters.core_tubing_length
        self._core_tubing_diameter = parameters.core_tubing_diameter
        self._outer_electrode_length = parameters.outer_electrode_length
        self._outer_electrode_diameter = parameters.outer_electrode_diameter
        self._outer_tubing_diameter = parameters.outer_tubing_diameter
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
        direction = self._direction
        distance_1 = self._core_electrode_diameter * 0.5
        point_1 = tuple(np.array(self._direction) * distance_1)
        radius_1 = self._core_electrode_diameter * 0.5 + thickness
        part_0 = occ.Sphere(c=point_1, r=radius_1)
        height_1 = self._core_electrode_length - distance_1
        part_1 = occ.Cylinder(p=point_1, d=direction, r=radius_1, h=height_1)

        distance_2 = self._core_electrode_length - thickness
        point_2 = tuple(np.array(self._direction) * distance_2)
        radius_2 = self._core_tubing_diameter * 0.5 + thickness
        height_2 = self._core_tubing_length
        part_2 = occ.Cylinder(p=point_2, d=direction, r=radius_2, h=height_2)

        distance_3 = distance_2 + self._core_tubing_length
        point_3 = tuple(np.array(self._direction) * distance_3)
        radius_3 = self._outer_electrode_diameter * 0.5 + thickness
        height_3 = self._outer_electrode_length
        part_3 = occ.Cylinder(p=point_3, d=direction, r=radius_3, h=height_3)

        distance_4 = distance_3 + self._outer_electrode_length
        point_4 = tuple(np.array(self._direction) * distance_4)
        radius_4 = self._outer_tubing_diameter * 0.5 + thickness
        height_4 = self._total_length - distance_4
        part_4 = occ.Cylinder(p=point_4, d=direction, r=radius_4, h=height_4)

        capsule = part_0 + part_1 + part_2 + part_3 + part_4
        capsule.bc('Capsule')
        capsule.mat('Capsule')
        return capsule.Move(v=self._position) - self.geometry()

    def geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of electrode.

        Returns
        -------
        netgen.libngpy._NgOCC.TopoDS_Shape
        """
        return self.__construct_geometry()

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

    def __construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        electrode = occ.Glue([self.__body(), self.__contacts()])
        return electrode.Move(v=self._position)

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        direction = self._direction
        distance_1 = self._core_electrode_length
        point_1 = tuple(np.array(self._direction) * distance_1)
        radius_1 = self._core_tubing_diameter * 0.5
        height_1 = self._core_tubing_length
        body_pt1 = occ.Cylinder(p=point_1, d=direction, r=radius_1, h=height_1)

        distance_2 = (self._core_electrode_length
                      + self._core_tubing_length
                      + self._outer_electrode_length)
        point_2 = tuple(np.array(direction) * distance_2)
        radius_2 = self._outer_tubing_diameter * 0.5
        height_2 = self._total_length - distance_2
        body_pt2 = occ.Cylinder(p=point_2, d=direction, r=radius_2, h=height_2)
        body = occ.Fuse([body_pt1, body_pt2])
        body.bc(self._boundaries['Body'])
        return body

    def __contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        direction = self._direction
        radius_1 = self._core_electrode_diameter * 0.5
        center = tuple(np.array(self._direction) * radius_1)
        contatct_tip = occ.Sphere(c=center, r=radius_1)
        height = self._core_electrode_length - radius_1
        contatct = occ.Cylinder(p=center, d=direction, r=radius_1, h=height)
        contact_1 = contatct_tip + contatct

        distance = self._core_electrode_length + self._core_tubing_length
        point = tuple(np.array(self._direction) * distance)
        radius_2 = self._outer_electrode_diameter * 0.5
        height_2 = self._outer_electrode_length
        contact_2 = occ.Cylinder(p=point, d=direction, r=radius_2, h=height_2)

        contact_1.bc(self._boundaries['Contact_1'])
        contact_2.bc(self._boundaries['Contact_2'])

        for edge in contact_1.edges:
            edge.name = self._boundaries['Contact_1']

        for edge in contact_2.edges:
            edge.name = self._boundaries['Contact_2']

        return netgen.occ.Glue([contact_1, contact_2])
