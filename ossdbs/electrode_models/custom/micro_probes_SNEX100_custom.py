# MicroProbes SNEX 100 Concentric Bipolar Electrodes
from ossdbs.electrode_models.electrode import ElectrodeModel
import netgen
import netgen.occ as occ
import numpy as np
import os
import json


class MicroProbesSNEX_100Custom(ElectrodeModel):
    """MicroProbes SNEX 100 Concentric Bipolar electrode.

    Attributes
    ----------
    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    translation : tuple
        Translation vector (x,y,z) of electrode.
    """
    # dimensions [m]
    CORE_ELECTRODE_LENGTH = 0.25
    CORE_ELECTRODE_DIAMETER = 0.1
    CORE_TUBING_LENGTH = 0.5
    CORE_TUBING_DIAMETER = 0.140
    OUTER_ELECTRODE_LENGTH = 0.25
    OUTER_ELECTRODE_DIAMETER = 0.330
    OUTER_TUBING_DIAMETER = 0.411
    TOTAL_LENGTH = 100.0

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        self.__position = position
        norm = np.linalg.norm(direction)
        self.__direction = tuple(direction / norm) if norm else (0, 0, 1)
        self.__load_parameters()
        self.__boundaries = {'Body': 'Body',
                             'Contact_1': 'Contact_1',
                             'Contact_2': 'Contact_2',
                             'Contact_3': 'Contact_3',
                             'Contact_4': 'Contact_4',
                             'Contact_5': 'Contact_5',
                             'Contact_6': 'Contact_6',
                             'Contact_7': 'Contact_7',
                             'Contact_8': 'Contact_8'}

    def rename_boundaries(self, boundaries: dict) -> None:
        self.__boundaries.update(boundaries)

    def geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of electrode.

        Returns
        -------
        geometry : netgen.libngpy._NgOCC.TopoDS_Shape
        """
        return self.__construct_geometry()

    def __construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        electrode = occ.Glue([self.__body(), self.__contacts()])
        return electrode.Move(v=self.__position)

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        direction = self.__direction
        distance_1 = self.CORE_ELECTRODE_LENGTH
        point_1 = tuple(np.array(self.__direction) * distance_1)
        radius_1 = self.CORE_TUBING_DIAMETER * 0.5
        height_1 = self.CORE_TUBING_LENGTH
        body_pt1 = occ.Cylinder(p=point_1, d=direction, r=radius_1, h=height_1)

        distance_2 = (self.CORE_ELECTRODE_LENGTH
                      + self.CORE_TUBING_LENGTH
                      + self.OUTER_ELECTRODE_LENGTH)
        point_2 = tuple(np.array(direction) * distance_2)
        radius_2 = self.OUTER_TUBING_DIAMETER * 0.5
        height_2 = self.TOTAL_LENGTH - distance_2
        body_pt2 = occ.Cylinder(p=point_2, d=direction, r=radius_2, h=height_2)
        body = occ.Fuse([body_pt1, body_pt2])
        body.bc(self.__boundaries['Body'])
        return body

    def __contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        direction = self.__direction
        radius_1 = self.CORE_ELECTRODE_DIAMETER * 0.5
        center = tuple(np.array(self.__direction) * radius_1)
        contatct_tip = occ.Sphere(c=center, r=radius_1)
        height = self.CORE_ELECTRODE_LENGTH - radius_1
        contatct_pt2 = occ.Cylinder(p=center, d=direction, r=radius_1, h=height)
        contact_1 = contatct_tip + contatct_pt2

        distance = self.CORE_ELECTRODE_LENGTH + self.CORE_TUBING_LENGTH
        point = tuple(np.array(self.__direction) * distance)
        radius_2 = self.OUTER_ELECTRODE_DIAMETER * 0.5
        height_2 = self.OUTER_ELECTRODE_LENGTH
        contact_2 = occ.Cylinder(p=point, d=direction, r=radius_2, h=height_2)

        contact_1.bc(self.__boundaries['Contact_1'])
        contact_2.bc(self.__boundaries['Contact_2'])

        for edge in contact_1.edges:
            edge.name = self.__boundaries['Contact_1']

        for edge in contact_2.edges:
            edge.name = self.__boundaries['Contact_2']

        return netgen.occ.Glue([contact_1, contact_2])

    def capsule_geometry(self, thickness: float, max_h: float = 0.1) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:
        direction = self.__direction
        distance_1 = self.CORE_ELECTRODE_DIAMETER * 0.5
        point_1 = tuple(np.array(self.__direction) * distance_1)
        radius_1 = self.CORE_ELECTRODE_DIAMETER * 0.5 + thickness
        part_0 = occ.Sphere(c=point_1, r=radius_1)
        height_1 = self.CORE_ELECTRODE_LENGTH - distance_1
        part_1 = occ.Cylinder(p=point_1, d=direction, r=radius_1, h=height_1)

        distance_2 = self.CORE_ELECTRODE_LENGTH - thickness
        point_2 = tuple(np.array(self.__direction) * distance_2)
        radius_2 = self.CORE_TUBING_DIAMETER * 0.5 + thickness
        height_2 = self.CORE_TUBING_LENGTH
        part_2 = occ.Cylinder(p=point_2, d=direction, r=radius_2, h=height_2)

        distance_3 = distance_2 + self.CORE_TUBING_LENGTH
        point_3 = tuple(np.array(self.__direction) * distance_3)
        radius_3 = self.OUTER_ELECTRODE_DIAMETER * 0.5 + thickness
        height_3 = self.OUTER_ELECTRODE_LENGTH
        part_3 = occ.Cylinder(p=point_3, d=direction, r=radius_3, h=height_3)

        distance_4 = distance_3 + self.OUTER_ELECTRODE_LENGTH
        point_4 = tuple(np.array(self.__direction) * distance_4)
        radius_4 = self.OUTER_TUBING_DIAMETER * 0.5 + thickness
        height_4 = self.TOTAL_LENGTH - distance_4
        part_4 = occ.Cylinder(p=point_4, d=direction, r=radius_4, h=height_4)

        capsule = part_0 + part_1 + part_2 + part_3 + part_4
        capsule.bc('Capsule')
        capsule.mat('Capsule')
        capsule.maxh = max_h
        return capsule.Move(v=self.__position)

    def __load_parameters(self) -> None:
        dir_name = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(dir_name, 'micro_probes_SNEX100_custom.json')

        with open(path, 'r') as json_file:
            parameters = json.load(json_file)

        self.CORE_ELECTRODE_LENGTH = parameters['CoreElectrodeLength[mm]']
        self.CORE_ELECTRODE_DIAMETER = parameters['CoreElectrodeDiameter[mm]']
        self.CORE_TUBING_LENGTH = parameters['CoreTubingLength[mm]']
        self.CORE_TUBING_DIAMETER = parameters['CoreTubingDiameter[mm]']
        self.OUTER_ELECTRODE_LENGTH = parameters['OuterElectrodeLength[mm]']
        self.OUTER_ELECTRODE_DIAMETER = parameters['OuterElectrodeDiameter[mm]']
        self.OUTER_TUBING_DIAMETER = parameters['OuterTubingDiameter[mm]']
        self.TOTAL_LENGTH = parameters['LeadLength[mm]']
