# MicroProbes SNEX 100 Concentric Bipolar Electrodes
from ossdbs.electrodes.electrode import Electrode
import netgen
import numpy as np


class MicroProbesSNEX_100(Electrode):
    """MicroProbes SNEX 100 Concentric Bipolar electrode.

    Attributes
    ----------
    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    translation : tuple
        Translation vector (x,y,z) of electrode.

    Methods
    -------
    generate_geometry()
        Generate geometry of electrode.
    """
    # dimensions [mm]
    CORE_ELECTRODE_LENGTH = 0.25
    CORE_ELECTRODE_DIAMETER = 0.1
    CORE_TUBING_LENGTH = 0.5
    CORE_TUBING_DIAMETER = 0.140
    OUTER_ELECTRODE_LENGTH = 0.25
    OUTER_ELECTRODE_DIAMETER = 0.330
    OUTER_TUBING_DIAMETER = 0.411
    TOTAL_LENGTH = 10.0

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 translation: tuple = (0, 0, 0)) -> None:
        self.__translation = translation
        norm = np.linalg.norm(direction)
        self.__direction = tuple(direction / norm) if norm else (0, 0, 1)
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

    def generate_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of electrode.

        Returns
        -------
        geometry : netgen.libngpy._NgOCC.TopoDS_Shape
        """
        return self.__construct_geometry()

    def __construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        contacts = self.__contacts()
        body = self.__body()
        electrode = netgen.occ.Glue([body, contacts])
        moved_electrode = electrode.Move(self.__translation)
        return moved_electrode

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:

        distance_1 = self.CORE_ELECTRODE_LENGTH
        point_1 = tuple(np.array(self.__direction) * distance_1)
        body_1 = netgen.occ.Cylinder(p=point_1,
                                     d=self.__direction,
                                     r=self.CORE_TUBING_DIAMETER * 0.5,
                                     h=self.CORE_TUBING_LENGTH)

        distance_2 = (self.CORE_ELECTRODE_LENGTH
                      + self.CORE_TUBING_LENGTH
                      + self.OUTER_ELECTRODE_LENGTH)
        point_2 = tuple(np.array(self.__direction) * distance_2)
        body_2 = netgen.occ.Cylinder(p=point_2,
                                     d=self.__direction,
                                     r=self.OUTER_TUBING_DIAMETER * 0.5,
                                     h=self.TOTAL_LENGTH - distance_2)

        body = netgen.occ.Fuse([body_1, body_2])
        body.bc(self.__boundaries['Body'])
        return body

    def __contacts(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        inner_radius = self.CORE_ELECTRODE_DIAMETER * 0.5
        center = tuple(np.array(self.__direction) * inner_radius)
        inner_contatct_tip = netgen.occ.Sphere(c=center, r=inner_radius)
        height = self.CORE_ELECTRODE_LENGTH - inner_radius
        inner_contatct_pt2 = netgen.occ.Cylinder(p=center,
                                                 d=self.__direction,
                                                 r=inner_radius,
                                                 h=height)

        distance = self.CORE_ELECTRODE_LENGTH + self.CORE_TUBING_LENGTH
        point = tuple(np.array(self.__direction) * distance)
        outer_radius = self.OUTER_ELECTRODE_DIAMETER * 0.5
        outer_contact = netgen.occ.Cylinder(p=point,
                                            d=self.__direction,
                                            r=outer_radius,
                                            h=self.OUTER_ELECTRODE_LENGTH)

        inner_contact = inner_contatct_tip + inner_contatct_pt2
        inner_contact.bc(self.__boundaries['Contact_1'])
        outer_contact.bc(self.__boundaries['Contact_2'])

        return netgen.occ.Glue([inner_contact, outer_contact])
