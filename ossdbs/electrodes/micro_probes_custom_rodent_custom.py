
from ossdbs.electrodes.electrode import Electrode
import netgen
import numpy as np
import os
import json


class MicroProbesCustomRodentCustom(Electrode):
    """MicroProbes Custom Rodent electrode.

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
    CONTACT_LENGTH = 0.01125e-3
    LEAD_DIAMETER = 0.225e-3
    TOTAL_LENGHTH = 13.3e-3
    TUBE_THICKNESS = .01e-3

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0)) -> None:
        self.__position = tuple(position)
        norm = np.linalg.norm(direction)
        self.__direction = tuple(direction / norm) if norm else (0, 0, 1)
        self.__load_parameters()
        self.__boundaries = {'Body': 'Body', 'Contact_1': 'Contact_1'}

    def rename_boundaries(self, boundaries: dict) -> None:
        self.__boundaries.update(boundaries)

    def generate_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of electrode.

        Returns
        -------
        geometry : netgen.libngpy._NgOCC.TopoDS_Shape
        """
        electrode = netgen.occ.Glue([self.__contact(), self.__body()])
        moved_electrode = electrode.Move(self.__position)
        return moved_electrode

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self.LEAD_DIAMETER * 0.5
        tube_radius = radius + self.TUBE_THICKNESS
        center = tuple(np.array(self.__direction) * radius)
        tip = netgen.occ.Sphere(c=center, r=tube_radius)
        lead = netgen.occ.Cylinder(p=center,
                                   d=self.__direction,
                                   r=tube_radius,
                                   h=self.TOTAL_LENGHTH - radius)
        point = tuple(np.array(self.__direction) * self.CONTACT_LENGTH)
        space = netgen.occ.HalfSpace(p=point, n=self.__direction)
        body = tip + lead - space
        body.bc(self.__boundaries['Body'])
        return body

    def __contact(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self.LEAD_DIAMETER * 0.5
        point = tuple(np.array(self.__direction) * self.CONTACT_LENGTH)
        space = netgen.occ.HalfSpace(p=point, n=self.__direction)
        center = tuple(np.array(self.__direction) * radius)
        tip = netgen.occ.Sphere(c=center, r=radius) * space
        lead = netgen.occ.Cylinder(p=center,
                                   d=self.__direction,
                                   r=radius,
                                   h=self.CONTACT_LENGTH - radius)

        if self.CONTACT_LENGTH <= radius:
            contact = tip
        else:
            contact = tip + lead
        contact.bc(self.__boundaries['Contact_1'])
        return contact

    def __load_parameters(self) -> None:
        dir_name = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(dir_name, 'medtronic_custom.json')

        with open(path, 'r') as json_file:
            parameters = json.load(json_file)

        self.CONTACT_LENGTH = parameters['ContactLength[m]']
        self.LEAD_DIAMETER = parameters['LeadDiameter[m]']
        self.TOTAL_LENGHTH = parameters['LeadLength[m]']
        self.TUBE_THICKNESS = parameters['TubeThickness[m]']
