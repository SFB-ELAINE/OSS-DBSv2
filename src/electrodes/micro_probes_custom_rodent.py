from src.electrodes.abstract_electrode import AbstractElectrode
import netgen
import numpy as np


class MicroProbesCustomRodent(AbstractElectrode):
    """MicroProbes Custom Rodent electrode.

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
        Generate mesh of electrode.
    """

    # dimensions [mm]
    CONTACT_LENGTH = 0.01125
    LEAD_DIAMETER = 0.225
    TOTAL_LENGHTH = 1.3
    TUBE_THICKNESS = .01

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 translation: tuple = (0, 0, 0),
                 boundaries: list = None) -> None:
        self.__boundaries = boundaries
        self.__translation = translation
        norm = np.linalg.norm(direction)
        self.__direction = tuple(direction / norm) if norm else (0, 0, 1)

    def generate_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of electrode.

        Returns
        -------
        geometry : netgen.libngpy._NgOCC.TopoDS_Shape
        """
        electrode = netgen.occ.Glue([self.__contact(), self.__body()])
        moved_electrode = electrode.Move(self.__translation)
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
        body.bc("Body")
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
        contact.bc(self.__boundaries[0][0])
        return contact

    def boundary_values(self):
        values = {key: value for key, value in self.__boundaries}
        values['Body'] = 0.0
        return values
