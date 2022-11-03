from src.electrodes.abstract_electrode import AbstractElectrode
import netgen
import numpy as np


class Rodden(AbstractElectrode):
    """Rodden electrode.

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
    TIP_LENGTH = 0.1125
    CONTACT_LENGTH = 13
    LEAD_DIAMETER = 0.225
    TOTAL_LENGHTH = 100.0
    TUBE_THICKNESS = 0.0
    TUBE_FREE_LENGTH = 100.0

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 translation: tuple = (0, 0, 0)) -> None:
        self.__translation = translation
        norm = np.linalg.norm(direction)
        self.__direction = tuple(direction / norm) if norm else (0, 0, 1)

    def generate_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of electrode.

        Returns
        -------
        geometry : netgen.libngpy._NgOCC.TopoDS_Shape
        """
        return self.__construct_geometry()

    def __construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        tube = self.__tube()

        contacts = self.__contacts(diameter=self.LEAD_DIAMETER)
        body = self.__body(diameter=self.LEAD_DIAMETER) - contacts
        electrode = netgen.occ.Glue([body, contacts]) - tube

        diameter = self.LEAD_DIAMETER + 2 * self.TUBE_THICKNESS
        contacts_tube = self.__contacts(diameter=diameter)
        body_tube = self.__body(diameter=diameter) - contacts_tube
        electrode_tube = netgen.occ.Glue([body_tube, contacts_tube]) * tube

        electrode = netgen.occ.Glue([electrode, electrode_tube])
        moved_electrode = electrode.Move(self.__translation)

        return moved_electrode

    def __body(self, diameter: float) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        center = tuple(np.array(self.__direction) * self.LEAD_DIAMETER * 0.5)
        tip = netgen.occ.Sphere(c=center, r=diameter * 0.5)
        height = self.TOTAL_LENGHTH - self.TIP_LENGTH
        lead = netgen.occ.Cylinder(p=center,
                                   d=self.__direction,
                                   r=diameter * 0.5,
                                   h=height)
        body = tip + lead
        body.bc("Body")
        return body

    def __contacts(self, diameter: float) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:
        point = tuple(np.array(self.__direction) * self.TIP_LENGTH)
        contact = netgen.occ.Cylinder(p=point,
                                      d=self.__direction,
                                      r=diameter * 0.5,
                                      h=self.CONTACT_LENGTH)
        contact.bc("Contact_1")

        return contact

    def __tube(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self.LEAD_DIAMETER * 0.5 + self.TUBE_THICKNESS
        height = self.TOTAL_LENGHTH - self.TUBE_FREE_LENGTH
        point = tuple(np.array(self.__direction) * self.TUBE_FREE_LENGTH)
        tube = netgen.occ.Cylinder(p=point,
                                   d=self.__direction,
                                   r=radius,
                                   h=height)
        tube.bc('Body')
        lower_limit = self.TIP_LENGTH
        upper_limit = lower_limit + self.CONTACT_LENGTH

        if lower_limit < self.TUBE_FREE_LENGTH < upper_limit:
            tube.bc('Contact_1')

        return tube
