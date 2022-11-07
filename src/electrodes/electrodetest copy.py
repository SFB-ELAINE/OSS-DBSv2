# Boston Scientific (Marlborough, Massachusetts, USA) vercise
# from .electrode import Electrode
import netgen
import numpy as np


class TestElectrode():

    # dimensions [mm]
    TIP_LENGTH = 1.1
    CONTACT_LENGTH = 1.5
    CONTACT_SPACING = 0.5
    LEAD_DIAMETER = 1.3
    TOTAL_LENGHTH = 20.0
    TUBE_THICKNESS = 0.0
    TUBE_FREE_LENGTH = 20.0
    CONTACT_SPACING_RADIAL = 0.25

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 translation: tuple = (0, 0, 0)) -> None:
        self.__translation = translation
        norm = np.linalg.norm(direction)
        self.__direction = tuple(direction / norm) if norm else (0, 0, 1)

    def generate_geometry(self) -> netgen.libngpy._meshing.Mesh:
        """Generate geometry of electrode.

        Returns
        -------
        geometry : netgen.libngpy._NgOCC.TopoDS_Shape
        """
        return self.__construct_geometry()

    def __construct_geometry(self) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:
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
        height = self.TOTAL_LENGHTH - self.TIP_LENGTH
        body = netgen.occ.Cylinder(p=center,
                                   d=self.__direction,
                                   r=diameter * 0.5,
                                   h=height)
        body.bc("Body")
        return body

    def __contacts(self, diameter: float) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:

        point = tuple(np.array(self.__direction) * self.TIP_LENGTH)
        space = netgen.occ.HalfSpace(p=point, n=self.__direction)
        center = tuple(np.array(self.__direction) * self.LEAD_DIAMETER * 0.5)
        active_tip_pt1 = netgen.occ.Sphere(c=center, r=diameter * 0.5) * space
        height = self.TIP_LENGTH - self.LEAD_DIAMETER * 0.5
        active_tip_pt2 = netgen.occ.Cylinder(p=center,
                                             d=self.__direction,
                                             r=diameter * 0.5,
                                             h=height)
        active_tip = netgen.occ.Fuse([active_tip_pt1, active_tip_pt2])
        contact_8 = netgen.occ.Cylinder(p=(0, 0, 0),
                                        d=self.__direction,
                                        r=diameter * 0.5,
                                        h=self.CONTACT_LENGTH)

        distance = (self.TIP_LENGTH +
                    3 * self.CONTACT_SPACING +
                    2 * self.CONTACT_LENGTH)

        contact_8 = contact_8.Move(tuple(np.array(self.__direction) *
                                         distance))

        axis = netgen.occ.Axis((0, 0, 0), self.__direction)
        distance = self.TIP_LENGTH + self.CONTACT_SPACING
        point = tuple(np.array(self.__direction) * distance)
        contact = self.__contact(diameter=diameter).Move(point)
        distance = self.CONTACT_LENGTH + self.CONTACT_SPACING
        point = tuple(np.array(self.__direction) * distance)
        contacts = [active_tip,
                    contact,
                    contact.Rotate(axis, 120),
                    contact.Rotate(axis, 240),
                    contact.Move(point),
                    contact.Rotate(axis, 120).Move(point),
                    contact.Rotate(axis, 240).Move(point),
                    contact_8
                    ]

        for index, contact in enumerate(contacts, 1):
            contact.bc("Contact_{}".format(index))

        return netgen.occ.Fuse(contacts)

    def __contact(self, diameter: float) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = diameter * 0.5
        body = netgen.occ.Cylinder(p=(0, 0, 0),
                                   d=self.__direction,
                                   r=radius,
                                   h=self.CONTACT_LENGTH)

        direction2 = self.__secound_direction()
        new_direction = tuple(np.cross(direction2, self.__direction))
        eraser = netgen.occ.HalfSpace(p=(0, 0, 0), n=new_direction)
        delta = self.CONTACT_SPACING_RADIAL / radius * (180 / np.pi)
        angle = 30 + delta * 0.5
        axis = netgen.occ.Axis((0, 0, 0), self.__direction)
        return body - eraser.Rotate(axis, angle) - eraser.Rotate(axis, -angle)

    def __secound_direction(self):
        x, y, z = self.__direction

        if not x and not y:
            return (0, 1, 0)

        if not x and not z:
            return (0, 0, 1)

        if not y and not z:
            return (0, 1, 0)

        return (x, y, not z)

    def __tube(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:

        if self.TUBE_FREE_LENGTH - self.TIP_LENGTH >= 0:
            radius = self.LEAD_DIAMETER * 0.5 + self.TUBE_THICKNESS
            height = self.TOTAL_LENGHTH - self.TUBE_FREE_LENGTH
            point = tuple(np.array(self.__direction) * self.TUBE_FREE_LENGTH)
            tube = netgen.occ.Cylinder(p=point,
                                       d=self.__direction,
                                       r=radius,
                                       h=height)
            tube.bc('Contact_1')

        if (self.TIP_LENGTH + self.CONTACT_SPACING <
            self.TUBE_FREE_LENGTH <
            self.TIP_LENGTH + self.CONTACT_SPACING + self.CONTACT_LENGTH):

            axis = netgen.occ.Axis((0, 0, 0), self.__direction)
            diameter = self.LEAD_DIAMETER + 2 * self.TUBE_THICKNESS
            contact = self.__contact(diameter=diameter)
            height = self.TOTAL_LENGHTH - self.TUBE_FREE_LENGTH
            distance = self.CONTACT_LENGTH + self.CONTACT_SPACING
            contacts = [contact,
                        contact.Rotate(axis, 120),
                        contact.Rotate(axis, 240),
                        ]

            for index, contact in enumerate(contacts, 2):
                contact.bc("Contact_{}".format(index))

        if (self.TIP_LENGTH + 3 * self.CONTACT_SPACING + 
            2 * self.CONTACT_LENGTH <
            self.TUBE_FREE_LENGTH <
            self.TIP_LENGTH + 3 * self.CONTACT_SPACING + 3 * self.CONTACT_LENGTH):
            
            axis = netgen.occ.Axis((0, 0, 0), self.__direction)
            diameter = self.LEAD_DIAMETER + 2 * self.TUBE_THICKNESS
            contact = self.__contact(diameter=diameter)
            height = self.TOTAL_LENGHTH - self.TUBE_FREE_LENGTH
            distance = self.CONTACT_LENGTH + self.CONTACT_SPACING
            contacts = [contact,
                        contact.Rotate(axis, 120),
                        contact.Rotate(axis, 240),
                        ]

            for index, contact in enumerate(contacts, 5):
                contact.bc("Contact_{}".format(index))

    
        if (self.TIP_LENGTH + 3 * self.CONTACT_SPACING + 
            2 * self.CONTACT_LENGTH <
            self.TUBE_FREE_LENGTH <
            self.TIP_LENGTH + 3 * self.CONTACT_SPACING + 3 * self.CONTACT_LENGTH):
            radius = self.LEAD_DIAMETER * 0.5 + self.TUBE_THICKNESS
            height = self.TOTAL_LENGHTH - self.TUBE_FREE_LENGTH
            point = tuple(np.array(self.__direction) * self.TUBE_FREE_LENGTH)
            tube = netgen.occ.Cylinder(p=point,
                                    d=self.__direction,
                                    r=radius,
                                    h=height)
            tube.bc('Contact_8')
        else:
            radius = self.LEAD_DIAMETER * 0.5 + self.TUBE_THICKNESS
            height = self.TOTAL_LENGHTH - self.TUBE_FREE_LENGTH
            point = tuple(np.array(self.__direction) * self.TUBE_FREE_LENGTH)
            tube = netgen.occ.Cylinder(p=point,
                                    d=self.__direction,
                                    r=radius,
                                    h=height)
            tube.bc('Body')

        return tube



class AbbottStjudeDirected6172():
    """Abbott/St Jude Directed 6172 electrode.

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
    TIP_LENGTH = 1.1
    CONTACT_LENGTH = 1.5
    CONTACT_SPACING = 0.5
    LEAD_DIAMETER = 1.3
    TOTAL_LENGHTH = 20.0
    TUBE_THICKNESS = 1.0
    TUBE_FREE_LENGTH = 1.7
    N_CONTACTS = 4

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
        point = tuple(np.array(self.__direction) * self.LEAD_DIAMETER * 0.5)
        height = self.TOTAL_LENGHTH - self.TIP_LENGTH
        body = netgen.occ.Cylinder(p=point,
                                   d=self.__direction,
                                   r=diameter * 0.5,
                                   h=height)
        body.bc("Body")
        return body

    def __contacts(self, diameter: float)\
            -> netgen.libngpy._NgOCC.TopoDS_Shape:

        point = tuple(np.array(self.__direction) * self.TIP_LENGTH)
        space = netgen.occ.HalfSpace(p=point, n=self.__direction)
        center = tuple(np.array(self.__direction) * self.LEAD_DIAMETER * 0.5)
        tip = netgen.occ.Sphere(c=center, r=diameter * 0.5) * space
        height = max(self.TIP_LENGTH - self.LEAD_DIAMETER * 0.5, 0)

        if height > 0:
            contact_1 = netgen.occ.Cylinder(p=center,
                                            d=self.__direction,
                                            r=diameter * 0.5,
                                            h=height)
            active_tip = netgen.occ.Fuse([tip, contact_1])
        else:
            active_tip = tip

        active_tip.bc('Contact_1')

        contact = netgen.occ.Cylinder(p=(0, 0, 0),
                                      d=self.__direction,
                                      r=diameter * 0.5,
                                      h=self.CONTACT_LENGTH)
        length = (self.CONTACT_LENGTH + self.CONTACT_SPACING)
        offset = self.TIP_LENGTH + self.CONTACT_SPACING
        contacts = [contact.Move(tuple(np.array(self.__direction) * distance))
                    for distance
                    in np.arange(self.N_CONTACTS - 1) * length + offset]

        for index, contact in enumerate(contacts, 2):
            contact.bc("Contact_{}".format(index))

        return netgen.occ.Fuse([active_tip] + contacts)

    def __tube(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self.LEAD_DIAMETER * 0.5 + self.TUBE_THICKNESS
        height = self.TOTAL_LENGHTH - self.TUBE_FREE_LENGTH
        point = tuple(np.array(self.__direction) * self.TUBE_FREE_LENGTH)
        tube = netgen.occ.Cylinder(p=point,
                                   d=self.__direction,
                                   r=radius,
                                   h=height)
        tube.bc('Body')

        offset = self.TIP_LENGTH + self.CONTACT_SPACING
        distance = self.CONTACT_LENGTH + self.CONTACT_SPACING
        distances = offset + np.arange(self.N_CONTACTS - 1) * distance
        lower_limit = np.append(0, distances)
        upper_limit = np.append(self.TIP_LENGTH,
                                distances + self.CONTACT_LENGTH)
        intersection = np.logical_and(lower_limit < self.TUBE_FREE_LENGTH,
                                      upper_limit > self.TUBE_FREE_LENGTH)
        if np.any(intersection):
            tube.bc('Contact_{}'.format(np.argmax(intersection) + 1))

        return tube


import ngsolve

#electrode = AbbottStjudeDirected6172(direction=(0, 0, 1))
electrode = TestElectrode(direction=(0, 0, 1))
geometry = electrode.generate_geometry()

with ngsolve.TaskManager():
    mesh = ngsolve.Mesh(netgen.occ.OCCGeometry(geometry).GenerateMesh())
print('Boundaries:', set(mesh.GetBoundaries()))
print()
print('Materials:', mesh.GetMaterials())
bnd_dict = {"Contact_{}".format(i): float(i) for i in range(1, 9)}
bndcf = mesh.BoundaryCF(bnd_dict, default=-1)
ngsolve.Draw(bndcf, mesh, "BND")
