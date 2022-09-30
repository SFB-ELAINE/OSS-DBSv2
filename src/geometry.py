
import netgen.occ as occ

class Electrode:

    def __init__(self) -> None:
        radius = 225e-6 / 2 * 1.1 # [m]   
        direction = occ.Z
        height = 0.02
        shaft = occ.Cylinder(p=(-0.005, 0 , radius), 
                                d=direction,
                                r=radius, 
                                h=height)
        shaft.bc('shaft')
        shaft.maxh = radius * 2 * 50
        contact = occ.Sphere(c=(-0.005, 0, 0), r=radius) - shaft
        contact.bc('contact')
        contact.maxh = 50 * radius / 20 
        self.solid = occ.Glue([contact, shaft])



class Wire:

    def __init__(self) -> None:
        r_wire = 250e-6 / 2
        height = 0.02
        wire = occ.Cylinder((0.005,0,0), occ.Z, r=r_wire, h=height)
        wire.bc('wire')
        wire.maxh = 50 * r_wire
        self.solid = wire


class Brain:

    def __init__(self) -> None:
        sphere = occ.Sphere((0,0,0), 0.02)
        sphere.bc('brain')
        sphere.mat('saline')
        self.solid = sphere


class SimpleGeometry:

    def __init__(self) -> None:
        electrode = Electrode()
        wire = Wire()
        brain = Brain()
        model = brain.solid - electrode.solid - wire.solid
        self.__geometry = occ.OCCGeometry(model)

    def generate_mesh(self):
        return self.__geometry.GenerateMesh()