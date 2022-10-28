# Boston Scientific (Marlborough, Massachusetts, USA) vercise
# from .electrode import Electrode
import netgen.occ as occ
import netgen
import numpy as np

class BostonScientificVercise():

    def __init__(self,
                 rotation: float = 0,
                 direction: tuple = (0, 0, 1),
                 translation: tuple = (0, 0, 0)):
        self.__rotation = rotation
        self.__direction = direction
        self.__translation = translation

    def build(self):
        tip_length = 1.1
        contact_length = 1.5
        contact_spacing = 0.5
        lead_diameter = 1.3
        length = 50

        lead_radius = lead_diameter / 2

        d_norm = np.array(self.__direction) / np.linalg.norm(self.__direction)

        center = tuple(d_norm * lead_radius)
        tip = occ.Sphere(c=center, r=lead_radius)
        lead = occ.Cylinder(p=(0, 0, lead_radius),
                            d=occ.Z,
                            r=lead_radius,
                            h=length)
        body = tip + lead
        body.bc("body")

        contact = occ.Cylinder(p=(0, 0, 0),
                               d=occ.Z,
                               r=lead_radius,
                               h=contact_length)

        h1 = tip_length
        contact_1 = contact.Move((0, 0, h1))
        contact_1.bc("Contact_1")

        h2 = h1 + contact_length + contact_spacing
        contact_2 = contact.Move((0, 0, h2))
        contact_2.bc("Contact_2")

        h3 = h2 + contact_length + contact_spacing
        contact_3 = contact.Move((0, 0, h3))
        contact_3.bc("Contact_3")

        h4 = h3 + contact_length + contact_spacing
        contact_4 = contact.Move((0, 0, h4))
        contact_4.bc("Contact_4")

        h5 = h4 + contact_length + contact_spacing
        contact_5 = contact.Move((0, 0, h5))
        contact_5.bc("Contact_5")

        h6 = h5 + contact_length + contact_spacing
        contact_6 = contact.Move((0, 0, h6))
        contact_6.bc("Contact_6")

        h7 = h6 + contact_length + contact_spacing
        contact_7 = contact.Move((0, 0, h7))
        contact_7.bc("Contact_7")

        h8 = h7 + contact_length + contact_spacing
        contact_8 = contact.Move((0, 0, h8))
        contact_8.bc("Contact_8")

        electrode = occ.Glue([body,
                              contact_1,
                              contact_2,
                              contact_3,
                              contact_4,
                              contact_5,
                              contact_6,
                              contact_7,
                              contact_8])

        moved_electrode = electrode.Move(self.__translation)

        axis = netgen.libngpy._NgOCC.Axis(p=(0, 0, 0), d=self.__direction)

        rotated_electrode = moved_electrode.Rotate(axis=axis,
                                                   ang=self.__rotation)
    
        return rotated_electrode



import ngsolve

order = 2
electrode = BostonScientificVercise()
geo = occ.OCCGeometry(electrode.build((10,5,10)))
with ngsolve.TaskManager():
    mesh = ngsolve.Mesh(geo.GenerateMesh())
    mesh.Curve(order)
print(mesh.GetBoundaries())
print(mesh.GetMaterials())
bnd_dict = {}
for i in range(1, 9):
    bnd_dict["Contact_{}".format(i)] = float(i)
bndcf = mesh.BoundaryCF(bnd_dict, default=-1)
ngsolve.Draw(bndcf, mesh, "BND")