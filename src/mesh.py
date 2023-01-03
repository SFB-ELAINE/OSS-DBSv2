
from typing import List
from src.brain_imaging.magnetic_resonance_imaging import MagneticResonanceImage
from src.brainsubstance import Material
from src.voxels import Voxels
import ngsolve
import numpy as np


class Mesh:

    def __init__(self,
                 geometry,
                 order: int) -> None:

        self.__mesh = ngsolve.Mesh(ngmesh=geometry.GenerateMesh())
        self.__mesh.Curve(order=order)
        self.__order = order

    def boundaries(self, name: str) -> ngsolve.comp.Region:
        return self.__mesh.Boundaries(pattern=name)

    def get_boundaries(self):
        return list(set(self.__mesh.GetBoundaries()) - set(['default']))

    def boundary_coefficients(self, boundaries) \
            -> ngsolve.fem.CoefficientFunction:
        return self.__mesh.BoundaryCF(values=boundaries)

    def flux_space(self, complex: bool = True) -> ngsolve.comp.HDiv:
        return ngsolve.HDiv(mesh=self.__mesh,
                            order=self.__order-1,
                            complex=complex)

    def materials(self) -> tuple:
        return self.__mesh.GetMaterials()

    def ngsolvemesh(self) -> ngsolve.comp.Mesh:
        return self.__mesh

    def refine(self) -> None:
        self.__mesh.Refine()
        self.__mesh.Curve(order=self.__order)

    def elements_with_error(self,
                            error: ngsolve.fem.CoefficientFunction) -> List:
        errors = ngsolve.Integrate(cf=error,
                                   mesh=self.__mesh,
                                   VOL_or_BND=ngsolve.VOL,
                                   element_wise=True).real
        limit = 0.5 * max(errors)
        return [errors[el.nr] > limit for el in self.__mesh.Elements()]

    def __elements_at_position(self, position: Voxels) -> None:
        space = ngsolve.L2(self.__mesh, order=0)
        grid_function = ngsolve.GridFunction(space=space)
        cf = ngsolve.VoxelCoefficient(start=tuple(position.start),
                                      end=tuple(position.end),
                                      values=position.data.astype(float),
                                      linear=False)
        grid_function.Set(cf)
        return grid_function.vec.FV().NumPy()

    def element_sizes(self) -> list:
        cf = ngsolve.CoefficientFunction(1)
        volumes = ngsolve.Integrate(cf=cf, mesh=self.__mesh, element_wise=True)
        return (6 * volumes.NumPy()) ** (1 / 3)

    def sobolev_space(self, complex: bool = False) -> ngsolve.comp.H1:
        dirichlet = '|'.join(boundary for boundary in self.get_boundaries())
        return ngsolve.H1(mesh=self.__mesh,
                          order=self.__order,
                          dirichlet=dirichlet,
                          complex=complex,
                          wb_withedges=False)

    def set_volume_refinement_flags(self, flags: List[bool]) -> None:
        for element, flag in zip(self.__mesh.Elements(ngsolve.VOL), flags):
            self.__mesh.SetRefinementFlag(ei=element, refine=flag)

    def refine_by_mri(self, mri: MagneticResonanceImage) -> None:
        maximum_size = min(mri.voxel_size())
        csf_voxel = mri.material_distribution(Material.CSF)
        flags = np.logical_and(self.__elements_at_position(csf_voxel),
                               self.element_sizes() > maximum_size)
        while np.any(flags) and self.sobolev_space().ndof < 1e5:
            self.set_volume_refinement_flags(flags)
            self.refine()
            csf_voxel = mri.material_distribution(Material.CSF)
            flags = np.logical_and(self.__elements_at_position(csf_voxel),
                                   self.element_sizes() > maximum_size)

    def refine_by_boundaries(self, boundaries: list) -> None:
        elements = self.__mesh.Elements(ngsolve.BND)
        flags = [element.mat in boundaries for element in elements]
        for element, flag in zip(self.__mesh.Elements(ngsolve.BND), flags):
            self.__mesh.SetRefinementFlag(ei=element, refine=flag)
        self.refine()
