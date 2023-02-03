
from typing import List
from ossdbs.voxels import Voxels
import ngsolve
import numpy as np


class Mesh:
    """Class for interacting with the mesh for FEM.

    Parameters
    ----------
    geometry : netgen.libngpy._NgOCC.OCCGeometry

    order : int
        Order of mesh elements.

    complex_datatype : bool
            True for complex data type, False otherwise.
    """

    def __init__(self,
                 ngsolve_mesh: ngsolve.comp.Mesh,
                 order: int,
                 complex_datatype: bool = False) -> None:
        self.__mesh = ngsolve_mesh
        self.__mesh.Curve(order=order)
        self.__order = order
        self.__complex = complex_datatype

    def boundary_coefficients(self, boundaries) \
            -> ngsolve.fem.CoefficientFunction:
        """Return a boundary coefficient function.

        Returns
        -------
        ngsolve.fem.CoefficientFunction
        """

        return self.__mesh.BoundaryCF(values=boundaries)

    def flux_space(self) -> ngsolve.comp.HDiv:
        """Return a flux space based on the mesh.

        Returns
        -------
        ngsolve.comp.HDiv
        """

        return ngsolve.HDiv(mesh=self.__mesh,
                            order=self.__order-1,
                            complex=self.__complex)

    def ngsolvemesh(self) -> ngsolve.comp.Mesh:
        """Return mesh as a ngsolve object.

        Returns
        -------
        ngsolve.comp.Mesh
        """

        return self.__mesh

    def point_in_mesh(self, point: tuple) -> bool:
        """Check if point lays in mesh.

        Parameters
        ----------
        point : tuple
            Represents the coordinates x, y, z of 3D point.

        Returns
        -------
        bool
            True if point is inside mesh, False otherwise.
        """
        x, y, z = point
        return self.__mesh(x=x, y=y, z=z).nr != -1

    def included_points(self, points: np.ndarray) -> np.ndarray:
        """Check if point lays in mesh.

        Parameters
        ----------
        points : np.ndarray
            Collection of 3D point coordinates (x, y, z).

        Returns
        -------
        np.ndarray
            Collection of bool values:
            True if point is inside mesh, False otherwise.
        """
        x, y, z = points.T
        mips = self.__mesh(x, y, z)
        return points[np.array([mip[5] != -1 for mip in mips])]

    def refine(self) -> None:
        """Refine the mesh."""

        self.__mesh.Refine()
        self.__mesh.Curve(order=self.__order)

    def save(self, file_name: str) -> None:
        """Save netgen mesh.

        Parameters
        ----------
        file_name : str
            File name of the mesh data.
        """

        self.__mesh.ngmeshSave(file_name)

    def h1_space(self, boundaries: List[str]) -> ngsolve.comp.H1:
        """Return a h1 space based on the mesh.

        Returns
        -------
        ngsolve.comp.H1
        """

        dirichlet = '|'.join(boundary for boundary in boundaries)
        return ngsolve.H1(mesh=self.__mesh,
                          order=self.__order,
                          dirichlet=dirichlet,
                          complex=self.__complex,
                          wb_withedges=False)

    def number_space(self) -> ngsolve.comp.NumberSpace:
        """Return a number space based on the mesh.

        Returns
        -------
        ngsolve.comp.NumberSpace
            Space with only one single (global) DOF.
        """

        return ngsolve.NumberSpace(mesh=self.__mesh,
                                   order=0,
                                   complex=self.__complex)

    def refine_at_voxel(self, marked_voxels: Voxels) -> None:
        """Refine the mesh at the marked locations.

        Parameters
        ----------
        marked_locations : Voxels
            Representation of the locations which are to be refined:
            True if specific position has to be refined, False otherwise.
        """

        space = ngsolve.L2(self.__mesh, order=0)
        grid_function = ngsolve.GridFunction(space=space)
        values = marked_voxels.data
        cf = ngsolve.VoxelCoefficient(start=marked_voxels.start,
                                      end=marked_voxels.end,
                                      values=values.astype(float),
                                      linear=False)
        grid_function.Set(cf)
        flags = grid_function.vec.FV().NumPy()

        for element, flag in zip(self.__mesh.Elements(ngsolve.VOL), flags):
            self.__mesh.SetRefinementFlag(ei=element, refine=flag)
        self.refine()

    def refine_by_boundaries(self, boundaries: list) -> None:
        """Refine the mesh by the boundaries.

        Parameters
        ----------
        boundaries : list of str
            Collection of boundary names.
        """

        for element in self.__mesh.Elements(ngsolve.BND):
            to_refine = element.mat in boundaries
            self.__mesh.SetRefinementFlag(ei=element, refine=to_refine)
        self.refine()

    def refine_by_error(self, gridfunction: ngsolve.GridFunction) -> List:
        """Refine the mesh by the error at each mesh element.

        Parameters
        ----------
        gridfunction : ngsolve.GridFunction
        """

        flux = ngsolve.grad(gridfunction)
        space = self.mesh.flux_space()
        flux_potential = ngsolve.GridFunction(space=space)
        flux_potential.Set(coefficient=flux)
        difference = flux - flux_potential
        error = difference * ngsolve.Conj(difference)

        element_errors = ngsolve.Integrate(cf=error,
                                           mesh=self.__mesh,
                                           VOL_or_BND=ngsolve.VOL,
                                           element_wise=True).real
        limit = 0.5 * max(element_errors)
        for element in self.__mesh.Elements(ngsolve.BND):
            to_refine = element_errors[element.nr] > limit
            self.__mesh.SetRefinementFlag(ei=element, refine=to_refine)
        self.refine()

    def surfacel2_space(self, boundaries: List[str]) -> ngsolve.comp.SurfaceL2:
        """Return a number SurfaceL2 on the mesh.

        Returns
        -------
        ngsolve.comp.SurfaceL2
            SurfaceL2 space with minimum order of 1.
        """

        dirichlet = '|'.join(boundary for boundary in boundaries)
        return ngsolve.SurfaceL2(mesh=self.__mesh,
                                 order=max(1, self.__order - 1),
                                 dirichlet=dirichlet,
                                 complex=self.__complex)
