
from typing import List
from ossdbs.voxels import Voxels
import ngsolve
import re


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

    def get_boundaries(self) -> List:
        """Return all boundary names

        Returns
        -------
        list
            Collection of strings.
        """

        return list(self.__mesh.GetBoundaries())

    def get_floating_electrodes(self) -> List:
        """Return boundary names for floating electrodes.

        Returns
        -------
        list
            Collection of strings.

        Notes
        -----

        The naming has to be `floating_IDX` where
        `IDX` is an integer.

        """
        electrodes = list(set(self.get_boundaries()) - set(["default"]))
        floating_electrodes = []
        pattern = 'floating_[0-9]+'
        for e in electrodes:
            match = re.match(pattern, e)
            if match is None:
                matched = False
            else:
                matched = match.end() == match.endpos()
            if matched:
                floating_electrodes.append(e)
        return floating_electrodes

    def get_floating_impedance_electrodes(self) -> List:
        """Return boundary names for floating impedance electrodes.

        Returns
        -------
        list
            Collection of strings.

        Notes
        -----

        The naming has to be `floating_impedance_IDX` where
        `IDX` is an integer.
        """
        electrodes = list(set(self.get_boundaries()) - set(["default"]))
        floating_impedance_electrodes = []
        pattern = 'floating_impedance_[0-9]+'
        for e in electrodes:
            match = re.match(pattern, e)
            if match is None:
                matched = False
            else:
                matched = match.end() == match.endpos()
            if matched:
                floating_impedance_electrodes.append(e)
        return floating_impedance_electrodes

    def get_dirichlet_boundaries(self) -> List:
        """Return boundary names for Dirichlet BC.

        Returns
        -------
        list
            Collection of strings.
        """
        all_boundaries = self.get_boundaries()
        floating_electrodes = self.get_floating_electrodes()
        floating_impedance_electrodes = self.get_floating_impedance_electrodes()
        return list(set(all_boundaries) - set(['default']) - set(floating_electrodes) - set(floating_impedance_electrodes))

    def get_not_floating_boundaries(self) -> List:
        """Return boundary names that are not floating electrodes.

        Returns
        -------
        list
            Collection of strings.
        """

        return list(set(self.__mesh.GetBoundaries()) - set(['floating']))

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
                            order=self.__order - 1,
                            complex=self.__complex)

    def ngsolvemesh(self) -> ngsolve.comp.Mesh:
        """Return mesh as a ngsolve object.

        Returns
        -------
        ngsolve.comp.Mesh
        """

        return self.__mesh

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

    def is_complex(self) -> bool:
        """Check complex data type.

        Returns
        -------
        bool
            True if complex, False otherwise.
        """

        return self.__complex

    def h1_space(self) -> ngsolve.comp.H1:
        """Return a h1 space based on the mesh.

        Returns
        -------
        ngsolve.comp.H1
        """

        dirichlet = '|'.join(boundary for boundary in self.get_dirichlet_boundaries())
        return ngsolve.H1(mesh=self.__mesh,
                          order=self.__order,
                          dirichlet=dirichlet,
                          complex=self.__complex,
                          wb_withedges=False)

    def surfacel2_space(self) -> ngsolve.comp.SurfaceL2:
        """Return a surface l2 space based on the mesh.

        Returns
        -------
        ngsolve.comp.SurfaceL2

        Notes
        -----

        To implement a floating conductor, Lagrange multipliers
        have to be defined on the surfaces of the floating
        conductors.
        For that, a SurfaceL2 space is used. To keep only
        the DOFs on the conductor surface, all other surfaces
        have to be declared as `dirichlet` surfaces when
        constructing the function space.

        """

        not_floating = '|'.join(boundary for boundary in self.get_not_floating_boundaries())

        if self.__order > 1:
            order_lam = self.__order - 1
        else:
            order_lam = self.__order

        return ngsolve.SurfaceL2(mesh=self.__mesh,
                                 order=order_lam,
                                 dirichlet=not_floating,
                                 complex=self.__complex)

    def number_space(self) -> ngsolve.comp.NumberSpace:
        """Return a number space based on the mesh.

        Returns
        -------
        ngsolve.comp.NumberSpace

        Notes
        -----

        This space contains just one single (global) DOF.
        """

        return ngsolve.NumberSpace(mesh=self.__mesh,
                                   order=0,  # do not change!
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
        values = marked_voxels.data.astype(float)
        cf = ngsolve.VoxelCoefficient(start=marked_voxels.start,
                                      end=marked_voxels.end,
                                      values=values,
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
