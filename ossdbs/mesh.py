
from typing import List
from ossdbs.voxels import Voxels
import ngsolve


class Mesh:
    """Class for interacting with the mesh for FEM.

    Parameters
    ----------
    geometry : netgen.libngpy._NgOCC.OCCGeometry

    order : int
        Order of mesh elements.
    """

    def __init__(self, mesh: ngsolve.comp.Mesh, order: int) -> None:
        self.__mesh = mesh
        self.__mesh.Curve(order=order)
        self.__order = order
        self.__complex = False

    def get_boundaries(self) -> List:
        """Return all boundary names.

        Returns
        -------
        list
            Collection of strings.
        """

        return list(set(self.__mesh.GetBoundaries()) - set(['default']))

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

    def set_complex(self, state: bool) -> None:
        """Set the data type to complex.

        Parameters
        ----------
        state : bool
            True for complex data type, False otherwise.
        """

        self.__complex = state

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

        dirichlet = '|'.join(boundary for boundary in self.get_boundaries())
        return ngsolve.H1(mesh=self.__mesh,
                          order=self.__order,
                          dirichlet=dirichlet,
                          complex=self.__complex,
                          wb_withedges=False)

    def refine_at_location(self, marked_locations: Voxels) -> None:
        """Refine the mesh at the marked locations.

        Parameters
        ----------
        marked_locations : Voxels
            Representation of the locations which are to be refined:
            True if specific position has to be refined, False otherwise.
        """

        space = ngsolve.L2(self.__mesh, order=0)
        grid_function = ngsolve.GridFunction(space=space)
        values = marked_locations.data.astype(float)
        cf = ngsolve.VoxelCoefficient(start=marked_locations.start,
                                      end=marked_locations.end,
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

        elements = self.__mesh.Elements(ngsolve.BND)
        flags = [element.mat in boundaries for element in elements]
        for element, flag in zip(self.__mesh.Elements(ngsolve.BND), flags):
            self.__mesh.SetRefinementFlag(ei=element, refine=flag)
        self.refine()

    def refine_by_error(self, error: ngsolve.fem.CoefficientFunction) -> List:
        """Refine the mesh by the error at each mesh element.

        Parameters
        ----------
        error : ngsolve.fem.CoefficientFunction
            Function holding all errors for each mesh element.
        """

        errors = ngsolve.Integrate(cf=error,
                                   mesh=self.__mesh,
                                   VOL_or_BND=ngsolve.VOL,
                                   element_wise=True).real
        limit = 0.5 * max(errors)
        flags = [errors[el.nr] > limit for el in self.__mesh.Elements()]
        for element, flag in zip(self.__mesh.Elements(ngsolve.BND), flags):
            self.__mesh.SetRefinementFlag(ei=element, refine=flag)
        self.refine()
