from typing import List
import netgen.occ
import netgen.meshing
import ngsolve
import numpy as np


class Mesh:
    """Class for interacting with the mesh for FEM.

    Parameters
    ----------
    geometry : netgen.occ.OCCGeometry

    order : int
        Order of mesh elements.

    complex_datatype : bool
            True for complex data type, False otherwise.
    """

    def __init__(self,
                 geometry: netgen.occ.OCCGeometry,
                 order: int) -> None:
        self._geometry = geometry
        self._order = order
        self._mesh = None

    def generate_mesh(self, meshing_parameters: dict) -> None:
        netgen_mp = self._meshing_parameters(meshing_parameters)
        self._mesh = ngsolve.Mesh(self.geometry.GenerateMesh(mp=netgen_mp))
        self._mesh.Curve(order=self.order)

    def load_mesh(self, filename: str) -> None:
        self._mesh = ngsolve.Mesh(filename=filename)
        self._mesh.ngmesh.SetGeometry(self._geometry)
        self._mesh.Curve(order=self.order)

    def _meshing_parameters(self, mesh_parameters: dict):
        mesh_type = mesh_parameters['Type']

        if mesh_type == "Custom":
            if "CustomParameters" in mesh_parameters:
                custom_parameters = mesh_parameters
                if not isinstance(custom_parameters, dict):
                    raise ValueError("CustomParameters have to passed as a dict.")
                return netgen.meshing.MeshingParameters(**custom_parameters)
            else:
                raise ValueError("You need to specific CustomParameters if you want to generate a custom mesh.")

        return {'Coarse': netgen.meshing.meshsize.coarse,
                'Fine': netgen.meshing.meshsize.fine,
                'Moderate': netgen.meshing.meshsize.moderate,
                'VeryCoarse': netgen.meshing.meshsize.very_coarse,
                'VeryFine': netgen.meshing.meshsize.very_fine,
                'Default': netgen.meshing.MeshingParameters(),
                }[mesh_type]

    @property
    def order(self) -> int:
        return self._order

    @order.setter
    def order(self, value: int) -> None:
        self._order = value

    @property
    def geometry(self) -> netgen.occ.OCCGeometry:
        return self._geometry

    @property
    def boundaries(self) -> List:
        return self.ngsolvemesh.GetBoundaries()

    @property
    def materials(self) -> List:
        return self.ngsolvemesh.GetMaterials()

    def boundary_coefficients(self, boundaries: dict) \
            -> ngsolve.fem.CoefficientFunction:
        """Return a boundary coefficient function.

        Returns
        -------
        ngsolve.fem.CoefficientFunction
        """

        return self._mesh.BoundaryCF(values=boundaries)

    def material_coefficients(self, materials: dict) \
            -> ngsolve.fem.CoefficientFunction:
        """Return a boundary coefficient function.

        Returns
        -------
        ngsolve.fem.CoefficientFunction
        """

        return self._mesh.MaterialCF(values=materials)

    @property
    def ngsolvemesh(self) -> ngsolve.Mesh:
        """Return mesh as a ngsolve object.

        Returns
        -------
        ngsolve.comp.Mesh
        """
        return self._mesh

    def is_included(self, points: np.ndarray) -> np.ndarray:
        """Check each point in collection for collision with geometry.
        True if point is included in geometry, false otherwise.

        Parameters
        ----------
        points: np.ndarray
            Array of point coordinates (x, y, z).

        Returns
        -------
        np.ndarray
            Array representing the state of collision for each point.
            True if point is included in geometry, False otherwise.
        """
        x, y, z = points.T
        mips = self._mesh(x, y, z)
        return np.array([mip[5] != -1 for mip in mips])

    def refine(self, at_surface: bool = False) -> None:
        """Refine the mesh."""

        self._mesh.Refine(mark_surface_elements=at_surface)
        self._mesh.Curve(order=self._order)

    def curve(self, order: int) -> None:
        self._order = order
        self._mesh.Curve(order=order)

    def save(self, file_name: str) -> None:
        """Save netgen mesh.

        Parameters
        ----------
        file_name : str
            File name of the mesh data.
        """

        self._mesh.ngmesh.Save(file_name)

    def refine_at_voxel(self,
                        start: tuple,
                        end: tuple,
                        data: np.ndarray) -> None:
        """Refine the mesh at the marked locations.

        Parameters
        ----------
        start : tuple
            Lower coordinates of voxel space.

        end : tuple
            Upper coordinates of voxel space.

        data : np.ndarray
            Voxelvalues.
        """

        space = ngsolve.L2(self._mesh, order=0)
        grid_function = ngsolve.GridFunction(space=space)
        cf = ngsolve.VoxelCoefficient(start=start,
                                      end=end,
                                      values=data.astype(float),
                                      linear=False)
        grid_function.Set(cf)
        flags = grid_function.vec.FV().NumPy()

        for element, flag in zip(self._mesh.Elements(ngsolve.VOL), flags):
            self._mesh.SetRefinementFlag(ei=element, refine=flag)
        self.refine()

    def refine_at_materials(self, materials: List[str]) -> None:
        """Refine the mesh by the boundaries.

        Parameters
        ----------
        boundaries : List[str]
            Collection of material names.

        """
        for element in self._mesh.Elements(ngsolve.VOL):
            to_refine = element.mat in materials
            self._mesh.SetRefinementFlag(ei=element, refine=to_refine)

    def refine_at_boundaries(self, boundaries: list) -> None:
        """Refine the mesh by the boundaries.

        Parameters
        ----------
        boundaries : list of str
            Collection of boundary names.
        """

        for element in self._mesh.Elements(ngsolve.VOL):
            self._mesh.SetRefinementFlag(ei=element, refine=False)
        for element in self._mesh.Elements(ngsolve.BND):
            to_refine = element.mat in boundaries
            self._mesh.SetRefinementFlag(ei=element, refine=to_refine)
        self.refine(at_surface=True)

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
                                           mesh=self._mesh,
                                           VOL_or_BND=ngsolve.VOL,
                                           element_wise=True).real
        limit = 0.5 * max(element_errors)
        for element in self._mesh.Elements(ngsolve.BND):
            to_refine = element_errors[element.nr] > limit
            self._mesh.SetRefinementFlag(ei=element, refine=to_refine)
        self.refine()
