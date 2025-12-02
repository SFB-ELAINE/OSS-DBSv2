# Copyright 2023, 2024 Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import logging
import os

import netgen.meshing
import netgen.occ
import ngsolve
import numpy as np

_logger = logging.getLogger(__name__)


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

    def __init__(self, geometry: netgen.occ.OCCGeometry, order: int) -> None:
        self._geometry = geometry
        self._order = order
        self._mesh = None

    def generate_mesh(self, meshing_parameters: dict) -> None:
        """Generate NGSolve mesh."""
        netgen_hypothesis = self.get_mesh_hypothesis(
            meshing_parameters["MeshingHypothesis"]
        )
        netgen_mp = self.get_meshing_parameters(meshing_parameters["MeshingHypothesis"])
        _logger.debug(f"Calling GenerateMesh with {netgen_hypothesis} and {netgen_mp}")
        self._mesh = ngsolve.Mesh(
            self.geometry.GenerateMesh(netgen_hypothesis, **netgen_mp)
        )
        if (
            "HPRefinement" in meshing_parameters
            and meshing_parameters["HPRefinement"]["Active"]
        ):
            _logger.info("Applying HP Refinement")
            self._mesh.RefineHP(
                levels=meshing_parameters["HPRefinement"]["Levels"],
                factor=meshing_parameters["HPRefinement"]["Factor"],
            )
        self._mesh.Curve(order=self.order)

    def load_mesh(self, filename: str) -> None:
        """Load NGSolve mesh from file."""
        if not os.path.isfile(filename):
            raise ValueError(
                f"Provide a correct filename to load the mesh,could not find {filename}"
            )
        self._mesh = ngsolve.Mesh(filename=filename)
        self._mesh.ngmesh.SetGeometry(self._geometry)
        self._mesh.Curve(order=self.order)

    def get_mesh_hypothesis(self, mesh_parameters: dict):
        """Get meshing hypothesis from Netgen/NGSolve."""
        mesh_type = mesh_parameters["Type"]

        if mesh_type == "Custom":
            if "CustomParameters" in mesh_parameters:
                custom_parameters = mesh_parameters
                if not isinstance(custom_parameters, dict):
                    raise ValueError("CustomParameters have to passed as a dict.")
                return netgen.meshing.MeshingParameters(**custom_parameters)
            else:
                raise ValueError(
                    """You need to specify CustomParameters
                    if you want to generate a custom mesh."""
                )

        return {
            "Coarse": netgen.meshing.meshsize.coarse,
            "Fine": netgen.meshing.meshsize.fine,
            "Moderate": netgen.meshing.meshsize.moderate,
            "VeryCoarse": netgen.meshing.meshsize.very_coarse,
            "VeryFine": netgen.meshing.meshsize.very_fine,
            "Default": netgen.meshing.MeshingParameters(),
        }[mesh_type]

    def get_meshing_parameters(self, mesh_parameters: dict):
        """Prepare NGSolve meshing parameters deviating from default."""
        meshing_hypothesis = {}
        if "MaxMeshSize" in mesh_parameters:
            meshing_hypothesis["maxh"] = mesh_parameters["MaxMeshSize"]
        if "CurvatureSafety" in mesh_parameters:
            meshing_hypothesis["curvaturesafety"] = mesh_parameters["CurvatureSafety"]
        if "Grading" in mesh_parameters:
            meshing_hypothesis["grading"] = mesh_parameters["Grading"]
        if "MeshSizeFilename" in mesh_parameters:
            meshing_hypothesis["meshsizefilename"] = mesh_parameters["MeshSizeFilename"]
        return meshing_hypothesis

    @property
    def order(self) -> int:
        """Order of curved mesh."""
        return self._order

    @order.setter
    def order(self, value: int) -> None:
        self._order = value

    @property
    def geometry(self) -> netgen.occ.OCCGeometry:
        """Underlying CAD geometry of mesh."""
        return self._geometry

    @property
    def boundaries(self) -> list:
        """Get list of boundary names."""
        return self.ngsolvemesh.GetBoundaries()

    @property
    def materials(self) -> list:
        """Get list of material names."""
        return self.ngsolvemesh.GetMaterials()

    def boundary_coefficients(
        self, boundaries: dict
    ) -> ngsolve.fem.CoefficientFunction:
        """Return a boundary coefficient function.

        Returns
        -------
        ngsolve.fem.CoefficientFunction
        """
        return self._mesh.BoundaryCF(boundaries)

    def material_coefficients(self, materials: dict) -> ngsolve.fem.CoefficientFunction:
        """Return a boundary coefficient function.

        Returns
        -------
        ngsolve.fem.CoefficientFunction
        """
        return self._mesh.MaterialCF(materials)

    @property
    def ngsolvemesh(self) -> ngsolve.Mesh:
        """Return mesh as a ngsolve object.

        Returns
        -------
        ngsolve.comp.Mesh
        """
        return self._mesh

    def not_included(self, points: np.ndarray) -> np.ndarray:
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
        return np.array([mip[5] == -1 for mip in mips])

    def refine(self, at_surface: bool = False) -> None:
        """Refine the mesh.

        Parameters
        ----------
        at_surface: bool
            Whether to mark surface elements.
        """
        self._mesh.Refine(mark_surface_elements=at_surface)
        self._mesh.Curve(order=self._order)

    def curve(self, order: int) -> None:
        """Curve mesh and overwrite mesh order.

        Parameters
        ----------
        order: int
            Order of curved mesh
        """
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

    def refine_at_voxel(self, start: tuple, end: tuple, data: np.ndarray) -> None:
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
        cf = ngsolve.VoxelCoefficient(
            start=start, end=end, values=data.astype(float), linear=False
        )
        grid_function.Set(cf)
        flags = grid_function.vec.FV().NumPy()

        for element, flag in zip(self._mesh.Elements(ngsolve.VOL), flags):
            self._mesh.SetRefinementFlag(ei=element, refine=flag)
        self.refine()

    def refine_at_materials(self, materials: list[str]) -> None:
        """Refine the mesh by the boundaries.

        Parameters
        ----------
        materials : list[str]
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

    def refine_by_error_cf(self, error_cf: ngsolve.GridFunction) -> list:
        """Refine the mesh by the error at each mesh element.

        Parameters
        ----------
        error_cf : ngsolve.GridFunction
            Estimated error from
        """
        elementwise_error = ngsolve.Integrate(
            cf=error_cf.real, mesh=self._mesh, VOL_or_BND=ngsolve.VOL, element_wise=True
        )
        max_error = max(elementwise_error)
        # convolved but here the Netgen (!) not NGsolve mesh is called
        self._mesh.ngmesh.Elements3D().NumPy()["refine"] = (
            elementwise_error.NumPy() > 0.25 * max_error
        )
        self.refine()

    def refine_by_material_cf(self, material_cf: ngsolve.GridFunction) -> list:
        """Refine the mesh by checking the material cf.
        Each element-wise integral divided by area should yield an integer value.
        If not, it needs to be refined.

        Parameters
        ----------
        material_cf : ngsolve.GridFunction
            CF with material indices (need to be integers)
        """
        element_idx = ngsolve.Integrate(
            cf=material_cf, mesh=self._mesh, VOL_or_BND=ngsolve.VOL, element_wise=True
        )
        element_area = ngsolve.Integrate(
            cf=ngsolve.CF(1.0),
            mesh=self._mesh,
            VOL_or_BND=ngsolve.VOL,
            element_wise=True,
        )

        # get the average element index per element
        element_idx = element_idx.NumPy() / element_area.NumPy()

        # check the modulo to get integers
        element_mods = np.mod(element_idx, 1)
        # modulos of integers are either close to 0
        # or to 1 (actually 0.9999...) due to rounding errors
        # if it is not integer, both are False
        to_refine = np.equal(
            np.isclose(element_mods, 0.0, atol=1e-4),
            np.isclose(element_mods, 1.0, atol=1e-4),
        )
        # label elements
        self._mesh.ngmesh.Elements3D().NumPy()["refine"] = to_refine
        # refine
        self.refine()

    @property
    def n_elements(self) -> int:
        """Number of elements."""
        return self.ngsolvemesh.ne
