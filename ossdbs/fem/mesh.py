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
        self._hp_refinement_params = None
        self._hp_refinement_applied = False
        # Cache of located points, see locate_points().
        self._point_location_cache = {}

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
        self.invalidate_point_location_cache()
        if (
            "HPRefinement" in meshing_parameters
            and meshing_parameters["HPRefinement"]["Active"]
        ):
            # Store HP refinement parameters for deferred application.
            # HP refinement must be applied after any bisection-based
            # refinement (e.g. material refinement), because RefineHP
            # introduces element types that Netgen's bisection cannot handle.
            self._hp_refinement_params = meshing_parameters["HPRefinement"]
        self._mesh.Curve(order=self.order)

    def set_hp_refinement_params(self, hp_params: dict) -> None:
        """Store HP refinement parameters for deferred application.

        This is used when loading a pre-existing mesh so that
        apply_hp_refinement() can still be called later.
        Existing stored parameters are overwritten.
        """
        if not self._hp_refinement_applied:
            self._hp_refinement_params = hp_params

    @property
    def hp_refinement_applied(self) -> bool:
        """Whether HP refinement has been applied to this mesh."""
        return self._hp_refinement_applied

    def apply_hp_refinement(self) -> None:
        """Apply deferred HP refinement.

        Must be called after all bisection-based refinement steps
        (e.g. material refinement) are complete.
        """
        if self._hp_refinement_params is not None:
            _logger.info("Applying HP Refinement")
            self._mesh.RefineHP(
                levels=self._hp_refinement_params["Levels"],
                factor=self._hp_refinement_params["Factor"],
            )
            self._mesh.Curve(order=self._order)
            self._hp_refinement_params = None
            self._hp_refinement_applied = True
            self.invalidate_point_location_cache()

    def load_mesh(self, filename: str) -> None:
        """Load NGSolve mesh from file."""
        if not os.path.isfile(filename):
            raise ValueError(
                f"Provide a correct filename to load the mesh,could not find {filename}"
            )
        self._mesh = ngsolve.Mesh(filename=filename)
        self._mesh.ngmesh.SetGeometry(self._geometry)
        self._mesh.Curve(order=self.order)
        self.invalidate_point_location_cache()

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
        # Map input parameter names to NGSolve parameter names
        param_mapping = {
            "MaxMeshSize": "maxh",
            "CurvatureSafety": "curvaturesafety",
            "Grading": "grading",
            "MeshSizeFilename": "meshsizefilename",
        }
        return {
            ngsolve_key: mesh_parameters[input_key]
            for input_key, ngsolve_key in param_mapping.items()
            if input_key in mesh_parameters
        }

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

    def _mesh_token(self) -> tuple:
        """Return a value that changes whenever the mesh topology changes.

        Used to validate cached point locations. Any refinement adds
        elements and vertices, so a token built from those counts cannot
        miss a change.
        """
        return (id(self._mesh), self._mesh.ne, self._mesh.nv)

    def invalidate_point_location_cache(self) -> None:
        """Drop cached point locations, e.g. after the mesh has changed."""
        self._point_location_cache.clear()

    def locate_points(self, points: np.ndarray):
        """Return the points located on the mesh, cached.

        Locating points is a spatial search over the whole mesh and costs
        orders of magnitude more than evaluating a CoefficientFunction at
        the result. The same lattice is located for the CSF and
        encapsulation masks and again for the potential and the field at
        every frequency, so the result is cached and reused.

        The cache is keyed on the identity of the point array -- the array
        itself is stored alongside the result, which keeps it alive and
        therefore keeps its ``id`` from being reused by another object.
        Entries are additionally validated against the mesh token, so a
        refinement invalidates them even if nothing calls
        ``invalidate_point_location_cache``.

        Parameters
        ----------
        points: np.ndarray
            Nx3 array of point coordinates (x, y, z).
        """
        token = self._mesh_token()
        key = id(points)
        cached = self._point_location_cache.get(key)
        if cached is not None:
            cached_token, cached_points, mapping = cached
            if cached_token == token and cached_points is points:
                return mapping

        x, y, z = points.T
        mapping = self._mesh(x, y, z)
        self._point_location_cache[key] = (token, points, mapping)
        return mapping

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
        # Not cached: this is called once, on the unfiltered coordinates,
        # while the cached lattice is the filtered subset -- a different
        # array, so there would be nothing to reuse and the result would
        # only be retained for the lifetime of the mesh.
        x, y, z = points.T
        mips = self._mesh(x, y, z)
        # mips is a numpy structured array, "nr" is the element number
        # and -1 marks a point that could not be located in the mesh
        return mips["nr"] == -1

    def refine(self, at_surface: bool = False) -> None:
        """Refine the mesh.

        Parameters
        ----------
        at_surface: bool
            Whether to mark surface elements.
        """
        self._mesh.Refine(mark_surface_elements=at_surface)
        self._mesh.Curve(order=self._order)
        self.invalidate_point_location_cache()

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

    def get_boundary_areas(self) -> dict:
        """Integrate over boundaries to obtain surface areas."""
        surface_area_dict = {}
        for boundary in self.boundaries:
            area = ngsolve.Integrate(
                ngsolve.CF(1.0) * ngsolve.ds(boundary), mesh=self.ngsolvemesh
            )
            surface_area_dict[boundary] = area
        return surface_area_dict

    @property
    def n_elements(self) -> int:
        """Number of elements."""
        return self.ngsolvemesh.ne
