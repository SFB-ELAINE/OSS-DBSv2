# Copyright 2023, 2024 Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import logging

import netgen.occ
import numpy as np
import numpy.linalg as npl
from scipy.spatial.transform import Rotation

from .bounding_box import BoundingBox

_logger = logging.getLogger(__name__)


class BrainGeometry:
    """CAD model of the brain geometry.

    Parameters
    ----------
    shape: str
        Choose between "Sphere", "Ellipsoid", "Box", and "Custom"
    bounding_box : BoundingBox
        Bounding box of the geometry in voxel space, can be none for custom geometry
    trafo_matrix: np.ndarray
        Matrix for affine transformation
    rotate_initial_geo: bool
        for spheres and ellipsoids, the initial geometry might have an edge
        on the surface that causes problems when merging with the electrode.
        Without lead axes the sphere is turned 180 degrees about z, otherwise
        its poles are turned 90 degrees about the seam direction.

    Notes
    -----
    The geometry will be initialized once.

    OpenCascade spheres have two poles (degenerate edges) and a seam
    meridian. Where a lead leaves the brain close to one of them, the
    boolean operations leave micrometre edges and sliver faces that Netgen
    fails to mesh. Once the lead axes are known (see `set_lead_axes`), the
    sphere is rotated so that its poles lie perpendicular to the exit points
    and its seam on the opposite side.
    A custom shape can be loaded from a STEP or BREP file using
    `import_geometry`.

    The brain has a single surface named `BrainSurface` and
    a single volume named `Brain`.
    An exception can be made by adding the geometry via the
    `set_geometry` method.


    Examples
    --------
    >>> from ossdbs.model_geometry import BrainGeometry, BoundingBox
    >>> bbox = BoundingBox((-20, -20, -20), (20, 20, 20))
    >>> brain = BrainGeometry("Sphere", bounding_box=bbox)

    """

    def __init__(
        self,
        shape: str,
        bounding_box: BoundingBox = None,
        trafo_matrix: np.ndarray = None,
        translation: np.ndarray = None,
        rotate_initial_geo: bool = False,
    ) -> None:
        self._bbox = bounding_box
        self._geometry = None
        self._shape = shape
        self._rotate_initial_geo = rotate_initial_geo
        self._lead_axes = []
        self._custom_geometry = False
        if trafo_matrix is None:
            self._trafo_matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        else:
            self._trafo_matrix = trafo_matrix
        if translation is None:
            self._translation = np.array([0, 0, 0])
        else:
            self._translation = translation
        self._affine_trafo = netgen.occ.gp_GTrsf(
            mat=self._trafo_matrix.ravel(), vec=self._translation
        )

    @property
    def brain_region(self):
        """Return box-shaped region covered by brain."""
        start, end = self.geometry.bounding_box
        return BoundingBox(start, end)

    def set_lead_axes(self, lead_axes: list[tuple]) -> None:
        """Orient a sphere or ellipsoid brain away from the leads.

        Parameters
        ----------
        lead_axes: list[tuple]
            (tip position, direction) of each lead in world coordinates.
        """
        if self._custom_geometry or self._shape not in ("Sphere", "Ellipsoid"):
            return
        self._lead_axes = [
            (np.asarray(position, float), np.asarray(direction, float))
            for position, direction in lead_axes
        ]
        self._geometry = None

    def get_surface_areas(self) -> dict:
        """Get areas of surfaces in geometry."""
        surface_areas = {}
        for face in self.geometry.faces:
            if face.name not in surface_areas:
                surface_areas[face.name] = face.mass
            else:
                surface_areas[face.name] += face.mass
        return surface_areas

    def get_surface_names(self) -> list[str]:
        """Get names of surfaces in geometry."""
        surface_list = []
        for face in self.geometry.faces:
            if face.name not in surface_list:
                surface_list.append(face.name)
        return surface_list

    @property
    def geometry(self) -> netgen.occ.Solid:
        """Return OCC shape of the brain geometry.

        Returns
        -------
        netgen.occ.Solid
        """
        if self._geometry is None:
            _logger.debug(f"Inititalize geo of shape {self._shape}")
            self._geometry = self._create_shape()
            self._geometry.bc("BrainSurface")
            self._geometry.mat("Brain")

        return self._geometry

    def _create_shape(self) -> netgen.occ.Solid:
        if self._shape == "Sphere":
            return self._create_sphere()
        elif self._shape == "Box":
            return self._create_box()
        elif self._shape == "Ellipsoid":
            return self._create_ellipsoid()
        elif self._shape == "Custom":
            raise RuntimeError(
                """Please load the custom geometry from a CAD file using the
                `import_geometry` or add it via the `set_geometry` method."""
            )
        else:
            raise NotImplementedError(
                f"""The shape {self._shape} is not implemented.
                Please choose among these shapes: Box, Sphere, Ellipsoid and Custom."""
            )

    def _create_ellipsoid(self) -> netgen.occ.Solid:
        half = np.subtract(self._bbox.end, self._bbox.start) / 2
        x, y, z = half
        transformator = netgen.occ.gp_GTrsf(mat=[x, 0, 0, 0, y, 0, 0, 0, z])
        center = netgen.occ.Pnt(1, 1, 1)
        sphere = netgen.occ.Sphere(c=center, r=1)
        sphere = self._orient_sphere(
            sphere,
            center,
            to_world=self._trafo_matrix @ np.diag(half),
            world_center=self._trafo_matrix @ np.add(self._bbox.start, half)
            + self._translation,
        )
        ellipsoid = transformator(sphere).Move(self._bbox.start)
        return self._affine_trafo(ellipsoid)

    def _create_sphere(self) -> netgen.occ.Solid:
        x, y, z = np.subtract(self._bbox.end, self._bbox.start) / 2
        center = np.add(self._bbox.start, (x, y, z))
        radius = np.min([x, y, z])
        sphere = netgen.occ.Sphere(c=netgen.occ.Pnt(center), r=radius)
        sphere = self._orient_sphere(
            sphere,
            netgen.occ.Pnt(center),
            to_world=radius * self._trafo_matrix,
            world_center=self._trafo_matrix @ center + self._translation,
        )
        return self._affine_trafo(sphere)

    def _orient_sphere(
        self,
        sphere: netgen.occ.Solid,
        center: netgen.occ.Pnt,
        to_world: np.ndarray,
        world_center: np.ndarray,
    ) -> netgen.occ.Solid:
        """Rotate the sphere about its centre, see `_sphere_rotation`.

        A unit vector u on the sphere maps to the world point
        ``to_world @ u + world_center``.
        """
        if not self._lead_axes:
            if self._rotate_initial_geo:
                _logger.debug("Rotate initial geometry")
                # we rotate 180 degrees around z, that should usually fix the issue
                sphere = sphere.Rotate(netgen.occ.Axis(center, (0, 0, 1)), 180)
            return sphere
        rotvec = self._sphere_rotation(to_world, world_center).as_rotvec()
        angle = npl.norm(rotvec)
        if np.isclose(angle, 0.0):
            return sphere
        axis = netgen.occ.Axis(center, tuple(rotvec / angle))
        return sphere.Rotate(axis, np.degrees(angle))

    def _sphere_rotation(
        self, to_world: np.ndarray, world_center: np.ndarray
    ) -> Rotation:
        """Rotation that moves the sphere poles and seam away from the leads.

        The exit points of the leads are mapped back to the unit sphere.
        The poles (local z) go perpendicular to them (least squares for
        several leads) and the seam meridian (local +x half plane) to the
        opposite side.
        """
        to_unit = npl.inv(to_world)
        exit_points = []
        for position, direction in self._lead_axes:
            tip = to_unit @ (position - world_center)
            axis = to_unit @ direction
            axis /= npl.norm(axis)
            # tip + t * axis on the unit sphere, t > 0
            b = tip @ axis
            c = tip @ tip - 1
            if c >= 0:
                _logger.warning("Lead tip outside of the brain geometry")
                continue
            exit_points.append(tip + (np.sqrt(b * b - c) - b) * axis)
        if not exit_points:
            return Rotation.identity()

        pole = npl.svd(np.array(exit_points))[2][-1]
        seam = -np.sum(exit_points, axis=0)
        seam -= (seam @ pole) * pole
        if np.isclose(npl.norm(seam), 0.0):
            seam = np.cross(pole, np.eye(3)[np.argmin(np.abs(pole))])
        seam /= npl.norm(seam)
        frame = np.column_stack([seam, np.cross(pole, seam), pole])
        if self._rotate_initial_geo:
            _logger.debug("Rotate initial geometry")
            # turn the poles 90 degrees about the seam direction
            frame = frame @ Rotation.from_euler("x", 90, degrees=True).as_matrix()
        return Rotation.from_matrix(frame)

    def _create_box(self) -> netgen.occ.Solid:
        if self._rotate_initial_geo:
            _logger.warning(
                "You chose to rotate the initial geometry."
                "This option is not available for a box geometry."
                "The geometry remains unchanged."
            )
        box = netgen.occ.Box(self._bbox.start, self._bbox.end)
        return self._affine_trafo(box)

    def import_geometry(self, path_to_geo_file: str):
        """Import brain geometry from a CAD file (STEP, IGES, or STL).

        The imported shape is assigned the boundary name ``BrainSurface``
        and the material name ``Brain``. The file is loaded via
        `Netgen's OCC interface <https://docu.ngsolve.org/latest/i-tutorials/unit-4a.2-csg2d/csg2d.html>`_.
        """
        _logger.debug(f"Import brain geometry from file: {path_to_geo_file}")
        occgeo = netgen.occ.OCCGeometry(path_to_geo_file)
        self._custom_geometry = True
        self._geometry = occgeo.shape
        self._geometry.bc("BrainSurface")
        self._geometry.mat("Brain")
        # get bounding box in voxel space
        bbox = self._geometry.bounding_box
        inv_trafo = npl.inv(self._trafo_matrix)
        inv_affine_trafo = netgen.occ.gp_GTrsf(
            mat=inv_trafo.ravel(), vec=-inv_trafo.dot(self._translation)
        )
        box = netgen.occ.Box(bbox[0], bbox[1])
        box_voxel = inv_affine_trafo(box)
        bbox_voxel = box_voxel.bounding_box
        self._bbox = BoundingBox(bbox_voxel[0], bbox_voxel[1])

    def set_geometry(self, geo: netgen.occ.Solid):
        """Set brain geometry from externally prepared OCC solid.

        Notes
        -----
        This function is the only possibility to use custom-defined
        names for surface parts of the brain!
        """
        self._custom_geometry = True
        self._geometry = geo
        bbox = self._geometry.bounding_box
        self._bbox = BoundingBox(bbox[0], bbox[1])
        self._geometry.mat("Brain")
        for face in self._geometry.faces:
            if face.name is None:
                face.name = "BrainSurface"
