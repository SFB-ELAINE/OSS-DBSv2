# Copyright 2023, 2024 Konstantin Butenko, Shruthi Chakravarthy
# Copyright 2023, 2024 Jan Philipp Payonk, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import logging
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass
from typing import Optional

import netgen
import netgen.occ as occ
import numpy as np
from ngsolve import BND, Mesh, VTKOutput

_logger = logging.getLogger(__name__)


class ElectrodeModel(ABC):
    """Deep Brain Simulation electrode.

    Attributes
    ----------
    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    position : tuple
        Position vector (x,y,z) of electrode tip.

    n_contacts: int
        Hard-coded number of electrode contacts.

    index: int
        Index of the electrode. Important for the model
        generation later and unambigous naming of boundaries.
    """

    _n_contacts: int

    def __init__(
        self,
        parameters: dataclass,
        rotation: float = 0,
        direction: tuple = (0, 0, 1),
        position: tuple = (0, 0, 0),
    ) -> None:
        self._position = position
        self._rotation = rotation
        norm = np.linalg.norm(direction)
        self._direction = tuple(direction / norm) if norm else (0, 0, 1)

        self._boundaries = {"Body": "Body"}
        for idx in range(1, self._n_contacts + 1):
            self._boundaries[f"Contact_{idx}"] = f"Contact_{idx}"

        self._parameters = parameters
        self.parameter_check()

        self._geometry = self._construct_geometry()
        self.check_initial_geometry()
        self._encapsulation_geometry = None
        self._encapsulation_thickness = 0.0
        self._index = 0

    def parameter_check(self):
        """Check electrode parameters."""
        # Check to ensure that all parameters are at least 0
        for param in asdict(self._parameters).values():
            if param < 0:
                raise ValueError("Parameter values cannot be less than zero")

    @property
    def n_contacts(self) -> int:
        """Returns number of contacts."""
        return self._n_contacts

    @property
    def boundaries(self) -> dict:
        """Returns names of boundaries."""
        return self._boundaries

    @property
    def geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Return geometry of electrode.

        Returns
        -------
        netgen.libngpy._NgOCC.TopoDS_Shape
        """
        return self._geometry

    @property
    def encapsulation_thickness(self) -> float:
        """Thickness of encapsulation layer."""
        return self._encapsulation_thickness

    @encapsulation_thickness.setter
    def encapsulation_thickness(self, thickness: float) -> None:
        if np.greater(thickness, 1e-3):
            self._encapsulation_geometry = self._construct_encapsulation_geometry(
                thickness
            )
        self._encapsulation_thickness = thickness

    def encapsulation_geometry(self, thickness: float) -> Optional[netgen.occ.Solid]:
        """Generate geometry of encapsulation layer around electrode.

        Parameters
        ----------
        thickness : float
            Thickness of encapsulation layer.

        Returns
        -------
        netgen.libngpy._NgOCC.TopoDS_Shape
        """
        if np.less(thickness, 1e-3):
            raise ValueError(
                "The specified thickness is too small. Choose a larger, positive value."
            )
        if not np.isclose(thickness, self._encapsulation_thickness):
            return self._construct_encapsulation_geometry(thickness)
        return self._encapsulation_geometry

    @abstractmethod
    def _construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        pass

    @abstractmethod
    def _construct_encapsulation_geometry(
        self, thickness: float
    ) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        pass

    def set_contact_names(self, boundaries: dict) -> None:
        """Set the names of electrode contacts.

        Parameters
        ----------
        boundaries : dict
            {'Body': 'body_name',
             'Contact_1': 'contact_name',
             'Contact_2': ...}
        """
        if self._boundaries == boundaries:
            _logger.info("Boundary names remain unchanged")
            return

        # TODO discuss if stricter checking required
        # currently, typos would not be catched, for example
        # checking that the keys are equivalent could help
        for face in self.geometry.faces:
            old_name = face.name
            if old_name in boundaries:
                face.name = boundaries[old_name]
        for edge in self.geometry.edges:
            old_name = edge.name
            if old_name in boundaries:
                edge.name = boundaries[old_name]

        self._boundaries.update(boundaries)
        _logger.info("Boundary names updated")

    @property
    def parameters(self) -> dataclass:
        """Electrode geometry parameters."""
        return self._parameters

    @property
    def index(self) -> int:
        """Index of electrode, relevant if multiple electrodes used."""
        return self._index

    @index.setter
    def index(self, index: int) -> None:
        self._index = index

    def get_max_mesh_size_contacts(self, ratio: float) -> float:
        """Use electrode's contact size to estimate maximal mesh size.

        Parameters
        ----------
        ratio: float
            Ratio between characteristic contact size and maximal mesh size.

        Notes
        -----
        For most of the electrodes, the electrode diameter is used.
        Exemptions are:
        * :class:`ossdbs.electrodes.MicroProbesSNEX100Model`

        """
        return self._parameters.lead_diameter / ratio

    def export_electrode(self, output_path, brain_dict, n_electrode) -> None:
        """Export electrode as Netgen and VTK file."""
        _logger.info("Export electrode as Netgen and VTK file")
        height = (
            np.amax(
                [
                    brain_dict["Dimension"]["x[mm]"],
                    brain_dict["Dimension"]["y[mm]"],
                    brain_dict["Dimension"]["z[mm]"],
                ]
            )
            / 2
        )

        radius = self._parameters.lead_diameter

        cylinder = netgen.occ.Cylinder(
            p=self._position,
            d=self._direction,
            r=radius,
            h=height,
        )

        occgeo = occ.OCCGeometry(cylinder * self.geometry)
        _logger.debug("Generating mesh")
        mesh_electrode = Mesh(occgeo.GenerateMesh())
        bnd_dict = {}
        for idx, contact in enumerate(self.boundaries):
            bnd_dict[contact] = idx
        boundary_cf = mesh_electrode.BoundaryCF(bnd_dict, default=-1)

        # export Netgen mesh
        mesh_electrode.ngmesh.Save(f"{output_path}/electrode_{n_electrode}.vol.gz")

        # export ParaView file
        VTKOutput(
            ma=mesh_electrode,
            coefs=[boundary_cf],
            names=["boundaries"],
            filename=f"{output_path}/electrode_{n_electrode}",
            subdivision=0,
        ).Do(vb=BND)

    def check_initial_geometry(self):
        """Check proper naming of electrode parts."""
        expected_names = [f"Contact_{i}" for i in range(1, self.n_contacts + 1)]
        # check faces
        face_names = []
        for face in self.geometry.faces:
            if face.name is None:
                continue
            if face.name not in face_names:
                face_names.append(face.name)
        # check edges
        edge_names = []
        for edge in self.geometry.edges:
            if edge.name is None:
                continue
            if edge.name not in edge_names:
                edge_names.append(edge.name)
        if not set(expected_names) == set(edge_names):
            raise RuntimeError("Edges have not been named correctly")
        expected_names.append("Body")
        if not set(expected_names) == set(face_names):
            raise RuntimeError("Faces have not been named correctly")

    def set_hp_flag(self, electrode_parameters: dict):
        """Set hp-flags only on active contacts."""
        if "Contacts" in electrode_parameters:
            for contact_info in electrode_parameters["Contacts"]:
                if contact_info["Active"]:
                    contact_idx = contact_info["Contact_ID"]
                    self._set_edge_hp_flag({f"Contact_{contact_idx}": 1})
                    self._set_vertex_hp_flag({f"Contact_{contact_idx}": 1})

    def _set_edge_hp_flag(self, edge_sizes: dict) -> None:
        """Set flags on edges."""
        for edge in self.geometry.edges:
            if edge.name in edge_sizes:
                edge.hpref = edge_sizes[edge.name]

    def _set_vertex_hp_flag(self, vertex_sizes) -> None:
        """Set flags on vertices."""
        for vertex in self.geometry.vertices:
            if vertex.name in vertex_sizes:
                vertex.hpref = vertex_sizes[vertex.name]
