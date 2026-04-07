# Copyright 2023, 2024 Konstantin Butenko, Jan Philipp Payonk
# Copyright 2023, 2024 Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import logging

import netgen.occ
import numpy as np

from ossdbs.electrodes import ElectrodeModel

from .brain_geometry import BrainGeometry
from .contacts import Contact
from .encapsulation_layers import EncapsulationLayer

_logger = logging.getLogger(__name__)


class ModelGeometry:
    """The final model geometry.

    bounding_box : BoundingBox
        Localization in 3D space of the geometry

    electrodes : Electrodes
        Collection of electrode models

    TODO refactor
    """

    def __init__(
        self,
        brain: BrainGeometry,
        electrodes: list[ElectrodeModel],
    ) -> None:
        self._brain = brain
        self._electrodes = electrodes
        self._contacts = []
        self._encapsulation_layers = []

        self._shape = self._construct_shape(self._brain, self._electrodes)
        self._geometry = None  # built lazily after mesh sizes are set
        self.update_contact_areas()

    @property
    def geometry(self) -> netgen.occ.OCCGeometry:
        """Return netgen geometry of the model.

        The OCCGeometry is built lazily so that mesh size properties
        (maxh on faces/edges) are captured by the constructor.

        Returns
        -------
        netgen.occ.OCCGeometry
        """
        if self._geometry is None:
            try:
                self._geometry = netgen.occ.OCCGeometry(self._shape)
            except netgen.occ.OCCException:
                _logger.error(
                    "The geometry couldn't be constructed. "
                    "Tip: reduce the size of the brain geometry "
                    "or remove the encapsulation layer."
                )
                raise
        return self._geometry

    def _construct_shape(self, brain, electrodes):
        """Create the OCC shape of this brain model.

        Returns the raw OCC shape (not yet wrapped in OCCGeometry)
        so that mesh size properties can be set before OCCGeometry
        construction captures them.
        """
        brain_geo = brain.geometry
        for idx, electrode in enumerate(electrodes, start=1):
            electrode.index = idx
            self._create_electrode_contacts(electrode)
            if not np.isclose(electrode.encapsulation_thickness, 0.0):
                encapsulation = electrode.encapsulation_geometry(
                    electrode.encapsulation_thickness
                )
                # because encapsulation can extend outside brain_geo
                encapsulation = encapsulation * brain_geo
                # TODO naming of surfaces does not work reliably, hack below
                # check future versions of NGSolve to enable following line
                # encapsulation.bc("EncapsulationLayerSurface_{}".format(idx))
                encapsulation_layer_name = f"EncapsulationLayer_{idx}"
                encapsulation.mat(encapsulation_layer_name)
                brain_geo = netgen.occ.Glue([brain_geo, encapsulation])
                # hack to name surfaces properly
                brain_geo = brain_geo - electrode.geometry
                for face in brain_geo.faces:
                    if face.name is None:
                        face.name = f"EncapsulationLayerSurface_{idx}"
                self._encapsulation_layers.append(
                    EncapsulationLayer(name=encapsulation_layer_name)
                )
            else:
                brain_geo = brain_geo - electrode.geometry
            built_correctly = self.check_brain_geo(brain_geo, electrode)
            if not built_correctly:
                raise RuntimeError("Geometry could not be built.")

        # add brain surface
        brain_surfaces = brain.get_surface_names()
        surface_areas = brain.get_surface_areas()
        for surface in brain_surfaces:
            self._contacts.append(Contact(name=surface, area=surface_areas[surface]))
        return brain_geo

    def update_contact_areas(self) -> None:
        """Update contact areas."""
        for contact in self.contacts:
            area_set = False
            for surface in self._shape.faces:
                if contact.name == surface.name:
                    if area_set:
                        _logger.warning(f"Trying to set area twice for {contact.name}")
                    area = surface.mass
                    area_set = True
                    contact.area = area
            if not area_set:
                raise RuntimeError(f"Area for {contact.name} not set")

    def check_brain_geo(
        self, brain_geo: netgen.occ.Solid, electrode: ElectrodeModel
    ) -> bool:
        """Check if brain geo has all contacts."""
        face_names = [face.name for face in brain_geo.faces]
        correct_geo = True
        for contact_index in range(1, electrode.n_contacts + 1):
            name = self.get_contact_name(electrode.index, contact_index)
            if name not in face_names:
                correct_geo = False
                _logger.error(f"Face {name} is not in final geometry.")
        return correct_geo

    @property
    def electrodes(self) -> list[ElectrodeModel]:
        """Return collection of electrodes.

        Returns
        -------
        list[ElectrodeModel]
        """
        return self._electrodes

    @property
    def contacts(self) -> list[Contact]:
        """Return collection of contacts.

        Returns
        -------
        list[Contact]
        """
        return self._contacts

    def get_contact_index(self, contact_name: str) -> int:
        """Get index of contact by name."""
        for idx, contact in enumerate(self._contacts):
            if contact.name == contact_name:
                return idx
        return -1

    def get_encapsulation_layer_index(self, encapsulation_layer_name: str) -> int:
        """Get index of encapsulation layer by name."""
        for idx, encapsulation_layer in enumerate(self._encapsulation_layers):
            if encapsulation_layer.name == encapsulation_layer_name:
                return idx
        return -1

    # ruff: noqa: C901
    def update_contact(self, idx: int, settings: dict) -> None:
        """Overwrite contact properties."""
        if idx >= len(self._contacts):
            raise ValueError(
                f"You want to access contact {idx} "
                "which is not in the geometry. "
                "The highest possible index is "
                f"{len(self._contacts) - 1}."
            )
        contact = self._contacts[idx]
        for setting, value in settings.items():
            if setting == "Active":
                contact.active = value
            elif setting == "Current[A]":
                contact.current = value
            elif setting == "Floating":
                contact.floating = value
            elif setting == "Voltage[V]":
                contact.voltage = value
            elif setting == "SurfaceImpedance":
                if "Model" not in value:
                    raise ValueError("No surface impedance model provided.")
                if "Parameters" not in value:
                    raise ValueError("No surface impedance model parameters provided.")
                contact.surface_impedance_model = value["Model"]
                contact.surface_impedance_parameters = value["Parameters"]
            elif setting == "MaxMeshSize":
                contact.max_h = value
                self.set_face_mesh_sizes({contact.name: value})
            elif setting == "MaxMeshSizeEdge":
                contact.edge_max_h = value
                self.set_edge_mesh_sizes({contact.name: value})
            elif setting in ["Contact_ID", "Name"]:
                continue
            else:
                raise ValueError(
                    f"Tried to update contact with invalid setting {setting}"
                )

        return

    def update_encapsulation_layer(self, idx: int, settings: dict) -> None:
        """Overwrite encapsulation layer properties."""
        encapsulation_layer = self._encapsulation_layers[idx]
        for setting, value in settings.items():
            if setting == "Material":
                encapsulation_layer.material = value
            elif setting == "DielectricModel":
                encapsulation_layer.dielectric_model = value
            elif setting == "DielectricParameters":
                encapsulation_layer.dielectric_parameters = value
            elif setting == "MaxMeshSize":
                encapsulation_layer.max_h = value
                self.set_volume_mesh_sizes({encapsulation_layer.name: value})
            elif setting == "Thickness[mm]":
                continue
            else:
                raise ValueError(
                    f"Tried to update encapsulation layer with setting {setting}"
                )

        return

    @property
    def encapsulation_layers(self) -> list:
        """Return collection of active contacts and contacts of property
        floating.
        """
        return self._encapsulation_layers

    def _create_electrode_contacts(self, electrode: ElectrodeModel) -> None:
        """Add contacts from electrode to geometry.


        Notes
        -----
        The properties of the electrode are overriden.

        """
        new_boundary_names = {}
        for contact_index in range(1, electrode.n_contacts + 1):
            name = self.get_contact_name(electrode.index, contact_index)
            new_boundary_names[f"Contact_{contact_index}"] = name
            self._contacts.append(Contact(name=name))
        electrode.set_contact_names(new_boundary_names)

    def get_contact_name(self, electrode_index: int, contact_index: int) -> str:
        """Return contact name."""
        return f"E{electrode_index}C{contact_index}"

    def get_floating_mode(self):
        """Check if floating and if yes, which mode.

        Returns
        -------
        str or None
            ``"FloatingImpedance"`` if any floating contact carries a surface
            impedance model, ``"Floating"`` if there are floating contacts but
            none has an impedance model, or ``None`` if no contacts are floating.
        """
        has_floating = False
        has_impedance = False
        for contact in self.contacts:
            if contact.floating:
                has_floating = True
                if contact.surface_impedance_model is not None:
                    has_impedance = True
        if has_impedance:
            return "FloatingImpedance"
        if has_floating:
            return "Floating"
        return None

    def set_mesh_sizes(self, mesh_sizes: dict) -> None:
        """Set mesh sizes on edges, faces, and volumes."""
        if "Edges" in mesh_sizes:
            self.set_edge_mesh_sizes(mesh_sizes["Edges"])
        if "Faces" in mesh_sizes:
            self.set_face_mesh_sizes(mesh_sizes["Faces"])
        if "Volumes" in mesh_sizes:
            self.set_volume_mesh_sizes(mesh_sizes["Volumes"])

    def set_edge_mesh_sizes(self, mesh_sizes: dict) -> None:
        """Set mesh sizes on edges."""
        for edge in self._shape.edges:
            if edge.name in mesh_sizes:
                edge.maxh = mesh_sizes[edge.name]
        self._invalidate_geometry()

    def set_face_mesh_sizes(self, mesh_sizes) -> None:
        """Set mesh sizes on faces."""
        for face in self._shape.faces:
            if face.name in mesh_sizes:
                face.maxh = mesh_sizes[face.name]
        self._invalidate_geometry()

    def set_volume_mesh_sizes(self, mesh_sizes: dict) -> None:
        """Set mesh sizes on volumes."""
        for solid in self._shape.solids:
            if solid.name in mesh_sizes:
                solid.maxh = mesh_sizes[solid.name]
        self._invalidate_geometry()

    def _invalidate_geometry(self) -> None:
        """Reset cached OCCGeometry so it is rebuilt with updated mesh sizes."""
        self._geometry = None
