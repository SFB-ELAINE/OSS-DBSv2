import logging
from typing import List

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
        electrodes: List[ElectrodeModel],
    ) -> None:
        self._brain = brain
        self._electrodes = electrodes
        self._contacts = []
        self._encapsulation_layers = []

        # TODO implement test for intersecting electrodes
        self._geometry = self._construct_geometry(self._brain, self._electrodes)

    @property
    def geometry(self) -> netgen.occ.OCCGeometry:
        """Return netgen geometry of the model.

        Returns
        -------
        netgen.occ.OCCGeometry
        """
        return self._geometry

    def _construct_geometry(self, brain, electrodes) -> netgen.occ.OCCGeometry:
        """Create a netgen geometry of this brain model.

        Returns
        -------
        netgen.occ.OCCGeometry
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

            brain_surfaces = brain.get_surface_names()
            for surface in brain_surfaces:
                self._contacts.append(Contact(name=surface))
        return netgen.occ.OCCGeometry(brain_geo)

    @property
    def electrodes(self) -> List[ElectrodeModel]:
        """Return collection of electrodes.

        Returns
        -------
        List[ElectrodeModel]
        """
        return self._electrodes

    @property
    def contacts(self) -> List[Contact]:
        """Return collection of contacts.

        Returns
        -------
        List[Contact]
        """
        return self._contacts

    def get_contact_index(self, contact_name: str) -> int:
        for idx, contact in enumerate(self._contacts):
            if contact.name == contact_name:
                return idx
        return -1

    def get_encapsulation_layer_index(self, encapsulation_layer_name: str) -> int:
        for idx, encapsulation_layer in enumerate(self._encapsulation_layers):
            if encapsulation_layer.name == encapsulation_layer_name:
                return idx
        return -1

    def update_contact(self, idx: int, settings: dict) -> None:
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
            elif setting == "SurfaceImpedance[Ohmm]":
                contact.surface_impedance = value
            elif setting == "MaxMeshSize":
                contact.max_h = value
                self.set_face_mesh_sizes({contact.name: value})
            elif setting == "MaxMeshSizeEdge":
                contact.max_h_edge = value
                self.set_edge_mesh_sizes({contact.name: value})
            elif setting == "Contact_ID":
                continue
            elif setting == "Name":
                continue
            else:
                raise ValueError(
                    f"Tried to update contact with invalid setting {setting}"
                )

        return

    def update_encapsulation_layer(self, idx: int, settings: dict) -> None:
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
    def encapsulation_layers(self) -> List:
        """Return collection of active contacts and contacts of property
        floating.

        Returns
        -------
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
            name = f"E{electrode.index}C{contact_index}"
            new_boundary_names[f"Contact_{contact_index}"] = name
            self._contacts.append(Contact(name=name))
        electrode.set_contact_names(new_boundary_names)

    def get_floating_mode(self):
        floating_mode = None
        for contact in self.contacts:
            if contact.floating:
                floating_mode = "Floating"
                if not np.isclose(contact.surface_impedance, 0.0j):
                    floating_mode = "FloatingImpedance"
                break
        if floating_mode == "Floating":
            for contact in self.contacts:
                if not np.isclose(contact.surface_impedance, 0.0j):
                    _logger.warning(
                        "Mode has been set to Floating but there is a nonzero surface impedance on contact {}".format(
                            contact.name
                        )
                    )

        return floating_mode

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
        for edge in self._geometry.shape.edges:
            if edge.name in mesh_sizes:
                edge.maxh = mesh_sizes[edge.name]

    def set_face_mesh_sizes(self, mesh_sizes) -> None:
        """Set mesh sizes on faces."""
        for face in self._geometry.shape.faces:
            if face.name in mesh_sizes:
                face.maxh = mesh_sizes[face.name]

    def set_volume_mesh_sizes(self, mesh_sizes: dict) -> None:
        """Set mesh sizes on volumes."""
        for solid in self._geometry.shape.solids:
            if solid.name in mesh_sizes:
                solid.maxh = mesh_sizes[solid.name]
