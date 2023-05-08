from typing import List
from .contacts import Contact, Contacts
from .brain_geometry import BrainGeometry
from ossdbs.electrodes import ElectrodeModel
import netgen
import numpy as np


class ModelGeometry:
    """The final model geometry.

    bounding_box : BoundingBox
        Localization in 3D space of the geometry

    electrodes : Electrodes
        Collection of electrode models

    TODO refactor
    """
    def __init__(self,
                 brain: BrainGeometry,
                 electrodes: List[ElectrodeModel],
                 ) -> None:
        self.__brain = brain
        self.__electrodes = electrodes
        self.__contacts = []

        # TODO implement test for intersecting electrodes
        self.__geometry = self._construct_geometry(self.__brain, self.__electrodes)

    @property
    def geometry(self) -> netgen.occ.OCCGeometry:
        """Return netgen geometry of the model.

        Returns
        -------
        netgen.occ.OCCGeometry
        """
        return self.__geometry

    def _construct_geometry(self, brain, electrodes) -> netgen.occ.OCCGeometry:
        """Create a netgen geometry of this brain model.

        Returns
        -------
        netgen.occ.OCCGeometry
        """
        self.__contacts.append(Contact(name="BrainSurface"))
        brain_geo = brain.geometry
        for idx, electrode in enumerate(electrodes, start=1):
            electrode.index = idx
            self._create_electrode_contacts(electrode)
            if not np.isclose(electrode.encapsulation_thickness, 0.0):
                encapsulation = electrode.encapsulation_geometry(electrode.encapsulation_thickness)
                # because encapsulation can extend outside brain_geo
                encapsulation = encapsulation * brain_geo
                # TODO naming of surfaces does not work reliably, hack below
                # check future versions of NGSolve to enable following line
                # encapsulation.bc("EncapsulationLayerSurface_{}".format(idx))
                encapsulation.mat("EncapsulationLayer_{}".format(idx))
                brain_geo = netgen.occ.Glue([brain_geo, encapsulation])
                # hack to name surfaces properly
                brain_geo = brain_geo - electrode.geometry
                for face in brain_geo.faces:
                    if face.name is None:
                        face.name = "EncapsulationLayerSurface_{}".format(idx)
            else:
                brain_geo = brain_geo - electrode.geometry
        return netgen.occ.OCCGeometry(brain_geo)

    @property
    def contacts(self) -> Contacts:
        """Return collection of active contacts and contacts of property
        floating

        Returns
        -------
        Contacts
        """
        return self.__contacts

    def _create_electrode_contacts(self, electrode: ElectrodeModel) -> None:
        """Add contacts from electrode to geometry


        Notes
        -----

        The properties of the electrode are overriden.

        """
        new_boundary_names = {}
        for contact_index in range(1, electrode.n_contacts + 1):
            name = 'E{}C{}'.format(electrode.index, contact_index)
            new_boundary_names["Contact_{}".format(contact_index)] = name
            self.__contacts.append(Contact(name=name))
        electrode.set_contact_names(new_boundary_names)
