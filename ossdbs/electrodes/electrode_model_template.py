from abc import ABC, abstractmethod
import netgen
import logging
import numpy as np
from dataclasses import dataclass

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

    def __init__(self,
                 parameters: dataclass,
                 rotation: float = 0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0),
                 ) -> None:
        self._position = position
        self._rotation = rotation
        norm = np.linalg.norm(direction)
        self._direction = tuple(direction / norm) if norm else (0, 0, 1)

        self._boundaries = {'Body': 'Body'}
        for idx in range(1, self._n_contacts + 1):
            self._boundaries['Contact_{}'.format(idx)] = 'Contact_{}'.format(idx)

        self._parameters = parameters
        self.parameter_check()
        self._geometry = self._construct_geometry()
        self._encapsulation_geometry = None
        self._encapsulation_thickness = 0.0
        self._index = 0

        pass

    @property
    def n_contacts(self) -> int:
        """Returns number of contacts.
        """
        return self._n_contacts

    @property
    def boundaries(self) -> dict:
        "Returns names of boundaries"
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
        return self._encapsulation_thickness

    @encapsulation_thickness.setter
    def encapsulation_thickness(self, thickness: float) -> None:
        self._encapsulation_geometry = self._construct_encapsulation_geometry(thickness)
        self._encapsulation_thickness = thickness

    def encapsulation_geometry(self, thickness: float) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:
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
            raise ValueError("The specified thickness is too small. Choose a larger, positive value.")
        if not np.isclose(thickness, self._encapsulation_thickness):
            return self._construct_encapsulation_geometry(thickness)
        return self._encapsulation_geometry

    @abstractmethod
    def parameter_check(self):
        pass
   
    @abstractmethod
    def _construct_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        pass

    @abstractmethod
    def _construct_encapsulation_geometry(self, thickness: float) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:
        pass

    def set_contact_names(self, boundaries: dict) -> None:
        """Set the names of electrode contacts.

        Parameters
        ----------
        contact_names : dict
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
    def index(self) -> int:
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
