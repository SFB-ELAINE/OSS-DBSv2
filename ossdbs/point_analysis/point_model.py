from abc import ABC, abstractmethod
from typing import Optional

import numpy as np

from ossdbs.fem import Mesh
from ossdbs.point_analysis.time_results import TimeResult


class PointModel(ABC):
    @property
    def coordinates(self) -> np.ndarray:
        return self._coordinates

    @abstractmethod
    def _initialize_coordinates(self) -> np.ndarray:
        pass

    @abstractmethod
    def save(self, time_result: TimeResult, path: str) -> None:
        pass

    @abstractmethod
    def set_location_names(self, names: np.ndarray) -> None:
        pass

    def points_in_mesh(self, mesh: Mesh):
        """Create a masked array of the points that are in the mesh."""
        mask = mesh.not_included(self.coordinates)
        # use copy to avoid allocation of new memory
        return np.ma.masked_array(
            self.coordinates, np.column_stack((mask, mask, mask)), copy=True
        )

    def filter_for_geometry(self, grid_pts: np.ma.MaskedArray) -> np.ndarray:
        """Return a lattice that NGSolve can process.

        Notes
        -----
        The masked array is expected to be constructed by
        :meth:`points_in_mesh`
        """
        x, y, z = grid_pts.T
        x_compressed = np.ma.compressed(x)
        y_compressed = np.ma.compressed(y)
        z_compressed = np.ma.compressed(z)
        if not (len(x_compressed) == len(y_compressed) == len(z_compressed)):
            raise RuntimeError(
                "The creation of the grid for the point analysis did not work"
            )
        lattice = np.ndarray(shape=(len(x_compressed), 3))
        lattice[:, 0] = x_compressed
        lattice[:, 1] = y_compressed
        lattice[:, 2] = z_compressed
        return lattice

    @abstractmethod
    def filter_csf_encap(self, inside_csf: np.ndarray, inside_encap: np.ndarray):
        """Remove points in CSF or encapsulation layer.


        Parameters
        ----------
        inside_csf: list of points in csf
        inside_encap: list of points in encapsulation layer
        """
        pass

    @abstractmethod
    def save_as_nifti(
        self,
        scalar_field: np.ndarray,
        filename: str,
        binarize: bool = False,
        activation_threshold: Optional[float] = None,
    ):
        """Save scalar field in abstract orthogonal space using nifti format.

        Parameters
        ----------
        scalar_field : Nx1 numpy.ndarray of scalar values on the lattice
        filename: str, name for the nifti file that should contain full path
        binarize: bool, thresholds the scalar field and saves the binarized result
        activation_threshold: float, threshold for the binarized result

        """
        pass

    @abstractmethod
    def save_hdf5(
        self,
        axon_mask: list,
        lattice: np.ndarray,
        potentials: np.ndarray,
        fields: np.ndarray,
        field_mags: np.ndarray,
        output_path: str,
    ) -> None:
        """Stores results for pathway analysis in hdf5 format.

        Parameters
        ----------
        axon_mask: list
        lattice: np.ndarray
        potentials: np.ndarray
        fields: np.ndarray
        field_mags: np.ndarray
        output_path: str
        """
        pass
