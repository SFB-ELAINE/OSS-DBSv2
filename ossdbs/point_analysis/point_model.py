from abc import ABC, abstractmethod
from typing import Optional

import h5py
import numpy as np

from ossdbs.fem import Mesh
from ossdbs.point_analysis.time_results import TimeResult


class PointModel(ABC):
    """Class to support evaluation of VCM at points."""

    @property
    def coordinates(self) -> np.ndarray:
        """Point coordinates."""
        return self._coordinates

    @abstractmethod
    def _initialize_coordinates(self) -> np.ndarray:
        """Create grid / list of points."""
        pass

    def save(self, data: TimeResult, file_name: str) -> None:
        """Save time-domain result to HDF5 file."""
        with h5py.File(file_name, "w") as file:
            self._write_file(data, file)

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

    def filter_csf_encap(self, inside_csf: np.ndarray, inside_encap: np.ndarray):
        """Remove points in CSF or encapsulation layer.

        Parameters
        ----------
        inside_csf: np.ndarray
            list of points in csf
        inside_encap: np.ndarray
            list of points in encapsulation layer
        """
        raise NotImplementedError("Filtering of points not implemented.")

    @abstractmethod
    def save_as_nifti(
        self,
        scalar_field: np.ndarray,
        filename: str,
        binarize: bool = False,
        activation_threshold: Optional[float] = None,
    ):
        """Save scalar field in abstract orthogonal space in nifti format.

        Parameters
        ----------
        scalar_field : numpy.ndarray
            Nx1 array of scalar values on the lattice
        filename: str
            Name for the nifti file that should contain full path
        binarize: bool
            Choose to threshold the scalar field and save the binarized result
        activation_threshold: float
            Activation threshold for VTA estimate
        """
        pass

    def save_hdf5(
        self,
        axon_mask: list,
        lattice: np.ndarray,
        potentials: np.ndarray,
        fields: np.ndarray,
        field_mags: np.ndarray,
        output_path: str,
    ) -> None:
        """Export result to HDF5 format."""
        raise NotImplementedError("Results can not be stored in HDF5 format.")

    def _write_file(self, data: TimeResult, file: h5py.File):
        """Create datasets in HDF5 file.

        Parameters
        ----------
        data: TimeResult
            Time-domain result to be exported.
        file: h5py.File
            HDF5 file that shall contain data.
        """
        file.create_dataset("TimeSteps[s]", data=data.time_steps)
        file.create_dataset("Points[mm]", data=data.points)
        file.create_dataset("InsideCSF", data=data.inside_csf)
        file.create_dataset("InsideEncap", data=data.inside_encap)
        file.create_dataset("Location", data=self._location.astype("S"))
        file.create_dataset("Potential[V]", data=data.potential)
        file.create_dataset(
            "Electric field magnitude[Vm^(-1)]", data=data.electric_field_magnitude
        )
        file.create_dataset(
            "Electric field vector x[Vm^(-1)]", data=data.electric_field_vector[0]
        )
        file.create_dataset(
            "Electric field vector y[Vm^(-1)]", data=data.electric_field_vector[1]
        )
        file.create_dataset(
            "Electric field vector z[Vm^(-1)]", data=data.electric_field_vector[2]
        )
