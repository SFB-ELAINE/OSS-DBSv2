import logging
import os
from abc import ABC, abstractmethod
from typing import Optional, Tuple

import h5py
import ngsolve
import numpy as np
import pandas as pd

from ossdbs.fem import Mesh
from ossdbs.point_analysis.time_results import TimeResult
from ossdbs.stimulation_signals import reconstruct_time_signals
from ossdbs.utils.collapse_vta import get_collapsed_VTA
from ossdbs.utils.field_computation import compute_field_magnitude

_logger = logging.getLogger(__name__)


class PointModel(ABC):
    """Class to support evaluation of VCM at points."""

    @property
    def collapse_VTA(self) -> bool:
        """Remove electrode from VTA."""
        return self._collapse_VTA

    @collapse_VTA.setter
    def collapse_VTA(self, value: bool):
        """Remove electrode from VTA."""
        if not isinstance(value, bool):
            raise ValueError("Provide a boolean value for VTA collapse")
        self._collapse_VTA = value

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

    def _write_file(self, data: TimeResult, file: h5py.File):
        """Create datasets in HDF5 file.

        Parameters
        ----------
        data: TimeResult
            Time-domain result to be exported.
        file: h5py.File
            HDF5 file that shall contain data.

        Notes
        -----
        TODO documentation
        """
        file.create_dataset("TimeSteps[s]", data=data.time_steps)
        file.create_dataset("Points[mm]", data=data.points)
        file.create_dataset("InsideCSF", data=data.inside_csf)
        file.create_dataset("InsideEncap", data=data.inside_encap)
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

    @property
    def lattice(self):
        """Lattice of points inside geometry."""
        return self._lattice

    @property
    def axon_index(self) -> np.ndarray:
        """Mapping between axons and points."""
        return self._axon_index

    @property
    def lattice_mask(self):
        """Mask of points inside geometry."""
        return self._lattice_mask

    @property
    def inside_csf(self):
        """Points inside csf."""
        return self._inside_csf

    @property
    def inside_encap(self):
        """Points inside encapsulation layer."""
        return self._inside_encap

    def prepare_VCM_specific_evaluation(self, mesh: Mesh, conductivity_cf):
        """Prepare data structure according to mesh.

        Parameters
        ----------
        mesh: Mesh
            Mesh object on which VCM is defined
        conductivity_cf: ConductivityCF
            Conductivity function that hold material info

        Notes
        -----
        Mask all points outside domain, filter CSF and
        encapsulation layer etc.

        Prepares data storage for all frequencies at all points.
        TODO type hints
        """
        grid_pts = self.points_in_mesh(mesh)
        self._lattice_mask = np.invert(grid_pts.mask)
        self._lattice = self.filter_for_geometry(grid_pts)
        self._axon_index = np.reshape(
            np.arange(len(self.lattice)), (len(self.lattice), 1)
        )
        self._inside_csf = self.get_points_in_csf(mesh, conductivity_cf)
        self._inside_encap = self.get_points_in_encapsulation_layer(mesh)

    def prepare_frequency_domain_data_structure(self, n_frequencies: int):
        """Allocate arrays that hold frequency domain data.

        Parameters
        ----------
        n_frequencies: int
            Number of frequencies of VCM
        """
        self.tmp_potential_freq_domain = np.zeros(
            shape=(n_frequencies, len(self.lattice)), dtype=complex
        )
        self.tmp_Ex_freq_domain = np.zeros(
            shape=(n_frequencies, len(self.lattice)), dtype=complex
        )
        self.tmp_Ey_freq_domain = np.zeros(
            shape=(n_frequencies, len(self.lattice)), dtype=complex
        )
        self.tmp_Ez_freq_domain = np.zeros(
            shape=(n_frequencies, len(self.lattice)), dtype=complex
        )

    def get_points_in_encapsulation_layer(self, mesh: Mesh) -> np.ndarray:
        """Return mask for points in encapsulation layer.

        Parameters
        ----------
        mesh: Mesh
            Mesh object on which VCM is defined
        """
        encap_cf = mesh.ngsolvemesh.RegionCF(
            ngsolve.VOL, {"EncapsulationLayer_*": 1.0}, default=0
        )
        ngmesh = mesh.ngsolvemesh
        x, y, z = self.lattice.T
        return np.isclose(encap_cf(ngmesh(x, y, z)), 1.0)

    def get_points_in_csf(self, mesh: Mesh, conductivity_cf) -> np.ndarray:
        """Return mask for points in CSF.

        Parameters
        ----------
        mesh: Mesh
            Mesh object on which VCM is defined
        conductivity_cf: ConductivityCF
            Conductivity function that hold material info

        Notes
        -----
        TODO Type hint
        """
        material_distribution = conductivity_cf.material_distribution(mesh)
        ngmesh = mesh.ngsolvemesh
        x, y, z = self.lattice.T
        return np.isclose(
            material_distribution(ngmesh(x, y, z)), conductivity_cf.materials["CSF"]
        )

    @property
    def output_path(self):
        """Path where to save results."""
        return self._output_path

    @output_path.setter
    def output_path(self, value):
        self._output_path = value

    def export_potential_at_frequency(
        self, frequency: float, frequency_index: float
    ) -> None:
        """Export potential at frequency to CSV.

        Parameters
        ----------
        frequency: float
            Frequency of exported solution
        frequency_index: int
            Index at which frequency is stored
        """
        potentials = self.tmp_potential_freq_domain[frequency_index]
        df_pot = pd.DataFrame(
            np.concatenate(
                [
                    self.axon_index,
                    self.lattice,
                    potentials.reshape((potentials.shape[0], 1)),
                    self.inside_csf,
                    self.inside_encap,
                ],
                axis=1,
                dtype=float,
            ),
            columns=[
                "index",
                "x-pt",
                "y-pt",
                "z-pt",
                "potential",
                "inside_csf",
                "inside_encap",
            ],
        )
        # save frequency
        df_pot["frequency"] = frequency
        df_pot.to_csv(os.path.join(self.output_path, "oss_potentials.csv"), index=False)

    def export_field_at_frequency(
        self,
        frequency: float,
        frequency_index: int,
        electrode=None,
        activation_threshold: Optional[float] = None,
    ):
        """Write field values to CSV and Nifti (if defined).

        Parameters
        ----------
        frequency: float
            Frequency of exported solution
        frequency_index: int
            Index at which frequency is stored
        activation_threshold: float
            Threshold to define VTA
        electrode:
            Electrode object with geometry details

        """
        Ex = self.tmp_Ex_freq_domain[frequency_index]
        Ey = self.tmp_Ey_freq_domain[frequency_index]
        Ez = self.tmp_Ez_freq_domain[frequency_index]
        fields = np.column_stack((Ex, Ey, Ez))
        field_mags = compute_field_magnitude(fields)
        df_field = pd.DataFrame(
            np.concatenate(
                [
                    self.axon_index,
                    self.lattice,
                    fields,
                    field_mags,
                    self.inside_csf,
                    self.inside_encap,
                ],
                axis=1,
                dtype=float,
            ),
            columns=[
                "index",
                "x-pt",
                "y-pt",
                "z-pt",
                "x-field",
                "y-field",
                "z-field",
                "magnitude",
                "inside_csf",
                "inside_encap",
            ],
        )
        # save frequency
        df_field["frequency"] = frequency

        if self.collapse_VTA:
            _logger.info("Collapse VTA by virtually removing the electrode")
            field_on_probed_points = np.concatenate(
                [self.lattice, fields, field_mags], axis=1, dtype=float
            )

            if electrode is None:
                raise ValueError(
                    "Electrode for exporting the collapsed VTA is missing."
                )
            implantation_coordinate = electrode._position
            lead_direction = electrode._direction
            lead_diam = electrode._parameters.lead_diameter

            field_on_probed_points_collapsed = get_collapsed_VTA(
                field_on_probed_points,
                implantation_coordinate,
                lead_direction,
                lead_diam,
            )

            df_collapsed_field = pd.DataFrame(
                np.concatenate(
                    [
                        self.axon_index,
                        field_on_probed_points_collapsed,
                        self.inside_csf,
                        self.inside_encap,
                    ],
                    axis=1,
                    dtype=float,
                ),
                columns=[
                    "index",
                    "x-pt",
                    "y-pt",
                    "z-pt",
                    "x-field",
                    "y-field",
                    "z-field",
                    "magnitude",
                    "inside_csf",
                    "inside_encap",
                ],
            )
            df_collapsed_field.to_csv(
                os.path.join(self.output_path, "E_field.csv"),
                index=False,
            )
        else:
            df_field.to_csv(
                os.path.join(self.output_path, "E_field.csv"),
                index=False,
            )

        # nifti exports
        field_mags_full = np.zeros(self.lattice_mask.shape[0])
        field_mags_full[self.lattice_mask[:, 0]] = np.real(field_mags[:, 0]) * 1000.0

        self.save_as_nifti(
            field_mags_full, os.path.join(self.output_path, "E_field_solution.nii")
        )
        self.save_as_nifti(
            field_mags_full,
            os.path.join(self.output_path, "VTA_solution.nii"),
            binarize=True,
            activation_threshold=activation_threshold,
        )

    def create_time_result(
        self,
        timesteps: np.ndarray,
        potential_in_time,
        field_in_time,
    ) -> TimeResult:
        """Prepare time result and save it to file.

        Notes
        -----
        TODO rethink structure for out-of-core processing.
        """
        time_result = TimeResult(
            time_steps=timesteps,
            points=self.lattice,
            potential=potential_in_time,
            electric_field_vector=field_in_time,
            electric_field_magnitude=compute_field_magnitude(field_in_time),
            inside_csf=self.inside_csf,
            inside_encap=self.inside_encap,
        )
        self.save(time_result, os.path.join(self.output_path, "oss_time_result.h5"))
        _logger.info("Created time result and saved to file")

    def compute_solutions_in_time_domain(
        self, signal_length: int
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute time-domain  solution for all proerties."""
        potential_in_time = reconstruct_time_signals(
            self.tmp_potential_freq_domain, signal_length
        )
        potential_in_time = np.swapaxes(potential_in_time, 0, 1)
        Ex_in_time = reconstruct_time_signals(self.tmp_Ex_freq_domain, signal_length)
        Ex_in_time = np.swapaxes(Ex_in_time, 0, 1)
        Ey_in_time = reconstruct_time_signals(self.tmp_Ey_freq_domain, signal_length)
        Ey_in_time = np.swapaxes(Ey_in_time, 0, 1)
        Ez_in_time = reconstruct_time_signals(self.tmp_Ez_freq_domain, signal_length)
        Ez_in_time = np.swapaxes(Ez_in_time, 0, 1)
        field_in_time = np.column_stack((Ex_in_time, Ey_in_time, Ez_in_time))
        return potential_in_time, field_in_time
