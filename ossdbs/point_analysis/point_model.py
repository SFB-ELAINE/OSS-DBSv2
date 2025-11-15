# Copyright 2023, 2024 Jan Philipp Payonk, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import logging
import os
from abc import ABC, abstractmethod
from typing import Optional

import h5py
import ngsolve
import numpy as np
import pandas as pd
from scipy.fft import ifft

from ossdbs.fem import Mesh
from ossdbs.point_analysis.time_results import TimeResult
from ossdbs.utils.collapse_vta import get_collapsed_VTA
from ossdbs.utils.field_computation import compute_field_magnitude_from_components

_logger = logging.getLogger(__name__)


class PointModel(ABC):
    """Class to support evaluation of VCM at points."""

    @property
    def name(self) -> str:
        """Name to distinguish model type."""
        return self._name

    @property
    def export_field(self) -> str:
        """Export electric field in time domain."""
        return self._export_field

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

    @property
    def time_domain_conversion(self) -> bool:
        """If time-domain signal shall be computed."""
        return self._time_domain_conversion

    @time_domain_conversion.setter
    def time_domain_conversion(self, value: bool):
        self._time_domain_conversion = value

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
        if data.electric_field_magnitude is not None:
            file.create_dataset(
                "Electric field magnitude[Vm^(-1)]", data=data.electric_field_magnitude
            )
        if not (
            data.electric_field_vector_x is None
            and data.electric_field_vector_y is None
            and data.electric_field_vector_z is None
        ):
            file.create_dataset(
                "Electric field vector x[Vm^(-1)]", data=data.electric_field_vector_x
            )
            file.create_dataset(
                "Electric field vector y[Vm^(-1)]", data=data.electric_field_vector_y
            )
            file.create_dataset(
                "Electric field vector z[Vm^(-1)]", data=data.electric_field_vector_z
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

    def prepare_frequency_domain_data_structure(
        self, signal_length: int, out_of_core: bool = False
    ) -> None:
        """Allocate arrays that hold frequency domain data.

        Parameters
        ----------
        signal_length: int
            Number of frequencies of FFT / input time vector
        out_of_core: bool
            If arrays should be stored on the hard drive
        """
        if out_of_core:
            # use unique file name
            self.tmp_hdf5_file = h5py.File(
                os.path.join(
                    self.output_path,
                    f"oss_freq_domain_tmp_{self.name}.hdf5",
                ),
                "w",
            )
            self.tmp_potential_freq_domain = self.tmp_hdf5_file.create_dataset(
                "Potential [V]",
                shape=(len(self.lattice), signal_length),
                dtype=complex,
                chunks=(len(self.lattice), 1),
            )
            self.tmp_Ex_freq_domain = self.tmp_hdf5_file.create_dataset(
                "Electric field vector x[Vm^(-1)]",
                shape=(len(self.lattice), signal_length),
                dtype=complex,
                chunks=(len(self.lattice), 1),
            )
            self.tmp_Ey_freq_domain = self.tmp_hdf5_file.create_dataset(
                "Electric field vector y[Vm^(-1)]",
                shape=(len(self.lattice), signal_length),
                dtype=complex,
                chunks=(len(self.lattice), 1),
            )
            self.tmp_Ez_freq_domain = self.tmp_hdf5_file.create_dataset(
                "Electric field vector z[Vm^(-1)]",
                shape=(len(self.lattice), signal_length),
                dtype=complex,
                chunks=(len(self.lattice), 1),
            )
        else:
            self.tmp_hdf5_file = None
            self.tmp_potential_freq_domain = np.zeros(
                shape=(len(self.lattice), signal_length), dtype=complex
            )
            self.tmp_Ex_freq_domain = np.zeros(
                shape=(len(self.lattice), signal_length), dtype=complex
            )
            self.tmp_Ey_freq_domain = np.zeros(
                shape=(len(self.lattice), signal_length), dtype=complex
            )
            self.tmp_Ez_freq_domain = np.zeros(
                shape=(len(self.lattice), signal_length), dtype=complex
            )

    def copy_frequency_domain_solution_from_vcm(
        self, freq_idx: int, potentials: np.ndarray, fields: Optional[np.ndarray] = None
    ) -> None:
        """Copy solution from volume conductor model."""
        signal_length = self.tmp_potential_freq_domain.shape[1]

        if signal_length % 2 == 0:
            # For even signals, highest frequency is not in positive frequencies
            if freq_idx > signal_length / 2 - 1:
                # copy potentials and fields
                self.tmp_potential_freq_domain[:, freq_idx] = np.conjugate(
                    potentials[:, 0]
                )
                if fields is not None:
                    self.tmp_Ex_freq_domain[:, freq_idx] = np.conjugate(fields[:, 0])
                    self.tmp_Ey_freq_domain[:, freq_idx] = np.conjugate(fields[:, 1])
                    self.tmp_Ez_freq_domain[:, freq_idx] = np.conjugate(fields[:, 2])
                return

        # copy potentials and fields
        self.tmp_potential_freq_domain[:, freq_idx] = potentials[:, 0]
        if fields is not None:
            self.tmp_Ex_freq_domain[:, freq_idx] = fields[:, 0]
            self.tmp_Ey_freq_domain[:, freq_idx] = fields[:, 1]
            self.tmp_Ez_freq_domain[:, freq_idx] = fields[:, 2]

        # if DC, there is no negative frequency
        if freq_idx == 0:
            return

        # get the negative frequency index
        negative_freq_idx = signal_length - freq_idx

        # Append the reverted signal without the DC frequency
        self.tmp_potential_freq_domain[:, negative_freq_idx] = np.conjugate(
            potentials[:, 0]
        )
        if fields is not None:
            self.tmp_Ex_freq_domain[:, negative_freq_idx] = np.conjugate(fields[:, 0])
            self.tmp_Ey_freq_domain[:, negative_freq_idx] = np.conjugate(fields[:, 1])
            self.tmp_Ez_freq_domain[:, negative_freq_idx] = np.conjugate(fields[:, 2])

    def close_output_file(self):
        """Close out-of-core file."""
        if self.tmp_hdf5_file is not None:
            self.tmp_hdf5_file.close()

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
        potentials = self.tmp_potential_freq_domain[:, frequency_index]
        df_pot = pd.DataFrame(
            np.concatenate(
                [
                    self.axon_index,
                    self.lattice,
                    potentials.reshape((potentials.shape[0], 1)).real,
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
        df_pot.to_csv(
            os.path.join(self.output_path, f"oss_potentials_{self.name}.csv"),
            index=False,
        )

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
        Ex = self.tmp_Ex_freq_domain[:, frequency_index]
        Ey = self.tmp_Ey_freq_domain[:, frequency_index]
        Ez = self.tmp_Ez_freq_domain[:, frequency_index]
        field_mags = compute_field_magnitude_from_components(Ex, Ey, Ez)
        # TODO find nicer solution
        df_field = pd.DataFrame(
            np.concatenate(
                [
                    self.axon_index,
                    self.lattice,
                    Ex.reshape((Ex.shape[0], 1)).real,
                    Ey.reshape((Ey.shape[0], 1)).real,
                    Ez.reshape((Ez.shape[0], 1)).real,
                    field_mags.reshape(field_mags.shape[0], 1),
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
            # TODO nicer solution
            field_on_probed_points = np.concatenate(
                [
                    self.lattice,
                    Ex.reshape((Ex.shape[0], 1)).real,
                    Ey.reshape((Ey.shape[0], 1)).real,
                    Ez.reshape((Ez.shape[0], 1)).real,
                    field_mags.reshape((field_mags.shape[0], 1)),
                ],
                axis=1,
                dtype=float,
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
                        field_on_probed_points_collapsed.real,
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
                os.path.join(self.output_path, f"E_field_{self.name}.csv"),
                index=False,
            )
        else:
            df_field.to_csv(
                os.path.join(self.output_path, f"E_field_{self.name}.csv"),
                index=False,
            )

        # nifti exports
        field_mags_full = np.zeros(self.lattice_mask.shape[0])
        # convert to V/m
        field_mags_full[self.lattice_mask[:, 0]] = field_mags * 1000.0

        self.save_as_nifti(
            field_mags_full,
            os.path.join(self.output_path, f"E_field_solution_{self.name}.nii"),
        )
        self.save_as_nifti(
            field_mags_full,
            os.path.join(self.output_path, f"VTA_solution_{self.name}.nii"),
            binarize=True,
            activation_threshold=activation_threshold,
        )

    def create_time_result(
        self,
        timesteps: np.ndarray,
        potential_in_time: np.ndarray,
        Ex_in_time: Optional[np.ndarray] = None,
        Ey_in_time: Optional[np.ndarray] = None,
        Ez_in_time: Optional[np.ndarray] = None,
        truncation_index: Optional[int] = None,
    ) -> TimeResult:
        """Prepare time result and save it to file.

        Parameters
        ----------
        timesteps: np.ndarray
            Array with timesteps related to solution
        potential_in_time: np.ndarray
            Solution array with electric potential
        Ex_in_time: np.ndarray
            Solution array with x-component of field
        Ey_in_time: np.ndarray
            Solution array with y-component of field
        Ez_in_time: np.ndarray
            Solution array with z-component of field
        truncation_index: int
            Index to truncate solution
        """
        _logger.info("Create time results")
        field_magnitude = None
        # if all field entries are defined
        if not (Ex_in_time is None and Ey_in_time is None and Ez_in_time is None):
            # truncate here otherwise they are none
            # first axis is points, second is time
            Ex_in_time = Ex_in_time[:, :truncation_index]
            Ey_in_time = Ey_in_time[:, :truncation_index]
            Ez_in_time = Ez_in_time[:, :truncation_index]
            field_magnitude = compute_field_magnitude_from_components(
                Ex_in_time, Ey_in_time, Ez_in_time
            )
        time_result = TimeResult(
            time_steps=timesteps[:truncation_index],
            points=self.lattice,
            # first axis is points, second is time
            potential=potential_in_time[:, :truncation_index],
            electric_field_vector_x=Ex_in_time,
            electric_field_vector_y=Ey_in_time,
            electric_field_vector_z=Ez_in_time,
            electric_field_magnitude=field_magnitude,
            inside_csf=self.inside_csf,
            inside_encap=self.inside_encap,
        )
        self.save(
            time_result,
            os.path.join(self.output_path, f"oss_time_result_{self.name}.h5"),
        )
        _logger.info("Created time result and saved to file")

    def compute_solutions_in_time_domain(
        self, signal_length: int, convert_field: bool = False
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Compute time-domain  solution for all properties."""
        out_of_core = self.tmp_hdf5_file is not None

        Ex_in_time = None
        Ey_in_time = None
        Ez_in_time = None

        """
        # TODO fix this part
        if out_of_core:
            # create dask array
            n_frequencies = self.tmp_potential_freq_domain.shape[1]
            n_points = self.tmp_potential_freq_domain.shape[0]
            array_size = self.tmp_potential_freq_domain.size
            # about 0.8 GB
            max_array_size = 5e7
            point_chunks = 1
            if array_size > max_array_size:
                max_points = int(max_array_size / n_frequencies)
                point_chunks = int(n_points / max_points)
            if point_chunks > 4:
                chunks = (point_chunks, n_frequencies)
                tmp_potential_freq_domain = da.from_array(
                    self.tmp_potential_freq_domain, chunks=chunks
                )
                potential_in_time = da.fft.ifft(tmp_potential_freq_domain, axis=1).real
                if convert_field:
                    _logger.warning(
                        "Out-of-core computation of field not yet implemented."
                    )
                return potential_in_time, Ex_in_time, Ey_in_time, Ez_in_time

            _logger.warning("Small array detected, will compute in core.")
        """

        if out_of_core:
            tmp_potential_freq_domain = self.tmp_potential_freq_domain[:]
        else:
            tmp_potential_freq_domain = self.tmp_potential_freq_domain

        potential_in_time = ifft(tmp_potential_freq_domain, axis=1, workers=-1).real
        if convert_field:
            Ex_in_time = ifft(self.tmp_Ex_freq_domain, axis=1, workers=-1).real
            Ey_in_time = ifft(self.tmp_Ey_freq_domain, axis=1, workers=-1).real
            Ez_in_time = ifft(self.tmp_Ez_freq_domain, axis=1, workers=-1).real
        return potential_in_time, Ex_in_time, Ey_in_time, Ez_in_time

    def export_point_model_information(self, filename: str) -> None:
        """Export all relevant information about the model to JSON."""
        raise NotImplementedError(
            "Point model information export has not yet been implemented."
        )

    def write_netgen_meshsize_file(self, meshsize: float, filename: str) -> None:
        """Use coordinates of point model to impose local mesh size.

        Notes
        -----
        Local mesh size for points is set.
        The file has the format (according to Netgen documentation):
          nr_points
          x1, y1, z1, meshsize
          x2, y2, z2, meshsize
          ...
          xn, yn, zn, meshsize

          nr_edges
          x11, y11, z11, x12, y12, z12, meshsize
          ...
          xn1, yn1, zn1, xn2, yn2, zn2, meshsize
        """
        points = self.coordinates
        with open(filename, "w") as fp:
            # write points to file
            fp.write(f"{len(points)}\n")
            fp.write("\n")
            for point in points:
                fp.write(f"{point[0]} {point[1]} {point[2]} {meshsize}\n")
            # we could also write lines but we do not
            fp.write("\n")
            fp.write("0\n")
