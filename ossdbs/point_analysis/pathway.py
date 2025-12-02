# Copyright 2023, 2024 Konstantin Butenko, Jan Philipp Payonk, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import logging
import os
from dataclasses import dataclass
from typing import Optional

import h5py
import numpy as np
import pandas as pd

from ossdbs.fem import Mesh
from ossdbs.utils.collapse_vta import get_collapsed_VTA
from ossdbs.utils.field_computation import (
    compute_field_magnitude_from_components,
)

from .lattice import PointModel
from .time_results import TimeResult

_logger = logging.getLogger(__name__)


class Pathway(PointModel):
    """Pathways comprise populations of axons."""

    @dataclass
    class Axon:
        """
        Attributes
        ----------
        name: str
            Naming of axons needs to be axon0, axon1, axon2, ... to be
            processed in the correct order.

        points: np.ndarray
            Contains 3D coordinates of each point within one axon.

        """

        name: str
        points: np.ndarray
        status: int  # 0 - normal, -1 - outside domain/encap, -2 - csf

    @dataclass
    class Population:
        """
        Attributes
        ----------
        name: str
            Name of neuronal population, e.g. a pathway.

        axons: list["Pathway.Axon"]
            List that contains all axons within one population.

        """

        name: str
        axons: list["Pathway.Axon"]

    def __init__(self, input_path: str, export_field: bool = False) -> None:
        # identifiers
        self._name = "PAM"
        self._export_field = export_field

        # path from where to read model
        self._path = input_path
        # never collapse VTA
        self.collapse_VTA = False
        # always compute time-domain signal
        self.time_domain_conversion = True

        with h5py.File(self._path, "r") as file:
            populations = [
                self.Population(group, self._create_axons(file, group))
                for group in file.keys()
            ]

        self._populations = populations
        n_points = sum(
            [
                len(axon.points)
                for population in self._populations
                for axon in sorted(population.axons, key=lambda x: int(x.name[4:]))
            ]
        )
        self._location = np.full(n_points, "")
        self._coordinates = self._initialize_coordinates()

        # will be set later
        self._lattice = None

    def _create_axons(self, file: h5py.File, group: str) -> list:
        """Create axons based on the input from the .h5 file.

        Parameters
        ----------
        file: h5py.File
            Loaded .h5 file, which contains structural information.

        group: str
            Name of the group, which contains the axons.

        Returns
        -------
        axons: list
            Returns list of all axons within one group.
        """
        return [
            self.Axon(sub_group, np.array(file[group][sub_group]), 0)
            for sub_group in file[group].keys()
        ]

    def _initialize_coordinates(self) -> np.ndarray:
        return np.concatenate(
            [
                axon.points
                for population in self._populations
                for axon in sorted(population.axons, key=lambda x: int(x.name[4:]))
            ]
        )

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
        Creates groups for each population in .h5 file.
        TODO Rename to 'create_data_export'or alike.

        """
        file.create_dataset("TimeSteps[s]", data=data.time_steps)
        start = 0
        idx = 0
        for population in self._populations:
            group = file.create_group(population.name)
            start, idx, status_list = self._create_datasets(
                data, start, idx, population, group
            )
            group.create_dataset("Status", data=status_list)

    def _create_datasets(self, data, start, idx, population, group):
        """Create datasets for each axon within the corresponding population.
        Axons are sorted numerically.
        """
        status_list = []
        for axon in sorted(population.axons, key=lambda x: int(x.name[4:])):
            sub_group = group.create_group(axon.name)
            sub_group.create_dataset("Points[mm]", data=axon.points)
            location = self._location[
                idx * len(axon.points) : (idx + 1) * len(axon.points)
            ]
            sub_group.create_dataset("Location", data=location.astype("S"))
            status_list.append(axon.status)
            if axon.status != -1:
                end = start + len(axon.points)
                # export potential
                potential = data.potential[start:end]
                sub_group.create_dataset("Potential[V]", data=potential)
                if data.electric_field_magnitude is not None:
                    # export field magnitude
                    electric_field_magnitude = data.electric_field_magnitude[start:end]
                    sub_group.create_dataset(
                        "Electric field magnitude[Vm^(-1)]",
                        data=electric_field_magnitude,
                    )
                if not (
                    data.electric_field_vector_x is None
                    and data.electric_field_vector_y is None
                    and data.electric_field_vector_z is None
                ):
                    # export field vector component-wise
                    electric_field_vector_x = data.electric_field_vector_x[start:end]
                    sub_group.create_dataset(
                        "Electric field vector x[Vm^(-1)]", data=electric_field_vector_x
                    )
                    electric_field_vector_y = data.electric_field_vector_y[start:end]
                    sub_group.create_dataset(
                        "Electric field vector y[Vm^(-1)]", data=electric_field_vector_y
                    )
                    electric_field_vector_z = data.electric_field_vector_z[start:end]
                    sub_group.create_dataset(
                        "Electric field vector z[Vm^(-1)]", data=electric_field_vector_z
                    )
                start = end
            idx = idx + 1

        return start, idx, status_list

    def filter_for_geometry(self, grid_pts: np.ma.MaskedArray) -> np.ndarray:
        """Check if any point of an axon is outside the geometry.
        If this is the case, the entire axon will be marked,
        and its points will be removed from further processing.

        Parameters
        ----------
        grid_pts: np.ma.MaskedArray
            Array containing points inside the mesh.

        Returns
        -------
        filtered_points: np.ndarray
            Returns filtered_points after removing axons that are (partially)
            outside the geometry.
        """
        x, y, z = grid_pts.T
        lattice_mask = np.invert(grid_pts.mask)[:, 0]
        idx_axon = 0

        total_points = sum(
            axon.points.shape[0]
            for population in self._populations
            for axon in sorted(population.axons, key=lambda x: int(x.name[4:]))
        )

        # Create an array of NaNs to filter the axons outside the domain
        all_points = np.full((total_points, 3), np.nan)

        point_idx = 0
        pop_axons_stats = []  # List of (population_name, n_axons, n_axons_inside)
        for population in self._populations:
            n_axons = len(population.axons)
            n_axons_inside = 0
            for axon in sorted(population.axons, key=lambda x: int(x.name[4:])):
                axon_length = axon.points.shape[0]
                axon_outside = False
                for idx in range(axon_length):
                    if not lattice_mask[idx_axon + idx]:
                        axon_outside = True
                        break
                if axon_outside:
                    axon.status = -1
                else:
                    # Fill all_points for this axon
                    all_points[point_idx : point_idx + axon_length, 0] = x.data[
                        idx_axon : idx_axon + axon_length
                    ]
                    all_points[point_idx : point_idx + axon_length, 1] = y.data[
                        idx_axon : idx_axon + axon_length
                    ]
                    all_points[point_idx : point_idx + axon_length, 2] = z.data[
                        idx_axon : idx_axon + axon_length
                    ]
                    point_idx += axon_length
                    n_axons_inside += 1
                idx_axon += axon_length
            pop_axons_stats.append((population.name, n_axons, n_axons_inside))

        filtered_points = all_points[~np.isnan(all_points).any(axis=1)]

        if np.isnan(filtered_points).any():
            raise RuntimeError(
                "NaN entries remain in filtered_points after filtering indicating a "
                "possible NumPy or logic bug."
            )

        if filtered_points.shape[0] == 0:
            raise ValueError("No points inside the computational domain.")

        for name, n_axons, n_axons_inside in pop_axons_stats:
            _logger.info(f"Total axons in {name}: {n_axons}")
            _logger.info(f"Outside the domain: {n_axons - n_axons_inside}")

        return filtered_points

    def filter_csf_encap(
        self, inside_csf: np.ndarray, inside_encap: np.ndarray
    ) -> None:
        """Change axon status if a single point of the axon is
        within the CSF or encapsulation layer.

        Parameters
        ----------
        inside_csf: np.ndarray
            The array contains 1 if the corresponding point is
            inside the CSF, 0 otherwise.

        inside_encap: np.ndarray
            The array contains 1 if the corresponding point is
            inside the encapsulation layer, 0 otherwise.
        """
        idx_axon = 0
        for population in self._populations:
            for axon in sorted(population.axons, key=lambda x: int(x.name[4:])):
                if axon.status != -1:
                    axon_length = axon.points.shape[0]
                    for idx in range(axon_length):
                        if inside_encap[idx_axon + idx]:
                            axon.status = -1  # set status -1 for inside encap
                            break
                        if inside_csf[idx_axon + idx]:
                            axon.status = -2  # set status -2 for inside csf
                            break
                    idx_axon += axon_length
        _logger.info("Marked axons inside CSF and encapsulation layer")
        return

    def create_index(self, lattice: np.ndarray) -> np.ndarray:
        """Create index for each point to the matching axon.

        Returns
        -------
        index: np.ndarray.
        """
        index = np.zeros(shape=len(lattice), dtype=int)
        axon_length = self.get_axon_length()
        for i in range(int(len(lattice) / axon_length)):
            index[i * axon_length : (i + 1) * axon_length] = int(i)
        return np.reshape(index, (len(index), 1))

    def get_axon_length(self) -> int:
        """
        Returns
        -------
        axon_length: int
            Number of points per axon

        Notes
        -----
        Assume the same length for all axons.
        """
        return self._populations[0].axons[0].points.shape[0]

    def get_population_names(self) -> list:
        """
        Returns
        -------
        population_names: list[str]
            Names of all populations defined
        """
        return [self._populations[idx].name for idx in range(len(self._populations))]

    def get_axon_names(self) -> list:
        """
        Returns
        -------
        axon_names: list[list,list,...]
            Names of axons in each population
        """
        axon_names = []
        for population in range(len(self._populations)):
            axon_names_in_population = []
            for axon in range(len(self._populations[population].axons)):
                axon_names_in_population.append(
                    self._populations[population].axons[axon].name
                )
            axon_names.append(axon_names_in_population)
        return axon_names

    def get_axon_numbers(self) -> list:
        """Get list of number of axons per population.

        Returns
        -------
        axon_number: list[int]
            Number of axons per population
        """
        return [
            len(self._populations[idx].axons) for idx in range(len(self._populations))
        ]

    def save_as_nifti(
        self, scalar_field, filename, binarize=False, activation_threshold=None
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
        raise NotImplementedError("Pathway results can not be stored in Nifti format.")

    def prepare_VCM_specific_evaluation(self, mesh: Mesh, conductivity_cf):
        """Prepare data structure according to mesh.

        Parameters
        ----------
        mesh: Mesh
            Mesh object on which VCM is defined
        conductivity_cf: ConductivityCF
            Conductivity function that holds material info

        Notes
        -----
        Mask all points outside domain, filter CSF and
        encapsulation layer etc.

        Prepares data storage for all frequencies at all points.
        """
        grid_pts = self.points_in_mesh(mesh)
        self._lattice_mask = np.invert(grid_pts.mask)
        self._lattice = self.filter_for_geometry(grid_pts)
        self._inside_csf = self.get_points_in_csf(mesh, conductivity_cf)
        self._inside_encap = self.get_points_in_encapsulation_layer(mesh)

        # mark complete axons and log how many axons were finally seeded
        self.filter_csf_encap(self.inside_csf, self.inside_encap)

        total_axons = sum(len(pop.axons) for pop in self._populations)
        seeded_axons = sum(
            sum(axon.status == 0 for axon in pop.axons) for pop in self._populations
        )
        _logger.info(f"Axons finally seeded: {seeded_axons} / {total_axons}")

        # create index for axons
        self._axon_index = self.create_index(self.lattice)

    def export_field_at_frequency(
        self,
        frequency: float,
        frequency_index: int,
        electrode=None,
        activation_threshold: Optional[float] = None,
    ):
        """Write field values to CSV.

        Parameters
        ----------
        frequency: float
            Frequency of exported solution
        frequency_index: int
            Index at which frequency is stored
        activation_threshold: float
            Threshold to define VTA
        electrode: ElectrodeModel
            electrode model that holds geometry information


        Notes
        -----
        No Nifti file is exported for a Pathway model.
        """
        if self.lattice is None:
            raise RuntimeError(
                "Please call first prepare_VCM_specific_evaluation "
                "to classify pathway points."
            )
        Ex = self.tmp_Ex_freq_domain[:, frequency_index]
        Ey = self.tmp_Ey_freq_domain[:, frequency_index]
        Ez = self.tmp_Ez_freq_domain[:, frequency_index]
        field_mags = compute_field_magnitude_from_components(Ex, Ey, Ez)
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
                [
                    self.lattice,
                    Ex.reshape((Ex.shape[0], 1)).real,
                    Ey.reshape((Ey.shape[0], 1)).real,
                    Ez.reshape((Ey.shape[0], 1)).real,
                    field_mags.reshape((field_mags.shape[0], 1)).real,
                ],
                axis=1,
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
        return
