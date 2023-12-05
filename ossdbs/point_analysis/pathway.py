import logging
import os
from dataclasses import dataclass
from typing import List

import h5py
import numpy as np

from .lattice import PointModel
from .time_results import TimeResult

_logger = logging.getLogger(__name__)


class Pathway(PointModel):
    @dataclass
    class Axon:
        name: str
        points: np.ndarray

    @dataclass
    class Population:
        name: str
        axons: List["Pathway.Axon"]

    def __init__(self, path) -> None:
        self._path = path
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
                for axon in population.axons
            ]
        )
        self._location = np.full(n_points, "")
        self._coordinates = self._initialize_coordinates()
        self._axon_mask = None

    def _create_axons(self, file, group) -> list:
        return [
            self.Axon(sub_group, np.array(file[group][sub_group]))
            for sub_group in file[group].keys()
        ]

    def _initialize_coordinates(self) -> np.ndarray:
        return np.concatenate(
            [
                axon.points
                for population in self._populations
                for axon in population.axons
            ]
        )

    def save(self, data: TimeResult, file_name: str):
        with h5py.File(file_name, "w") as file:
            self._write_file(data, file)

    def _write_file(self, data, file):
        file.create_dataset("TimeSteps[s]", data=data.time_steps)
        start = 0
        idx = 0
        for population in self._populations:
            group = file.create_group(population.name)
            # TODO use [0,-1,-2] instead of bool values for status
            group.create_dataset("Status", data=self._axon_mask[idx : idx + len(population.axons)])
            start, idx = self._create_datasets(data, start, idx, population, group)
           
    def _create_datasets(self, data, start, idx, population, group):
        # sort Axon instances numerically
        sorted_axons = sorted(population.axons, key=lambda axon_inst: int(axon_inst.name[4:]))
        for axon in sorted_axons:
            sub_group = group.create_group(axon.name)
            sub_group.create_dataset("Points[mm]", data=axon.points)
            location = self._location[idx * len(axon.points) : (idx + 1) * len(axon.points)]
            sub_group.create_dataset("Location", data=location.astype("S"))   
            if self._axon_mask[idx]:
                end = start + len(axon.points)
                potential = data.potential[start:end]
                sub_group.create_dataset("Potential[V]", data=potential)
                electric_field_magnitude = data.electric_field_magnitude[start:end]
                sub_group.create_dataset(
                    "Electric field magnitude[Vm^(-1)]", data=electric_field_magnitude
                )
                electric_field_vector_x = data.electric_field_vector[0][start:end]
                sub_group.create_dataset(
                    "Electric field vector x[Vm^(-1)]", data=electric_field_vector_x
                )
                electric_field_vector_y = data.electric_field_vector[1][start:end]
                sub_group.create_dataset(
                    "Electric field vector y[Vm^(-1)]", data=electric_field_vector_y
                )
                electric_field_vector_z = data.electric_field_vector[2][start:end]
                sub_group.create_dataset(
                    "Electric field vector z[Vm^(-1)]", data=electric_field_vector_z
                )
                start = end
            idx = idx + 1
        return start, idx

    def set_location_names(self, names: np.ndarray) -> None:
        self._location = names

    def filter_for_geometry(self, grid_pts: np.ma.MaskedArray) -> np.ndarray:
        """Return a lattice that NGSolve can process. If a single point from
        an axon is outside of the geometry, the whole axon will be removed.

        Parameters
        ----------
        grid_pts: np.ma.MaskedArray
        """
        axon_length = self._populations[0].axons[0].points.shape[0]

        if len(grid_pts) % axon_length != 0:
            raise RuntimeError(
                "The creation of the grid for the point analysis did not work"
            )

        total_number_axons = int(len(grid_pts) / axon_length)
        x, y, z = grid_pts.T

        lattice_mask = np.invert(grid_pts.mask)

        keep_axon = [True] * total_number_axons
        for i in range(0, total_number_axons):
            for j in range(axon_length):
                if not lattice_mask[i * axon_length + j][0]:
                    keep_axon[i] = False
                    j = axon_length - 1

        n_filtered_axons = 0
        for keep in keep_axon:
            if keep:
                n_filtered_axons = n_filtered_axons + 1

        _logger.info(f"Total amount of loaded axons: {total_number_axons}")
        _logger.info(
            f"Axons outside the computationl domain: {total_number_axons - n_filtered_axons}"
        )

        lattice = np.zeros(shape=(n_filtered_axons * axon_length, 3))

        idx = 0
        for i in range(total_number_axons):
            if keep_axon[i]:
                lattice[idx * axon_length: (idx + 1) * axon_length, 0] = x.data[
                    i * axon_length: (i + 1) * axon_length
                ]
                lattice[idx * axon_length: (idx + 1) * axon_length, 1] = y.data[
                    i * axon_length: (i + 1) * axon_length
                ]
                lattice[idx * axon_length: (idx + 1) * axon_length, 2] = z.data[
                    i * axon_length: (i + 1) * axon_length
                ]
                idx = idx + 1
        self._axon_mask = keep_axon
        return lattice

    def filter_csf_encap(
        self, inside_csf: np.ndarray, inside_encap: np.ndarray
    ) -> np.ndarray:
        """Marks whole axon in case a signle point of the axon is within the
        CSF or encapsulation layer.

        Parameters
        ----------
        inside_csf: np.ndarray

        inside_encap: np.ndarray
        """
        axon_length = self._populations[0].axons[0].points.shape[0]
        total_number_axons = int(len(inside_csf) / axon_length)

        axons_csf = [False] * len(inside_csf)
        axons_encap = [False] * len(inside_csf)

        for i in range(total_number_axons):
            for j in range(axon_length):
                if inside_csf[i * axon_length + j]:
                    axons_csf[i * axon_length: (i + 1) * axon_length] = [
                        True
                    ] * axon_length
                    j = axon_length - 1

        for i in range(total_number_axons):
            for j in range(axon_length):
                if inside_encap[i * axon_length + j]:
                    axons_encap[i * axon_length: (i + 1) * axon_length] = [
                        True
                    ] * axon_length
                    j = axon_length - 1
        axons_csf = np.asarray(axons_csf)
        axons_encap = np.asarray(axons_encap)

        return np.reshape(axons_csf, (len(axons_csf), 1)), np.reshape(
            axons_encap, (len(axons_encap), 1)
        )

    def save_hdf5(
        self,
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

        Notes
        -----
        TODO split in subfunctions
        """
        population_names = self.get_population_names()
        axon_names = self.get_axon_names()
        n_axons = self.get_axon_numbers()
        axon_length = self.get_axon_length()
        h5f_pts = h5py.File(os.path.join(output_path, "oss_pts.h5"), "w")
        idx = 0
        for i in range(len(population_names)):
            group = h5f_pts.create_group(population_names[i])
            for j in range(n_axons[i]):
                if self._axon_mask[sum(n_axons[:i]) + j]:
                    group.create_dataset(
                        axon_names[j],
                        data=lattice[idx * axon_length: (idx + 1) * axon_length, :],
                    )
                    idx = idx + 1
        h5f_pts.close()

        # Save potential evaluation
        h5f_pot = h5py.File(os.path.join(output_path, "oss_potentials.h5"), "w")
        idx = 0
        for i in range(len(population_names)):
            group = h5f_pot.create_group(population_names[i])
            for j in range(n_axons[i]):
                if self._axon_mask[sum(n_axons[:i]) + j]:
                    group.create_dataset(
                        axon_names[j],
                        data=lattice[idx * axon_length: (idx + 1) * axon_length, :],
                    )
                    group.create_dataset(
                        axon_names[j] + "_potentials",
                        data=potentials[idx * axon_length: (idx + 1) * axon_length, :],
                    )
                    idx = idx + 1
        h5f_pot.close()

        # Save electric field evaluation
        h5f_field = h5py.File(os.path.join(output_path, "oss_field.h5"), "w")
        idx = 0
        for i in range(len(population_names)):
            group = h5f_field.create_group(population_names[i])
            for j in range(n_axons[i]):
                if self._axon_mask[sum(n_axons[:i]) + j]:
                    group.create_dataset(
                        axon_names[j],
                        data=lattice[idx * axon_length: (idx + 1) * axon_length, :],
                    )
                    group.create_dataset(
                        axon_names[j] + "_field_vecs",
                        data=fields[idx * axon_length: (idx + 1) * axon_length, :],
                    )
                    group.create_dataset(
                        axon_names[j] + "_field_mags",
                        data=field_mags[idx * axon_length: (idx + 1) * axon_length, :],
                    )
                    idx = idx + 1
        h5f_field.close()

    def create_index(self, lattice: np.ndarray) -> np.ndarray:
        """Create index for each point to the matching axon.

        Returns
        -------
        index: np.ndarray.
        """
        index = np.zeros(shape=len(lattice), dtype=int)
        axon_length = self.get_axon_length()
        for i in range(int(len(lattice) / axon_length)):
            index[i * axon_length: (i + 1) * axon_length] = int(i)
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
        axon_names: list[str]
            Names of all axons defined
        """
        axon_names = []
        for population in range(len(self._populations)):
            for axon in range(len(self._populations[population].axons)):
                axon_names.append(self._populations[population].axons[axon].name)
        return axon_names

    def get_axon_numbers(self) -> list:
        """
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
        raise NotImplementedError("Pathway results can not be stored in Nifti format.")

    def collapse_VTA(self, field_on_points, implantation_coordinate, lead_direction, lead_diam):
        raise NotImplementedError("Collapse VTA for pathways not implemented")
