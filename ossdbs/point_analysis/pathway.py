import logging
from dataclasses import dataclass
from typing import List

import h5py
import numpy as np

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

        axons: List["Pathway.Axon"]
            List that contains all axons within one population.

        """

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
                for axon in sorted(population.axons, key=lambda x: int(x.name[4:]))
            ]
        )
        self._location = np.full(n_points, "")
        self._coordinates = self._initialize_coordinates()

    def _create_axons(self, file: h5py.File, group: str) -> list:
        """Create axons based on the input from the .h5 file.

        Parameters
        ----------
        file: h5py.File
            Loaded .h5 file, which contains structural information.

        group: str
            Name of the group, which contains the axons.

        Retruns
        -------
        axons: list
            Retruns list of all axons within one group.
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

        return start, idx, status_list

    def filter_for_geometry(self, grid_pts: np.ma.MaskedArray) -> np.ndarray:
        """Checks if a point of an axon is outside the geometry.
        If thats the case, the whole axon will be marked,
        and the points are removed for further processing.

        Parameters
        ----------
        grid_pts: np.ma.MaskedArray
            The array contains 1 if the corresponding point is inside
            the geometry, 0 otherwise.

        Returns
        -------
        filtered_points: np.ndarray
            Returns filtered_points which are inside the geometry.
        """
        x, y, z = grid_pts.T
        lattice_mask = np.invert(grid_pts.mask)[:, 0]
        idx_axon = 0
        n_points = 0

        for population in self._populations:
            for axon in sorted(population.axons, key=lambda x: int(x.name[4:])):
                axon_length = axon.points.shape[0]
                for idx in range(axon_length):
                    if not lattice_mask[idx_axon + idx]:
                        axon.status = -1
                if axon.status != -1:
                    n_points += axon_length
                idx_axon += axon_length

        if n_points == 0:
            raise ValueError("No points inside the computational domain.")
        filtered_points = np.zeros((n_points, 3))
        idx_points = 0
        idx_grid = 0
        for population in self._populations:
            counter = 0
            n_axons = len(population.axons)
            for axon in sorted(population.axons, key=lambda x: int(x.name[4:])):
                axon_length = axon.points.shape[0]
                if axon.status == 0:
                    counter += 1
                    filtered_points[idx_points : idx_points + axon_length, 0] = x.data[
                        idx_grid : idx_grid + axon_length
                    ]
                    filtered_points[idx_points : idx_points + axon_length, 1] = y.data[
                        idx_grid : idx_grid + axon_length
                    ]
                    filtered_points[idx_points : idx_points + axon_length, 2] = z.data[
                        idx_grid : idx_grid + axon_length
                    ]
                    idx_points += axon_length
                idx_grid += axon_length
            _logger.info(f"Total axons in {population.name}: {n_axons}")
            _logger.info(f"Outside the domain: {n_axons - counter}")
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
                        if inside_csf[idx_axon + idx]:
                            axon.status = -2  # set status -2 for inside csf
                        if inside_encap[idx_axon + idx]:
                            axon.status = -1  # set status -1 for inside encap
                    idx_axon = idx_axon + axon_length
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
