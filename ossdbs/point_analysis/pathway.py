from dataclasses import dataclass
from typing import List
import numpy as np
import logging
import h5py
import os

from .lattice import PointModel
from .time_results import TimeResult

_logger = logging.getLogger(__name__)


class Pathway(PointModel):
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
        status: int

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
        self._axon_mask = None



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
        for population in self._populations:
            print("population", population.name)
            for axon in sorted(population.axons, key=lambda x: int(x.name[4:])):
                if axon.name == "axon0":
                    print(np.array(axon.points).shape)
                # print(axon.name)
        return np.concatenate(
            [
                axon.points
                for population in self._populations
                for axon in sorted(population.axons, key=lambda x: int(x.name[4:]))
            ]
        )

    def save(self, data: TimeResult, file_name: str) -> None:   # creates "oss_time_result.h5"
        with h5py.File(file_name, "w") as file:
            self._write_file(data, file)

    def _write_file(self, data, file):
        file.create_dataset("TimeSteps[s]", data=data.time_steps)
        start = 0
        idx = 0
        for population in self._populations:
            group = file.create_group(population.name)
            # TODO use [0,-1,-2] instead of bool values for status
            #group.create_dataset("Status", data=self._axon_mask[idx : idx + len(population.axons)])
            start, idx = self._create_datasets(data, start, idx, population, group)
           
    def _create_datasets(self, data, start, idx, population, group):
        # sort Axon instances numerically

        for axon in sorted(population.axons, key=lambda x: int(x.name[4:])):
            sub_group = group.create_group(axon.name)
            sub_group.create_dataset("Points[mm]", data=axon.points)
            location = self._location[idx * len(axon.points) : (idx + 1) * len(axon.points)]
            sub_group.create_dataset("Location", data=location.astype("S"))   
            if axon.status != 1:
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


    def filter_for_geometry_new(self, grid_pts: np.ma.MaskedArray) -> np.ndarray:
        x, y, z = grid_pts.T
        lattice_mask = np.invert(grid_pts.mask)[:,0]
        temp = 0
        n_points = 0

        for population in self._populations:
            for axon in sorted(population.axons, key=lambda x: int(x.name[4:])):
                axon_length = axon.points.shape[0]
                for idx in range(axon_length):
                    if lattice_mask[temp + idx] == False:
                        axon.status = -1
                if lattice_mask[temp + idx] == True:
                        n_points = n_points + axon_length
                temp = temp + axon_length

        lattice = np.zeros(shape=(n_points, 3))

        temp_lattice = 0
        temp_grid = 0
        for population in self._populations:
            counter = 0
            n_axons = len(population.axons)
            for axon in sorted(population.axons, key=lambda x: int(x.name[4:])):
                axon_length = axon.points.shape[0]
                if axon.status == 0:
                    counter = counter + 1 
                    lattice[temp_lattice : temp_lattice + axon_length, 0] = x.data[temp_grid : temp_grid + axon_length]
                    lattice[temp_lattice : temp_lattice + axon_length, 1] = y.data[temp_grid : temp_grid + axon_length]
                    lattice[temp_lattice : temp_lattice + axon_length, 2] = z.data[temp_grid : temp_grid + axon_length]
                    temp_lattice = temp_lattice + axon_length
                    
                temp_grid = temp_grid + axon_length
            _logger.info(f"Total axons in {population.name}: {n_axons}")
            _logger.info(f"Outside the domain: {n_axons - counter}")

        return lattice

    # filter_geo_for_poulation as new function
    def filter_for_geometry(self, grid_pts: np.ma.MaskedArray) -> np.ndarray:
        """Return a lattice that NGSolve can process. If a single point from
        an axon is outside of the geometry, the whole axon will be removed.

        Parameters
        ----------
        grid_pts: np.ma.MaskedArray
        """
        # estimate which axons should be preserved
        lattice_mask = np.invert(grid_pts.mask)[:,0]
        x, y, z = grid_pts.T
        print("Lattice mask", lattice_mask)
        keep_axon = []
        size = 0
        last_population_index = 0
        for population in self._populations:
            axon_length = population.axons[0].points.shape[0]
            number_of_axons = len(population.axons)
            inside = [True] * number_of_axons

            print("Axon length:", axon_length)
            print("Number of axons", number_of_axons)
            index = []
            # shift based on population index whereever lattice is used
            for i in range(number_of_axons):
                for j in range(axon_length):
                    index.append(i * axon_length + j)
                    if not lattice_mask[last_population_index + i * axon_length + j]:
                        inside[i] = False
                        j = axon_length - 1
            
            n_filtered_axons = 0
            for keep in inside:
                if keep:
                    n_filtered_axons = n_filtered_axons + 1

            keep_axon = keep_axon + inside
            size = size + n_filtered_axons * axon_length

            _logger.info(f"Total amount of loaded axons: {number_of_axons}")
            _logger.info(
                f"Axons outside the computationl domain: {number_of_axons - n_filtered_axons}"
            )
            last_population_index = max(index)
            print("Max index", )
        self._axon_mask = keep_axon
        # create new lattice with filtered points
        lattice = np.zeros(shape=(size, 3))

        idx = 0
        i = 0
        for population in self._populations:
            axon_length = population.axons[0].points.shape[0]
            for axon_number in range(len(population.axons)):
                if self._axon_mask[i]:
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
                i = i + 1

        print("keep axons", keep_axon)
        print("Lattice filter for geo (in pw)", lattice)
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
    # stores oss_pts.h5, oss_potential.h5, oss_field.h5
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
                        axon_names[i][j],
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
                        axon_names[i][j],
                        data=lattice[idx * axon_length: (idx + 1) * axon_length, :],
                    )
                    group.create_dataset(
                        axon_names[i][j] + "_potentials",
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
                        axon_names[i][j],
                        data=lattice[idx * axon_length: (idx + 1) * axon_length, :],
                    )
                    group.create_dataset(
                        axon_names[i][j] + "_field_vecs",
                        data=fields[idx * axon_length: (idx + 1) * axon_length, :],
                    )
                    group.create_dataset(
                        axon_names[i][j] + "_field_mags",
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
        axon_names: list[list,list,...]
            Names of axons in each population
        """
        axon_names = []
        for population in range(len(self._populations)):
            axon_names_in_population = []
            for axon in range(len(self._populations[population].axons)):
                axon_names_in_population.append(self._populations[population].axons[axon].name)
            axon_names.append(axon_names_in_population)
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
