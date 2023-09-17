from dataclasses import dataclass
from typing import List

import h5py
import numpy as np

from .lattice import PointModel
from .time_results import TimeResult


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
        for population in self._populations:
            group = file.create_group(population.name)
            start = self._create_datasets(data, start, population, group)

    def _create_datasets(self, data, start, population, group):
        for axon in population.axons:
            end = start + len(axon.points)
            sub_group = group.create_group(axon.name)
            sub_group.create_dataset("Points[mm]", data=axon.points)
            location = self._location[start:end]
            sub_group.create_dataset("Location", data=location.astype("S"))
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
                "Electric field vector x[Vm^(-1)]", data=electric_field_vector_y
            )
            electric_field_vector_z = data.electric_field_vector[2][start:end]
            sub_group.create_dataset(
                "Electric field vector x[Vm^(-1)]", data=electric_field_vector_z
            )
            start = start + end
        return start

    def set_location_names(self, names: np.ndarray) -> None:
        self._location = names
