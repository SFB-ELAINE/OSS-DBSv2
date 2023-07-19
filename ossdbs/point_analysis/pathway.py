
from dataclasses import dataclass
from typing import List
import numpy as np
import h5py
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
        axons: List['Pathway.Axon']

    def __init__(self, path) -> None:
        self.__path = path
        with h5py.File(self.__path, "r") as file:
            populations = [self.Population(group,
                                           self.__create_axons(file, group))
                           for group in file.keys()]

        self.__populations = populations
        n_points = sum([len(axon.points)
                        for population in self.__populations
                        for axon in population.axons])
        self.__location = np.full(n_points, '')

    def __create_axons(self, file, group) -> list:
        return [self.Axon(sub_group, np.array(file[group][sub_group]))
                for sub_group in file[group].keys()]

    def coordinates(self) -> np.ndarray:
        return np.concatenate([axon.points
                               for population in self.__populations
                               for axon in population.axons])

    def save(self, data: TimeResult, file_name: str):
        with h5py.File(file_name, "w") as file:
            self.__write_file(data, file)

    def __write_file(self, data, file):
        file.create_dataset('TimeSteps[s]', data=data.time_steps)
        start = 0
        for population in self.__populations:
            group = file.create_group(population.name)
            start = self.__create_datasets(data, start, population, group)

    def __create_datasets(self, data, start, population, group):
        for axon in population.axons:
            end = start + len(axon.points)
            sub_group = group.create_group(axon.name)
            sub_group.create_dataset('Points[mm]', data=axon.points)
            location = self.__location[start:end]
            sub_group.create_dataset('Location', data=location.astype('S'))
            potential = data.potential[start:end]
            sub_group.create_dataset('Potential[V]', data=potential)
            electric_field_magnitude = data.electric_field_magnitude[start:end]
            sub_group.create_dataset('Electric field magnitude[Vm^(-1)]', data=electric_field_magnitude)
            electric_field_vector_x = data.electric_field_vector[0][start:end]
            sub_group.create_dataset('Electric field vector x[Vm^(-1)]', data=electric_field_vector_x)
            electric_field_vector_y = data.electric_field_vector[1][start:end]
            sub_group.create_dataset('Electric field vector x[Vm^(-1)]', data=electric_field_vector_y)
            electric_field_vector_z = data.electric_field_vector[2][start:end]
            sub_group.create_dataset('Electric field vector x[Vm^(-1)]', data=electric_field_vector_z)
            start = start + end
        return start

    def set_location_names(self, names: np.ndarray) -> None:
        self.__location = names
