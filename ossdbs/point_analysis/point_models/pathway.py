
from dataclasses import dataclass
from typing import List
import numpy as np
import h5py
from ossdbs.point_analysis.point_models.lattice import PointModel
from ossdbs.point_analysis.time_results import TimeResult


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
            populations = [self.Population(group, self.__axons(file, group))
                           for group in file.keys()]

        self.__populations = populations
        n_points = sum([len(axon.points)
                        for population in self.__populations
                        for axon in population])
        self.__location = np.full(n_points, '')

    def __axons(self, file, group) -> list:
        return [self.Axon(sub_group, np.array(file[group][sub_group]))
                for sub_group in file.keys[group]]

    def coordinates(self) -> np.ndarray:
        return np.concatenate([axon.points
                               for population in self.__populations
                               for axon in population.axons])

    def save(self, data: TimeResult, file_name: str):
        with h5py.File(file_name, "w") as file:
            self.__write_file(data, file)

    def __write_file(self, data, file):
        file.create_dataset('TimeSteps[s]', data=self.time_steps)
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
            current_density = data.current_density[start:end]
            sub_group.create_dataset('CurrentDensity[A|m2]',
                                     data=current_density)
            start = start + end
        return start
