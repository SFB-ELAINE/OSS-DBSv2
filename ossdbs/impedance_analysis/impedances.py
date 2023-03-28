from dataclasses import dataclass
import numpy as np
import pandas as pd


@dataclass
class Impedances:

    frequencies: np.ndarray
    imdedances: np.ndarray
    contact_sets: list

    def save(self, path: str) -> None:
        dataframe = pd.DataFrame(self.__create_data())
        dataframe.to_csv(path, index=False, sep=',')

    def __create_data(self) -> dict:
        data = {'frequencies [Hz]': self.frequencies}
        for index, contact_set in enumerate(self.contact_sets):
            name_1, name_2 = contact_set[0], contact_set[1]
            resistance = '_'.join(['resistance', name_1, name_2, '[Ohm]'])
            reactance = '_'.join(['reactance', name_1, name_2, '[Ohm]'])
            data.update({resistance: np.real(self.imdedances[:, index]),
                         reactance: np.imag(self.imdedances[:, index])})
        return data
