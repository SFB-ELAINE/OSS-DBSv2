
from ossdbs.point_analysis.time_results import TimeResult
from ossdbs.point_analysis.voltage_setting import VoltageSetter
from abc import ABC, abstractmethod


class SpectrumMode(ABC):

    SPACING_FACTOR = 1e4

    def __init__(self, voltage_setter: VoltageSetter) -> None:
        self._voltage_setting = voltage_setter

    @abstractmethod
    def compute(self, signal, volume_conductor, points) -> TimeResult:
        pass
