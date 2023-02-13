
from ossdbs.impedance_mode.spectrum import SpectrumMode, ResultImpedance
import numpy as np
import ngsolve

from ossdbs.signals import Signal
from ossdbs.volume_conductor.volume_conductor_model import VolumeConductor


class OctaveBand:

    SQRT2 = np.sqrt(2)

    def __init__(self, frequency: float) -> None:
        self.frequency = frequency

    def lower_limit(self):
        return self.frequency / self.SQRT2

    def upper_limit(self):
        return self.frequency * self.SQRT2


class OctaveBandModeImpedance(SpectrumMode):

    def compute(self, signal: Signal, volume_conductor: VolumeConductor):

        n_frequencies = (self.SPACING_FACTOR // 2) + 1
        frequencies = signal.frequency * np.arange(0, n_frequencies)
        n_octaves = int(np.log2(len(frequencies) - 1)) + 1
        octave_indices = 2 ** np.arange(0, n_octaves)
        octave_frequencies = signal.frequency * octave_indices
        octave_bands = [OctaveBand(freq) for freq in octave_frequencies]

        mesh = volume_conductor.mesh.ngsolvemesh()
        impedances = np.zeros((len(frequencies)))

        solution = volume_conductor.compute_solution(0)

        field = ngsolve.grad(solution.potential)
        curr_dens_conj = ngsolve.Conj(solution.current_density)
        power = ngsolve.Integrate(field * curr_dens_conj, mesh)
        impedance = 1 / power if power else 0
        impedances[0] = impedance

        for octave_band in octave_bands:
            frequency = octave_band.frequency
            solution = volume_conductor.compute_solution(frequency)

            field = ngsolve.grad(solution.potential)
            curr_dens_conj = ngsolve.Conj(solution.current_density)
            power = ngsolve.Integrate(field * curr_dens_conj, mesh)
            impedance = 1 / power if power else 0

            start = int(octave_band.lower_limit() / signal.frequency + 1)
            end = int(octave_band.upper_limit() / signal.frequency + 1)
            end = min(end, len(frequencies))
            impedances[start:end] = impedance

        return ResultImpedance(frequency=frequencies, imdedance=impedances)
