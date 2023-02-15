
from ossdbs.impedance_analysis.fourier_analysis.spectrum \
    import SpectrumMode, Impedances
from ossdbs.stimulation_signal import Signal
from ossdbs.volume_conductor import VolumeConductor
import numpy as np
import ngsolve


class OctaveBandMode(SpectrumMode):

    class OctaveBand:

        SQRT2 = np.sqrt(2)

        def __init__(self, frequency: float) -> None:
            self.frequency = frequency

        def lower_limit(self):
            return self.frequency / self.SQRT2

        def upper_limit(self):
            return self.frequency * self.SQRT2

    def compute(self,
                signal: Signal,
                volume_conductor: VolumeConductor) -> Impedances:

        n_frequencies = (self.SPACING_FACTOR // 2) + 1
        frequencies = signal.frequency * np.arange(0, n_frequencies)
        n_octaves = int(np.log2(len(frequencies) - 1)) + 1
        octave_indices = 2 ** np.arange(0, n_octaves)
        octave_frequencies = signal.frequency * octave_indices
        octave_bands = [self.OctaveBand(freq) for freq in octave_frequencies]

        impedances = np.zeros((len(frequencies)))
        impedances[0] = self.__compute_impedance(volume_conductor, 0)

        for octave_band in octave_bands:
            start = int(octave_band.lower_limit() / signal.frequency + 1)
            end_index = int(octave_band.upper_limit() / signal.frequency + 1)
            end = min(end_index, len(frequencies))
            impedance = self.__compute_impedance(volume_conductor,
                                                 octave_band.frequency)
            impedances[start:end] = impedance

        return Impedances(frequency=frequencies, imdedance=impedances)

    @staticmethod
    def __compute_impedance(volume_conductor, frequency):
        solution = volume_conductor.compute_solution(frequency)
        field = ngsolve.grad(solution.potential)
        curr_dens_conj = ngsolve.Conj(solution.current_density)
        mesh = volume_conductor.mesh.ngsolvemesh()
        power = ngsolve.Integrate(field * curr_dens_conj, mesh)
        voltage = 1
        return voltage / power if power else 0
