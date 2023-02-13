
from ossdbs.modes.spectrum import SpectrumMode, Result
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


class OctaveBandMode(SpectrumMode):

    def compute(self,
                signal: Signal,
                volume_conductor: VolumeConductor,
                points: np.ndarray):

        sample_spacing = 1 / (signal.frequency * self.SPACING_FACTOR)
        samples = signal.generate_samples(sample_spacing)
        complex_values = np.fft.rfft(samples)
        frequencies = np.fft.rfftfreq(len(samples), sample_spacing)
        n_octaves = int(np.log2(len(frequencies) - 1)) + 1
        octave_indices = 2 ** np.arange(0, n_octaves)
        octave_frequencies = signal.frequency * octave_indices
        octave_bands = [OctaveBand(freq) for freq in octave_frequencies]

        mesh = volume_conductor.mesh.ngsolvemesh()
        mips = [mesh(*point) for point in points]
        n_points = len(points)
        n_freq = len(frequencies)
        potentials = np.zeros((n_points, n_freq), dtype=complex)
        current_densities = np.zeros((n_points, n_freq, 3), dtype=complex)
        conductivities = np.zeros((n_points, n_freq), dtype=complex)

        solution = volume_conductor.compute_solution(0)
        potential = [solution.potential(mip) for mip in mips]
        current_density = [solution.current_density(mip) for mip in mips]
        conductivity = [solution.conductivity(mip) for mip in mips]
        potentials[:, 0] = np.array(potential) * complex_values[0]
        current_densities[:, 0] = np.array(current_density) * complex_values[0]
        conductivities[:, 0] = conductivity

        for octave_band in octave_bands:
            frequency = octave_band.frequency
            solution = volume_conductor.compute_solution(frequency)
            potential = [solution.potential(mip) for mip in mips]
            current_density = [solution.current_density(mip) for mip in mips]
            conductivity = [solution.conductivity(mip) for mip in mips]

            start = int(octave_band.lower_limit() / signal.frequency + 1)
            end = int(octave_band.upper_limit() / signal.frequency + 1)
            end = min(end, len(frequencies))

            for index in range(start, end):
                pointer = complex_values[index]
                potentials[:, index] = np.array(potential) * pointer
                current_densities[:, index] = np.array(current_density) * pointer
                conductivities[:, index] = conductivity

        return Result(points=points,
                      frequency=frequencies,
                      potential=potentials,
                      current_density=current_densities,
                      conductivity=conductivities)
