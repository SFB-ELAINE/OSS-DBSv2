
from ossdbs.modes.spectrum import SpectrumMode, Result
import numpy as np
import ngsolve

from ossdbs.signals import Signal
from ossdbs.volume_conductor.volume_conductor_model import VolumeConductor


class NoTruncation(SpectrumMode):

    def compute(self,
                signal: Signal,
                volume_conductor: VolumeConductor,
                points: np.ndarray):

        sample_spacing = 1 / (signal.frequency * self.SPACING_FACTOR)
        samples = signal.generate_samples(sample_spacing)
        complex_values = np.fft.rfft(samples)
        frequencies = np.fft.rfftfreq(len(samples), sample_spacing)

        mesh = volume_conductor.mesh.ngsolvemesh()
        mips = [mesh(*point) for point in points]

        n_points = len(points)
        n_freq = len(frequencies)
        potentials = np.zeros((n_points, n_freq), dtype=complex)
        current_densities = np.zeros((n_points, n_freq, 3), dtype=complex)
        conductivities = np.zeros((n_points, n_freq), dtype=complex)

        for index, frequency in enumerate(frequencies[:1]):
            solution = volume_conductor.compute_solution(frequency)
            potential = [solution.potential(mip) for mip in mips]
            current_density = [solution.current_density(mip) for mip in mips]
            conductivity = [solution.conductivity(mip) for mip in mips]

            pointer = complex_values[index]
            potentials[:, index] = np.array(potential) * pointer
            current_densities[:, index] = np.array(current_density) * pointer
            conductivities[:, index] = conductivity

        return Result(points=points,
                      frequency=frequencies,
                      potential=potentials,
                      current_density=current_densities,
                      conductivity=conductivities)
