
from ossdbs.impedance_mode.spectrum import SpectrumMode, ResultImpedance
import numpy as np
import ngsolve

from ossdbs.signals import Signal
from ossdbs.volume_conductor.volume_conductor_model import VolumeConductor


class NoTruncationImpedance(SpectrumMode):

    def compute(self, signal: Signal, volume_conductor: VolumeConductor,):
        n_fft_frequencies = int(self.SPACING_FACTOR // 2) + 1
        n_intervalls = int(np.log10(n_fft_frequencies)) + 1
        frequencies = [0]
        for index in range(n_intervalls):
            interval__freq = signal.frequency * 10 ** index * np.arange(1, 10)
            frequencies.extend(interval__freq)

        mesh = volume_conductor.mesh.ngsolvemesh()
        impedance = np.zeros(len(frequencies), dtype=complex)

        for index, frequency in enumerate(frequencies[:2]):
            solution = volume_conductor.compute_solution(frequency)
            field = ngsolve.grad(solution.potential)
            curr_dens_conj = ngsolve.Conj(solution.current_density)
            power = ngsolve.Integrate(field * curr_dens_conj, mesh)
            impedance[index] = 1 / power if power else 0

        return ResultImpedance(frequency=frequencies, imdedance=impedance)
