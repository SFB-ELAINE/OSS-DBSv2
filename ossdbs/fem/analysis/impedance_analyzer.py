# Copyright 2026 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import logging
import os
from collections.abc import Sequence

import ngsolve
import numpy as np
import pandas as pd

from ossdbs.fem.volume_conductor.volume_conductor_model import VolumeConductor

_logger = logging.getLogger(__name__)


class ImpedanceAnalyzer:
    """Standalone multicontact admittance / impedance analysis tool.

    Computes the full NxN admittance matrix Y via the superposition
    approach (see ``examples/MulticontactCurrents/Superposition-approach.ipynb``)
    across a list of frequencies, and the impedance matrix Z = Y^-1.

    This analyzer is decoupled from ``VolumeConductor.run_full_analysis``
    and has no effect on the stimulation pipeline (no coupling to
    ``_scale_factor``, ``_stimulation_variable``, or point models).
    It shares the volume conductor's mesh, conductivity, contacts, and
    solver.

    Parameters
    ----------
    volume_conductor : VolumeConductor
        Fully constructed volume conductor (mesh generated, HP applied).
    """

    def __init__(
        self,
        volume_conductor: VolumeConductor,
        include_floating: bool = True,
    ) -> None:
        self._vc = volume_conductor
        self._include_floating = include_floating
        self._frequencies: np.ndarray | None = None
        self._admittance_matrices: np.ndarray | None = None
        self._impedance_matrices: np.ndarray | None = None

    @property
    def contact_names(self) -> list[str]:
        """Ordered list of contact names participating in the analysis.

        Includes active contacts, and floating contacts if
        ``include_floating`` was True (default).
        """
        vc = self._vc
        contacts = vc.contacts.active
        if self._include_floating:
            contacts = contacts + vc.contacts.floating
        return [c.name for c in contacts]

    @property
    def frequencies(self) -> np.ndarray:
        """Frequencies [Hz] at which the analysis was performed."""
        if self._frequencies is None:
            raise RuntimeError("compute() must be called before accessing results.")
        return self._frequencies

    @property
    def admittance_matrices(self) -> np.ndarray:
        """Stacked admittance matrices with shape ``(n_freq, N, N)``."""
        if self._admittance_matrices is None:
            raise RuntimeError("compute() must be called before accessing results.")
        return self._admittance_matrices

    @property
    def impedance_matrices(self) -> np.ndarray:
        """Stacked impedance matrices ``Z = Y^-1`` with shape ``(n_freq, N, N)``."""
        if self._impedance_matrices is None:
            raise RuntimeError("compute() must be called before accessing results.")
        return self._impedance_matrices

    def compute(self, frequencies: Sequence[float]) -> None:
        """Compute Y and Z = Y^-1 for each requested frequency.

        Parameters
        ----------
        frequencies : sequence of float
            Frequencies [Hz] at which to evaluate the admittance matrix.
        """
        freqs = np.asarray(frequencies, dtype=float)
        names = self.contact_names
        N = len(names)
        if N < 2:
            raise ValueError(
                f"Admittance analysis requires at least 2 contacts, got {N}."
            )
        dtype = complex if self._vc.is_complex else float
        Y = np.zeros((len(freqs), N, N), dtype=dtype)
        Z = np.zeros_like(Y)
        for idx, freq in enumerate(freqs):
            Y[idx] = self._compute_Y_at_frequency(freq, names)
            Z[idx] = np.linalg.inv(Y[idx])
        self._frequencies = freqs
        self._admittance_matrices = Y
        self._impedance_matrices = Z

    def export(self, output_path: str) -> None:
        """Write ``admittance_matrix.csv`` and ``impedance_matrix.csv``.

        Both files share the flat schema ``freq,row,col,real,imag`` and
        are written into ``output_path``.
        """
        if self._admittance_matrices is None:
            raise RuntimeError("compute() must be called before export().")
        names = self.contact_names
        freqs = self._frequencies
        os.makedirs(output_path, exist_ok=True)
        for label, matrices in [
            ("admittance_matrix", self._admittance_matrices),
            ("impedance_matrix", self._impedance_matrices),
        ]:
            rows = []
            for freq_idx, freq in enumerate(freqs):
                mat = matrices[freq_idx]
                for i, row_name in enumerate(names):
                    for j, col_name in enumerate(names):
                        rows.append(
                            {
                                "freq": freq,
                                "row": row_name,
                                "col": col_name,
                                "real": mat[i, j].real,
                                "imag": mat[i, j].imag
                                if np.iscomplex(mat[i, j])
                                else 0.0,
                            }
                        )
            df = pd.DataFrame(rows)
            df.to_csv(os.path.join(output_path, f"{label}.csv"), index=False)
            _logger.info(f"Saved {label} to {label}.csv")

    def _compute_Y_at_frequency(
        self, frequency: float, contact_names: list[str]
    ) -> np.ndarray:
        vc = self._vc
        N = len(contact_names)
        dtype = complex if vc.is_complex else float
        Y = np.zeros((N, N), dtype=dtype)

        _logger.info(
            f"Computing {N}x{N} admittance matrix "
            f"({N * (N + 1) // 2} BVPs) at {frequency} Hz"
        )

        sigma = vc.conductivity_cf(vc.mesh, frequency)
        surface_impedances = vc.contacts.get_surface_impedances(
            frequency, is_complex=vc.is_complex
        )
        if not any(v is not None for v in surface_impedances.values()):
            surface_impedances = None

        # Include the outer brain surface as Dirichlet u=0 (ground reference)
        # so that Y is non-singular and all contacts appear in Z.
        brain_surfaces = vc.model_geometry.brain_surface_names
        dirichlet_boundaries = [*contact_names, *brain_surfaces]
        space = vc.h1_space(boundaries=dirichlet_boundaries, is_complex=vc.is_complex)
        gfu = ngsolve.GridFunction(space=space)

        # Diagonal entries: Y_ii = P / Vi^2
        # with P = 0.5 * integral of sigma|grad u|^2 dV.
        Vi = 1.0
        for i in range(N):
            boundary_values = dict.fromkeys(contact_names, 0.0)
            boundary_values[contact_names[i]] = Vi
            power = 0.5 * self._solve_bvp(
                space, gfu, boundary_values, sigma, surface_impedances
            )
            Y[i, i] = 2.0 * power / (Vi * Vi)
            _logger.debug(f"Y[{i},{i}] = {Y[i, i]}")

        # Off-diagonal entries via superposition.
        Vi = 1.0
        Vj = 2.0
        for i in range(N):
            for j in range(i + 1, N):
                boundary_values = dict.fromkeys(contact_names, 0.0)
                boundary_values[contact_names[i]] = Vi
                boundary_values[contact_names[j]] = Vj
                power = 0.5 * self._solve_bvp(
                    space, gfu, boundary_values, sigma, surface_impedances
                )
                Y[i, j] = power / (Vi * Vj) - 0.5 * (
                    Vi / Vj * Y[i, i] + Vj / Vi * Y[j, j]
                )
                Y[j, i] = Y[i, j]
                _logger.debug(f"Y[{i},{j}] = Y[{j},{i}] = {Y[i, j]}")

        _logger.info(f"Admittance matrix:\n{Y}")
        return Y

    def _solve_bvp(
        self,
        space: ngsolve.H1,
        gfu: ngsolve.GridFunction,
        boundary_values: dict,
        sigma: ngsolve.CoefficientFunction,
        surface_impedances: dict | None,
    ) -> complex:
        """Solve one Dirichlet BVP and return total dissipated power."""
        vc = self._vc
        coefficient = vc.mesh.boundary_coefficients(boundary_values)
        gfu.Set(coefficient, VOL_or_BND=ngsolve.BND)

        u = space.TrialFunction()
        v = space.TestFunction()

        bilinear_form = ngsolve.BilinearForm(space=space)
        bilinear_form += sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx

        linear_form = ngsolve.LinearForm(space=space)

        if surface_impedances is not None:
            for contact_name, voltage in boundary_values.items():
                zs = surface_impedances.get(contact_name)
                if zs is None:
                    continue
                ys = ngsolve.CF(1.0 / zs)
                bilinear_form += ys * u * v * ngsolve.ds(contact_name)
                linear_form += ys * ngsolve.CF(voltage) * v * ngsolve.ds(contact_name)

        vc.solver.bvp(bilinear_form, linear_form, gfu)

        E = -ngsolve.grad(gfu)
        J = sigma * E
        power = ngsolve.Integrate(ngsolve.Conj(E) * J, vc.mesh.ngsolvemesh)

        if surface_impedances is not None:
            for contact_name, voltage in boundary_values.items():
                zs = surface_impedances.get(contact_name)
                if zs is None:
                    continue
                ys = ngsolve.CF(1.0 / zs)
                diff = vc.mesh.boundary_coefficients({contact_name: voltage}) - gfu
                power += ngsolve.Integrate(
                    ys * diff * ngsolve.Conj(diff),
                    mesh=vc.mesh.ngsolvemesh,
                    definedon=vc.mesh.ngsolvemesh.Boundaries(contact_name),
                )

        return power
