import numpy as np
import pytest

import ossdbs
from ossdbs.fem import Mesh
from ossdbs.point_analysis import Lattice, Pathway, VoxelLattice
from ossdbs.stimulation_signals import FrequencyDomainSignal, RectangleSignal
from ossdbs.utils.nifti1image import MagneticResonanceImage
from ossdbs.utils.settings import Settings

# ruff: noqa: F401
from .test_fem import geometry_fixture, mri_fixture, settings_fixture


# ruff: noqa: F811
@pytest.fixture
def mesh_fixture(geometry_fixture, settings_fixture):
    geometry = geometry_fixture[2].geometry
    mesh = Mesh(geometry, settings_fixture["FEMOrder"])
    mesh.generate_mesh({"MeshingHypothesis": {"Type": "Moderate"}})
    return mesh


# ruff: noqa: F811
@pytest.fixture
def conductivity_fixture(mri_fixture, geometry_fixture, settings_fixture):
    mri_image, _ = mri_fixture

    brain_region, _, geometry = geometry_fixture
    dielectric_model = ossdbs.prepare_dielectric_properties(settings_fixture)
    materials = settings_fixture["MaterialDistribution"]["MRIMapping"]
    return ossdbs.ConductivityCF(
        mri_image,
        brain_region,
        dielectric_model,
        materials,
        geometry.encapsulation_layers,
        complex_data=settings_fixture["EQSMode"],
    )


@pytest.fixture
def settings():
    settings = Settings({}).complete_settings()
    settings["PointModel"]["Lattice"]["Active"] = True
    settings["PointModel"]["VoxelLattice"]["Active"] = True
    settings["MaterialDistribution"]["MRIPath"] = "input_files/homogeneous.nii.gz"
    settings["PointModel"]["Pathway"]["Active"] = True
    settings["PointModel"]["Pathway"]["FileName"] = "input_files/sub-John_Doe/data.h5"

    return settings


@pytest.fixture
def parameters(settings):
    shape_par = settings["PointModel"]["Lattice"]["Shape"]
    shape = shape_par["x"], shape_par["y"], shape_par["z"]
    center_par = settings["PointModel"]["Lattice"]["Center"]
    center = center_par["x[mm]"], center_par["y[mm]"], center_par["z[mm]"]
    dir_par = settings["PointModel"]["Lattice"]["Direction"]
    direction = dir_par["x[mm]"], dir_par["y[mm]"], dir_par["z[mm]"]
    distance = settings["PointModel"]["Lattice"]["PointDistance[mm]"]
    collapse_vta = settings["PointModel"]["Lattice"]["CollapseVTA"]
    return shape, center, distance, direction, collapse_vta


@pytest.fixture
def pathway_fixture(settings):
    return Pathway(input_path=settings["PointModel"]["Pathway"]["FileName"])


class TestPointAnalysis:
    def test_pathway(self, pathway_fixture):
        try:
            pathway = pathway_fixture
            assert pathway is not None
        except Exception:
            pytest.fail("Cannot be instantiated.")

    def test_pathway_signal_assignment(
        self, pathway_fixture, mesh_fixture, conductivity_fixture
    ):
        pathway = pathway_fixture
        signal = RectangleSignal(
            frequency=130,
            # larger than usual for lower cutoff and less computational complexity
            pulse_width=600e-6,
            counter_pulse_width=0.0,
            inter_pulse_width=0.0,
            counter_pulse_amplitude=0.0,
        )
        base_frequency = signal.frequency
        cutoff_frequency = 3e4
        fft_frequencies, fft_coefficients, signal_length = signal.get_fft_spectrum(
            cutoff_frequency
        )

        frequency_domain_signal = FrequencyDomainSignal(
            frequencies=fft_frequencies,
            amplitudes=fft_coefficients,
            base_frequency=base_frequency,
            cutoff_frequency=cutoff_frequency,
            signal_length=signal_length,
            current_controlled=False,
        )
        mesh = mesh_fixture
        conductivity = conductivity_fixture
        # prepare VCM specific data structure
        pathway.prepare_VCM_specific_evaluation(mesh, conductivity)
        pathway.prepare_frequency_domain_data_structure(
            len(frequency_domain_signal.frequencies)
        )

        amplitudes = np.linspace(
            start=0, stop=len(pathway.lattice), num=len(pathway.lattice)
        )

        # copy signal (emulates simulations by using different amplitudes)
        for (freq_idx, _), scale_factor in zip(
            enumerate(fft_frequencies), fft_coefficients, strict=True
        ):
            pots = np.expand_dims(scale_factor * amplitudes, axis=-1)
            pathway.copy_frequency_domain_solution_from_vcm(freq_idx, pots)
        potentials, _, _, _ = pathway.compute_solutions_in_time_domain(
            signal_length, convert_field=False
        )

        # compute original time domain signal td_signal
        cutoff_frequency = signal.get_adjusted_cutoff_frequency(cutoff_frequency)
        dt = 1.0 / cutoff_frequency
        td_signal = signal.get_time_domain_signal(dt=dt, timesteps=signal_length)

        # go through all lattice points and check that fft yields correct signal at all points
        test_values = []
        for amplitude, potential in zip(amplitudes, potentials, strict=True):
            test_values.append(np.all(np.isclose(potential, amplitude * td_signal)))
        assert all(test_values)

    def test_lattice(self, parameters):
        try:
            shape, center, distance, direction, collapse_vta = parameters

            lattice = Lattice(
                shape=shape,
                center=center,
                distance=distance,
                direction=direction,
                collapse_vta=collapse_vta,
            )
            assert lattice is not None
        except Exception:
            pytest.fail("Cannot be instantiated.")

    def test_voxelLattice(self, settings, parameters):
        try:
            _, center, _, _, _ = parameters
            mri_image = MagneticResonanceImage(
                settings["MaterialDistribution"]["MRIPath"]
            )
            affine = mri_image.affine
            header = mri_image.header
            voxel_shape_par = settings["PointModel"]["VoxelLattice"]["Shape"]
            voxel_shape = np.array(
                [
                    voxel_shape_par["x"] + 1,
                    voxel_shape_par["y"] + 1,
                    voxel_shape_par["z"] + 1,
                ]
            )

            voxelLattice = VoxelLattice(center, affine, voxel_shape, header)
            assert voxelLattice is not None
        except Exception:
            pytest.fail("Cannot be instantiated.")


class TestScaleFactorOnCopy:
    """``scale_factor`` is folded into the write rather than pre-multiplied.

    A solved octave band covers many spectrum indices that share the same
    potentials and fields and differ only by that scalar, so scaling the
    arrays before the call allocated a full-size temporary per index. The
    values written must be unchanged by that: the factor is complex in
    current-controlled mode, where dropping the imaginary part would rotate
    the phase of every harmonic in the band.
    """

    @staticmethod
    def _model(n_points=4, signal_length=3):
        from ossdbs.point_analysis.point_model import PointModel

        class _Model:
            copy_frequency_domain_solution_from_vcm = (
                PointModel.copy_frequency_domain_solution_from_vcm
            )

            def __init__(self):
                shape = (n_points, signal_length)
                self.tmp_potential_freq_domain = np.zeros(shape, dtype=complex)
                self.tmp_Ex_freq_domain = np.zeros(shape, dtype=complex)
                self.tmp_Ey_freq_domain = np.zeros(shape, dtype=complex)
                self.tmp_Ez_freq_domain = np.zeros(shape, dtype=complex)

        return _Model()

    @staticmethod
    def _inputs(n_points=4):
        rng = np.random.default_rng(20260910)
        potentials = rng.standard_normal((n_points, 1)) + 1j * rng.standard_normal(
            (n_points, 1)
        )
        fields = rng.standard_normal((n_points, 3)) + 1j * rng.standard_normal(
            (n_points, 3)
        )
        return potentials, fields

    @pytest.mark.parametrize("scale_factor", [1.0, 2.5, -0.75, 1e-3, 0.5 + 2.0j, -1.5j])
    def test_matches_pre_multiplication(self, scale_factor):
        """Folding the factor into the write equals scaling the array first."""
        potentials, fields = self._inputs()

        folded = self._model()
        folded.copy_frequency_domain_solution_from_vcm(
            1, potentials, fields, scale_factor=scale_factor
        )

        pre_multiplied = self._model()
        pre_multiplied.copy_frequency_domain_solution_from_vcm(
            1, scale_factor * potentials, scale_factor * fields
        )

        for component in ("potential", "Ex", "Ey", "Ez"):
            name = f"tmp_{component}_freq_domain"
            np.testing.assert_allclose(
                getattr(folded, name), getattr(pre_multiplied, name)
            )

    def test_complex_factor_is_not_truncated(self):
        """A real-only write would silently drop the phase."""
        potentials, fields = self._inputs()
        scale_factor = 0.5 + 2.0j

        model = self._model()
        model.copy_frequency_domain_solution_from_vcm(
            0, potentials, fields, scale_factor=scale_factor
        )

        np.testing.assert_allclose(
            model.tmp_potential_freq_domain[:, 0], scale_factor * potentials[:, 0]
        )
        np.testing.assert_allclose(
            model.tmp_Ez_freq_domain[:, 0], scale_factor * fields[:, 2]
        )
        assert np.any(model.tmp_potential_freq_domain[:, 0].imag != 0)

    def test_default_factor_leaves_values_untouched(self):
        """Voltage-controlled runs pass no factor at all."""
        potentials, fields = self._inputs()

        model = self._model()
        model.copy_frequency_domain_solution_from_vcm(2, potentials, fields)

        np.testing.assert_allclose(
            model.tmp_potential_freq_domain[:, 2], potentials[:, 0]
        )
        np.testing.assert_allclose(model.tmp_Ex_freq_domain[:, 2], fields[:, 0])

    def test_other_columns_are_not_written(self):
        """Each index writes exactly its own column of the band."""
        potentials, fields = self._inputs()

        model = self._model()
        model.copy_frequency_domain_solution_from_vcm(
            1, potentials, fields, scale_factor=3.0
        )

        assert np.all(model.tmp_potential_freq_domain[:, 0] == 0)
        assert np.all(model.tmp_potential_freq_domain[:, 2] == 0)

    def test_fields_none_writes_potential_only(self):
        """VTA runs without field export must not touch the field arrays."""
        potentials, _ = self._inputs()

        model = self._model()
        model.copy_frequency_domain_solution_from_vcm(
            0, potentials, None, scale_factor=2.0
        )

        np.testing.assert_allclose(
            model.tmp_potential_freq_domain[:, 0], 2.0 * potentials[:, 0]
        )
        assert np.all(model.tmp_Ex_freq_domain == 0)
        assert np.all(model.tmp_Ey_freq_domain == 0)
        assert np.all(model.tmp_Ez_freq_domain == 0)
