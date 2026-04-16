"""Tests for ossdbs.utils.materials.have_dielectric_properties_changed."""

import pytest

from ossdbs.dielectric_model import ColeCole4Model, default_cole_cole4_parameters
from ossdbs.utils import have_dielectric_properties_changed


class _FakeModel:
    """Minimal model with controllable real and imaginary conductivity."""

    def __init__(self, conductivity_func, complex_conductivity_func):
        self._cond = conductivity_func
        self._ccond = complex_conductivity_func

    def conductivity(self, omega):
        return self._cond(omega)

    def complex_conductivity(self, omega):
        return self._ccond(omega)


class TestHaveDielectricPropertiesChanged:
    """Tests for the frequency-change detection function."""

    def test_same_frequency_returns_false(self):
        """Identical frequencies should never report a change."""
        props = {
            "Gray matter": ColeCole4Model(default_cole_cole4_parameters["Gray matter"]),
        }
        assert not have_dielectric_properties_changed(
            props, is_complex=False, old_freq=130.0, new_freq=130.0, threshold=0.0
        )
        assert not have_dielectric_properties_changed(
            props, is_complex=True, old_freq=130.0, new_freq=130.0, threshold=0.0
        )

    def test_real_mode_detects_real_change(self):
        """Real-only mode detects changes in real conductivity."""
        props = {
            "Gray matter": ColeCole4Model(default_cole_cole4_parameters["Gray matter"]),
        }
        # 1 Hz vs 1 kHz — large change in real conductivity
        assert have_dielectric_properties_changed(
            props, is_complex=False, old_freq=1.0, new_freq=1000.0, threshold=0.01
        )

    def test_complex_mode_detects_imaginary_change(self):
        """Complex mode must detect changes in the imaginary part.

        This is the regression test for the bug where the is_complex branch
        called model.conductivity() (real-only) instead of
        model.complex_conductivity().
        """
        # Model where real part is constant but imaginary part changes
        model = _FakeModel(
            conductivity_func=lambda omega: 1.0,
            complex_conductivity_func=lambda omega: 1.0 + 1j * omega,
        )
        props = {"test": model}

        # Real mode: sees no change (real part is constant)
        assert not have_dielectric_properties_changed(
            props, is_complex=False, old_freq=1.0, new_freq=2.0, threshold=0.01
        )
        # Complex mode: must detect the imaginary-part change
        assert have_dielectric_properties_changed(
            props, is_complex=True, old_freq=1.0, new_freq=2.0, threshold=0.01
        )

    def test_complex_mode_below_threshold(self):
        """Complex mode returns False when change is below threshold."""
        # Nearly constant complex conductivity
        model = _FakeModel(
            conductivity_func=lambda omega: 1.0,
            complex_conductivity_func=lambda omega: 1.0 + 1j * 1.0,
        )
        props = {"test": model}
        assert not have_dielectric_properties_changed(
            props, is_complex=True, old_freq=1.0, new_freq=2.0, threshold=0.01
        )

    def test_multiple_materials_worst_case(self):
        """Detection triggers if *any* material exceeds the threshold."""
        # Material A: constant
        model_a = _FakeModel(
            conductivity_func=lambda omega: 1.0,
            complex_conductivity_func=lambda omega: 1.0 + 0j,
        )
        # Material B: changes a lot
        model_b = _FakeModel(
            conductivity_func=lambda omega: omega,
            complex_conductivity_func=lambda omega: omega + 0j,
        )
        props = {"stable": model_a, "changing": model_b}
        assert have_dielectric_properties_changed(
            props, is_complex=False, old_freq=1.0, new_freq=2.0, threshold=0.01
        )

    def test_with_cole_cole4_complex(self):
        """Integration test: ColeCole4 model with is_complex=True.

        Gray matter has significant imaginary conductivity changes across
        the frequency range used in DBS.
        """
        props = {
            "Gray matter": ColeCole4Model(default_cole_cole4_parameters["Gray matter"]),
        }
        # Distant frequencies — should detect change
        assert have_dielectric_properties_changed(
            props, is_complex=True, old_freq=130.0, new_freq=130000.0, threshold=0.01
        )

    @pytest.mark.parametrize("is_complex", [True, False])
    def test_high_threshold_returns_false(self, is_complex):
        """A very high threshold should not trigger for nearby frequencies."""
        props = {
            "White matter": ColeCole4Model(
                default_cole_cole4_parameters["White matter"]
            ),
        }
        # Nearby frequencies with a generous threshold
        assert not have_dielectric_properties_changed(
            props,
            is_complex=is_complex,
            old_freq=130.0,
            new_freq=131.0,
            threshold=0.1,
        )
