# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

"""Tests for the volume conductor formulation strategies."""

import pytest

formulations = pytest.importorskip(
    "ossdbs.fem.volume_conductor.formulations",
    reason="formulations module not yet implemented",
)
FloatingFormulation = formulations.FloatingFormulation
FloatingImpedanceFormulation = formulations.FloatingImpedanceFormulation
FormulationStrategy = formulations.FormulationStrategy
NonFloatingFormulation = formulations.NonFloatingFormulation
get_formulation = formulations.get_formulation


class TestFormulationStrategy:
    """Tests for FormulationStrategy ABC."""

    def test_abstract_class(self):
        """Test that FormulationStrategy cannot be instantiated."""
        with pytest.raises(TypeError):
            FormulationStrategy()


class TestNonFloatingFormulation:
    """Tests for NonFloatingFormulation."""

    def test_name(self):
        """Test formulation name."""
        formulation = NonFloatingFormulation()
        assert formulation.name == "NonFloating"

    def test_is_formulation_strategy(self):
        """Test that it's a FormulationStrategy subclass."""
        formulation = NonFloatingFormulation()
        assert isinstance(formulation, FormulationStrategy)


class TestFloatingFormulation:
    """Tests for FloatingFormulation."""

    def test_name(self):
        """Test formulation name."""
        formulation = FloatingFormulation()
        assert formulation.name == "Floating"

    def test_is_formulation_strategy(self):
        """Test that it's a FormulationStrategy subclass."""
        formulation = FloatingFormulation()
        assert isinstance(formulation, FormulationStrategy)


class TestFloatingImpedanceFormulation:
    """Tests for FloatingImpedanceFormulation."""

    def test_name(self):
        """Test formulation name."""
        formulation = FloatingImpedanceFormulation()
        assert formulation.name == "FloatingImpedance"

    def test_is_formulation_strategy(self):
        """Test that it's a FormulationStrategy subclass."""
        formulation = FloatingImpedanceFormulation()
        assert isinstance(formulation, FormulationStrategy)

    def test_initial_state(self):
        """Test initial state of formulation."""
        formulation = FloatingImpedanceFormulation()
        assert formulation._surface_impedance_boundaries == []
        assert formulation._surface_impedances == {}

    def test_surface_impedance_boundaries_returns_copy(self):
        """Test that surface_impedance_boundaries returns a copy."""
        formulation = FloatingImpedanceFormulation()
        formulation._surface_impedance_boundaries = ["contact1", "contact2"]
        boundaries = formulation.surface_impedance_boundaries
        boundaries.clear()
        assert len(formulation._surface_impedance_boundaries) == 2

    def test_create_bilinear_form_without_impedances_raises(self):
        """Test that create_bilinear_form raises if impedances not set."""
        formulation = FloatingImpedanceFormulation()
        with pytest.raises(RuntimeError, match="Surface impedances not set"):
            formulation.create_bilinear_form(None, None, None)


class TestGetFormulation:
    """Tests for get_formulation factory function."""

    def test_get_nonfloating(self):
        """Test getting non-floating formulation."""
        formulation = get_formulation("NonFloating")
        assert isinstance(formulation, NonFloatingFormulation)

    def test_get_floating(self):
        """Test getting floating formulation."""
        formulation = get_formulation("Floating")
        assert isinstance(formulation, FloatingFormulation)

    def test_get_default_is_nonfloating(self):
        """Test that unknown mode returns non-floating."""
        formulation = get_formulation("Unknown")
        assert isinstance(formulation, NonFloatingFormulation)

    def test_get_floating_impedance(self):
        """Test getting floating impedance formulation."""
        formulation = get_formulation("FloatingImpedance")
        assert isinstance(formulation, FloatingImpedanceFormulation)
