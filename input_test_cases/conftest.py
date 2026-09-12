# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

"""Pytest configuration for OSS-DBS simulation tests."""

import importlib
from pathlib import Path

import pytest

NEURON_AVAILABLE = importlib.util.find_spec("neuron") is not None


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", "simulation: full simulation tests (run on PR only)"
    )
    config.addinivalue_line("markers", "slow: slow-running tests (> 1 min)")
    config.addinivalue_line(
        "markers", "requires_neuron: tests that require NEURON installation"
    )
    config.addinivalue_line("markers", "vta: VTA-related tests")
    config.addinivalue_line("markers", "pam: Pathway Activation Modeling tests")
    config.addinivalue_line("markers", "floating: floating contact tests")
    config.addinivalue_line("markers", "surface_impedance: surface impedance tests")


def pytest_collection_modifyitems(config, items):
    """Skip NEURON tests if NEURON is not available."""
    skip_neuron = pytest.mark.skip(reason="NEURON not installed")
    for item in items:
        if "requires_neuron" in item.keywords and not NEURON_AVAILABLE:
            item.add_marker(skip_neuron)


@pytest.fixture(scope="session")
def base_path():
    """Path to input_test_cases directory."""
    return Path(__file__).parent
