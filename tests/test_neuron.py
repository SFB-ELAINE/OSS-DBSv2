import pytest

from ossdbs.axon_processing.neuron_model import NeuronSimulator

neuron = pytest.importorskip("neuron")
pathway_dict = {}


def test_neuron_model(tmp_path):
    with pytest.raises(TypeError):
        _ = NeuronSimulator(pathway_dict, tmp_path)
