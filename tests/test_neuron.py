import pytest

from ossdbs.axon_processing.neuron_model import McNeal1976, NeuronSimulator

neuron = pytest.importorskip("neuron")
pathway_dict = {}
pathway_dict_McNeal = {
    "n_Ranvier": [1],
    "axon_diams": [5.0],
    "Axon_Model_Type": "McNeal1976",
    "Name_prepared_neuron_array": "Allocated_axons.h5",
    "Neuron_model_array_prepared": True,
    "N_seeded_neurons": [1],
    "N_orig_neurons": [1],
    "connectome_name": "TestConnectome",
}


def test_neuron_model(tmp_path):
    with pytest.raises(TypeError):
        _ = NeuronSimulator(pathway_dict, tmp_path)


def test_mc_neal_model(tmp_path):
    _ = McNeal1976(pathway_dict_McNeal, tmp_path)
