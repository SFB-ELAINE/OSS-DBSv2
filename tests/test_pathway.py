import h5py
import numpy as np
import pytest

from ossdbs.point_analysis.pathway import Pathway


def create_dummy_h5(filepath):
    with h5py.File(filepath, "w") as f:
        grp1 = f.create_group("pop1")
        grp1.create_dataset("axon0", data=np.array([[0, 0, 0], [1, 1, 1]]))
        grp1.create_dataset("axon1", data=np.array([[2, 2, 2], [300, 3, 3]]))


def test_pathway_init_and_population_names(tmp_path):
    h5file = tmp_path / "dummy.h5"
    create_dummy_h5(h5file)
    pathway = Pathway(h5file)
    assert pathway.get_population_names() == ["pop1"]
    axon_names = pathway.get_axon_names()
    assert axon_names == [["axon0", "axon1"]]
    axon_numbers = pathway.get_axon_numbers()
    assert axon_numbers == [2]


def test_create_index_and_axon_length(tmp_path):
    h5file = tmp_path / "dummy.h5"
    create_dummy_h5(h5file)
    pathway = Pathway(h5file)
    lattice = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3]])
    index = pathway.create_index(lattice)
    assert index.shape == (4, 1)
    assert np.all(index[:2] == 0)
    assert np.all(index[2:] == 1)
    assert pathway.get_axon_length() == 2


def test_filter_for_geometry(tmp_path):
    h5file = tmp_path / "dummy.h5"
    create_dummy_h5(h5file)
    pathway = Pathway(h5file)
    # Simulate a masked array: first axon inside, second axon outside
    data = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3]])
    mask = np.array([[False], [False], [True], [True]])
    grid_pts = np.ma.MaskedArray(data, mask=np.broadcast_to(mask, data.shape))
    filtered = pathway.filter_for_geometry(grid_pts)
    # Only first axon should remain
    assert filtered.shape == (2, 3)
    assert np.allclose(filtered, np.array([[0, 0, 0], [1, 1, 1]]))


def test_filter_csf_encap(tmp_path):
    h5file = tmp_path / "dummy.h5"
    create_dummy_h5(h5file)
    pathway = Pathway(h5file)
    # 4 points, first two are axon0, last two axon1
    inside_csf = np.array([0, 1, 0, 0])
    inside_encap = np.array([0, 0, 1, 0])
    pathway.filter_csf_encap(inside_csf, inside_encap)
    pop = pathway._populations[0]
    assert pop.axons[0].status == -2  # axon0: one point in csf
    assert pop.axons[1].status == -1  # axon1: one point in encap


def test_save_as_nifti(tmp_path):
    h5file = tmp_path / "dummy.h5"
    create_dummy_h5(h5file)
    pathway = Pathway(h5file)
    with pytest.raises(NotImplementedError):
        pathway.save_as_nifti(np.array([1, 2, 3]), "dummy.nii")
