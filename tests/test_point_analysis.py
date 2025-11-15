import numpy as np
import pytest

import ossdbs
from ossdbs.fem import Mesh
from ossdbs.point_analysis import Lattice, Pathway, VoxelLattice
from ossdbs.stimulation_signals import FrequencyDomainSignal, RectangleSignal
from ossdbs.utils.nifti1image import MagneticResonanceImage
from ossdbs.utils.settings import Settings


@pytest.fixture
def mesh_fixture(geometry_fixture, settings_fixture):
    geometry = geometry_fixture[2].geometry
    mesh = Mesh(geometry, settings_fixture["FEMOrder"])
    mesh.generate_mesh({"MeshingHypothesis": {"Type": "Moderate"}})
    return mesh


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
    settings["PointModel"]["Pathway"]["FileName"] = "input_files/data.h5"

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
            pulse_width=60e-6,
            counter_pulse_width=0.0,
            inter_pulse_width=0.0,
            counter_pulse_amplitude=0.0,
        )
        base_frequency = signal.frequency
        cutoff_frequency = 5e5
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

        # TODO continue

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
