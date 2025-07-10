import numpy as np
import pytest

from ossdbs.point_analysis import Lattice, Pathway, VoxelLattice
from ossdbs.utils.nifti1image import MagneticResonanceImage
from ossdbs.utils.settings import Settings


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


class TestPointAnalysis:
    def test_pathway(self, settings):
        try:
            pathway = Pathway(input_path=settings["PointModel"]["Pathway"]["FileName"])
            assert pathway is not None
        except Exception:
            pytest.fail("Cannot be instantiated.")

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
            shape, center, distance, direction, collapse_vta = parameters
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
