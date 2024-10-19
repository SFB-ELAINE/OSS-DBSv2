"""Example script to obtain VTA volume."""

from ossdbs import VTAImage

vta_voxel = VTAImage("Results_VTA/VTA_solution_Lattice.nii")
vta_voxel_ooc = VTAImage("Results_VTA_OOC/VTA_solution_Lattice.nii")
print("VTA volume: ", vta_voxel.get_vta_volume())
print("VTA volume OOC: ", vta_voxel.get_vta_volume())
