import h5py
import numpy as np

f_1 = h5py.File("VTA_array.h5", "r")
a = f_1["VTA_array"]
points = np.array(a)
f_1.close()

with h5py.File("VTA_array.hdf5", "w") as f:
	group = f.create_group('population_1')
	group.create_dataset('axon_1', data=points)
	
f.close()
