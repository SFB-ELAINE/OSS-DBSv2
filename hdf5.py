import h5py
import numpy as np

# with h5py.File("mytestfile.hdf5", "w") as f:
#     dset = f.create_dataset("mydataset", (100,), dtype='int')

path = 'mytestfile.hdf5'
path = 'Johnson_Butenko_axon_McIntyre.h5'

f = h5py.File(path, 'r')
print(list(f.keys()))

data = f['CF_21']
print(data)

# print(data.shape)

# print(data.dtype)

# print(type(data))

# print(data[0])

# print(np.array(data))

# print(f.name)

# print(data.name)

f.close()