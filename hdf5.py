import h5py
import numpy as np

# with h5py.File("result.hdf5", "w") as f:
#     f.create_dataset("points", data=np.array([(1, 2, 3)]))
#     f.create_dataset("frequencies", data=np.array([1, 2, 3]))
#     f.create_dataset("potential_values", data=np.array([1 + 3j]))
#     f.create_dataset("current_density_values", data=np.array([1 + 3j]))


path = 'result.hdf5'
# path = 'Johnson_Butenko_axon_McIntyre.h5'
# path = 'solution_130.0.h5'
f = h5py.File(path, 'r')
print(list(f.keys()))
print(np.array(f['frequencies']))
print(np.array(f['points']))
print(np.array(f['potential_values']))
print(np.array(f['frequencies']).shape)
print(np.array(f['potential_values']).shape)
# data = f['mydataset']
# print(np.array(data))

# print(data.shape)

# print(data.dtype)

# print(type(data))

# print(data[0])

# print(np.array(data))

# print(f.name)

# print(data.name)

f.close()