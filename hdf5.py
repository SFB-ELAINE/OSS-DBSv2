import h5py
import numpy as np

# with h5py.File("test.hdf5", "w") as f:
#     f.create_dataset('CF_21', data=np.array([[14.5, -14.7, -9.1]]))
#     f.create_dataset('CbTh_31', data=np.array([[14.5, -14.7, -9.1]]))
#     f.create_dataset('EPN2VA_VL_21', data=np.array([[14.5, -14.7, -9.1]]))
#     f.create_dataset('HDP_21', data=np.array([[14.5, -14.7, -9.1],
#                                               [-100, -100, -100]]))


with h5py.File("test.hdf5", "w") as f:
    group = f.create_group('Population_1')
    sub_group = group.create_group('Axon')
    sub_group.create_dataset('Time', data=np.array([[1, 2, 3]]))
    group.create_dataset('Axon_1', data=np.array([[1, 2, 3]]))
    group.create_dataset('Axon_2', data=np.array([[4, 5, 6]]))
    group = f.create_group('Population_2')
    group.create_dataset('Axon_1', data=np.array([[1, 2, 3], [4, 5, 6]]))
    group.create_dataset('Axon_2', data=np.array([[4, 5, 6]]))


f = h5py.File("test.hdf5", 'r')
print(list(f.keys()))

print(list(f['Population_1']['Axon'].keys())[0], np.array(f['Population_1']['Axon']['Time']))

for group in f.keys():
    for sub_group in f[group].keys():
        print(sub_group, np.array(f[group][sub_group]))

# for group in f.keys():
#     if group == 'TimeSteps':
#         print(group, np.array(f[group]))
#         continue
    
#     print(group)
#     for sub_group in f[group].keys():
#         print(sub_group, np.array(f[group][sub_group]).shape)



# path = 'Johnson_Butenko_axon_McIntyre.h5'
# # path = 'Johnson_Butenko_axon_McIntyre.h5'
# # path = 'solution_130.0.h5'
# f = h5py.File(path, 'r')
# print(list(f.keys()))

# points = np.concatenate([f[key] for key in f.keys()])
# n_points = points.shape[0]

# print(points.shape)

f.close()