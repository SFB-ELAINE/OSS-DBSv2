import numpy as np
import nibabel


def main():
    x = 50
    y = 50
    z = 50
    box = np.full((x, y, z), 1)
    array = np.empty((2*x, 2*y, 2*z))
    array[:x, :y, :z] = box * 0
    array[x:, :y, :z] = box * 1
    array[:x, y:, :z] = box * 2
    array[x:, y:, :z] = box * 3
    array[:x, :y, z:] = box * 3
    array[x:, :y, z:] = box * 2
    array[:x, y:, z:] = box * 1
    array[x:, y:, z:] = box * 0

    affine = np.array([[0.1, 0, 0, 0],
                       [0, 0.1, 0, 0],
                       [0, 0, 0.1, 0],
                       [0, 0, 0, 1]])
    nii_image = nibabel.Nifti1Image(dataobj=array, affine=affine)
    nii_image.header['xyzt_units'] = 2
    nibabel.save(nii_image, 'TestMRI')

    array = np.full((x, y, z, 3, 3), np.eye(3))
    nii_image = nibabel.Nifti1Image(dataobj=array, affine=affine)
    nii_image.header['xyzt_units'] = 2
    nibabel.save(nii_image, 'TestDTI')


if __name__ == '__main__':
    main()
