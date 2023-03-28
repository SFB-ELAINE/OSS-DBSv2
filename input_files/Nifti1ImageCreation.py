import numpy as np
import nibabel


def main():
    x = 193
    y = 229
    z = 193
    array = np.full((x, y, z), 1)

    affine = np.array([[1., 0, 0, -96],
                       [0, 1., 0, -132],
                       [0, 0, 1., -74],
                       [0, 0, 0, 1]])
    nii_image = nibabel.Nifti1Image(dataobj=array, affine=affine)
 #  nii_image.header['xyzt_units'] = 2
    nibabel.save(nii_image, 'TestMRI')

    array = np.full((x, y, z, 3, 3), np.eye(3))
    nii_image = nibabel.Nifti1Image(dataobj=array, affine=affine)
    nii_image.header['xyzt_units'] = 2
    nibabel.save(nii_image, 'TestDTI')


if __name__ == '__main__':
    main()
