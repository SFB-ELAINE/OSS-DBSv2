from src.Nifti1Image import Nifti1Image


class DiffusionTensorImage(Nifti1Image):

    def diffusion_at(self, positions: list) -> list:
        return [[(xx, xy, xz), (xy, yy, yz), (xz, yz, zz)]
                        for xx, xy, xz, yy, yz, zz in self.values_at(positions)]
