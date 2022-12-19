import os
import ngsolve


class Output:

    def __init__(self, mesh, potential) -> None:
        self.__mesh = mesh
        self.__potential = potential

    def save(self, path: str = ''):

        file_name = os.path.basename(path)

        if not file_name:
            file_name = 'result'

        file_dir = os.path.dirname(path)

        if not file_dir:
            file_dir = 'result'

        if not os.path.exists(file_dir):
            os.makedirs(file_dir)

        ngsolve.VTKOutput(ma=self.__mesh,
                          coefs=[self.__potential.real, self.__potential.imag],
                          names=["potential_real", "potential_imag"],
                          filename=os.path.join(file_dir, file_name),
                          subdivision=0
                          ).Do()


class OutputTest:

    def __init__(self, mesh, potential, density) -> None:
        self.__mesh = mesh
        self.__potential = potential
        self.__density = density

    def save(self, path: str = ''):

        file_name = os.path.basename(path)

        if not file_name:
            file_name = 'result'

        file_dir = os.path.dirname(path)

        if not file_dir:
            file_dir = 'result'

        if not os.path.exists(file_dir):
            os.makedirs(file_dir)

        field = ngsolve.grad(self.__potential)
        Power = ngsolve.Integrate(field * ngsolve.Conj(self.__density),
                                  self.__mesh)

        print(1 / Power)
        ngsolve.VTKOutput(ma=self.__mesh,
                          coefs=[self.__potential.real,
                                 self.__potential.imag,
                                 field.real,
                                 field.imag,
                                 self.__density.real,
                                 self.__density.imag],
                          names=["potential_real",
                                 "potential_imag",
                                 "field_real",
                                 "field_imag",
                                 "current_density_real",
                                 "current_density_imag"],
                          filename=os.path.join(file_dir, file_name),
                          subdivision=0
                          ).Do()
