import netgen
import ngsolve
from ossdbs.model_geometry import BrainGeometry
from ossdbs.fem import Mesh


class MeshFactory:

    def __init__(self, geometry: BrainGeometry) -> None:
        self.__geometry = geometry

    def create(self, mesh_parameters: dict) -> Mesh:

        if mesh_parameters['LoadMesh']:
            return self.__load_mesh(mesh_parameters)

        return self.__create_mesh(mesh_parameters)

    def __create_mesh(self, mesh_parameters: dict) -> Mesh:
        netgen_geometry = self.__geometry.geometry()
        parameters = self.__meshing_parameters(mesh_parameters)
        ng_mesh = netgen_geometry.GenerateMesh(parameters)
        ngsolve_mesh = ngsolve.Mesh(ngmesh=ng_mesh)
        return Mesh(ngsolve_mesh, mesh_parameters["MeshElementOrder"])

    def __load_mesh(self, mesh_parameters: dict) -> Mesh:
        netgen_geometry = self.__geometry.geometry()
        file_path = mesh_parameters['LoadPath']
        ngsolve_mesh = ngsolve.Mesh(filename=file_path)
        ngsolve_mesh.ngmesh.SetGeometry(netgen_geometry)
        return Mesh(ngsolve_mesh, mesh_parameters["MeshElementOrder"])

    def __meshing_parameters(self, mesh_parameters: dict):
        mesh_type = mesh_parameters['MeshingHypothesis']['Type']
        max_h = mesh_parameters['MeshingHypothesis']['MaxMeshSizeHeight']
        return {'Coarse': netgen.meshing.meshsize.coarse,
                'Fine': netgen.meshing.meshsize.fine,
                'Moderate': netgen.meshing.meshsize.moderate,
                'VeryCoarse': netgen.meshing.meshsize.very_coarse,
                'VeryFine': netgen.meshing.meshsize.very_fine,
                'Default': netgen.meshing.MeshingParameters(),
                'Custom': netgen.meshing.MeshingParameters(max_h=max_h)
                }[mesh_type]
