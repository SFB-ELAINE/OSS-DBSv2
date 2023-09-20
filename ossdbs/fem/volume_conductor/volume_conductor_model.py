import logging
import os
import time
from abc import ABC, abstractmethod
from typing import List, Optional

import h5py
import ngsolve
import numpy as np
import pandas as pd

from ossdbs.fem.mesh import Mesh
from ossdbs.fem.solver import Solver
from ossdbs.model_geometry import Contacts, ModelGeometry
from ossdbs.point_analysis import Lattice, PointModel, VoxelLattice
from ossdbs.stimulation_signals import FrequencyDomainSignal
from ossdbs.utils.vtk_export import FieldSolution

from .conductivity import ConductivityCF

_logger = logging.getLogger(__name__)


class VolumeConductor(ABC):
    """Template class of a volume conductor.

    Parameters
    ----------
    mesh : Mesh
    conductivity : ConductivityCF
    model_geometry : ModelGeometry
    solver : Solver

    # TODO add more abstractmethod ?
    """

    def __init__(
        self,
        geometry: ModelGeometry,
        conductivity: ConductivityCF,
        solver: Solver,
        order: int,
        meshing_parameters: dict,
        frequency_domain_signal: FrequencyDomainSignal,
    ) -> None:
        self._solver = solver
        self._model_geometry = geometry
        # to update contacts later, shall not be changed
        self._base_contacts = Contacts(geometry.contacts)
        # can be changed at every simulation run
        self._contacts = Contacts(geometry.contacts)
        self._order = order

        self._conductivity_cf = conductivity
        self._complex = conductivity.is_complex
        self._signal = frequency_domain_signal

        self._impedances = None
        self._mesh = Mesh(self._model_geometry.geometry, self._order)
        if meshing_parameters["LoadMesh"]:
            # TODO check that LoadPath is there if LoadMesh
            self._mesh.load_mesh(meshing_parameters["LoadPath"])
        else:
            self._mesh.generate_mesh(meshing_parameters["MeshingHypothesis"])

        if meshing_parameters["SaveMesh"]:
            self._mesh.save(meshing_parameters["SavePath"])

        # to save previous solution and do post-processing
        self._frequency = None
        self._sigma = None

        # to write results
        self._output_path = None

    @abstractmethod
    def compute_solution(self, frequency: float) -> None:
        """Compute solution at frequency."""
        pass

    def run_full_analysis(
        self,
        compute_impedance: bool = False,
        export_vtk: bool = False,
        point_model: PointModel = None,
        template_space: bool = False,
        activation_threshold: Optional[float] = None,
    ) -> None:
        """Run entire volume conductor model.

        Notes
        -----
        TODO full documentation
        """
        timings: dict = {}
        grid_pts = None
        lattice_mask = None
        lattice = None
        if point_model is not None:
            grid_pts = point_model.points_in_mesh(self.mesh)
            lattice_mask = np.invert(grid_pts.mask)
            lattice = point_model.filtered_lattice(grid_pts)
            # TODO how to use these masks?
            self.get_points_in_csf(lattice)
            self.get_points_in_encapsulation_layer(lattice)
        if lattice is not None:
            timings["PotentialProbing"] = []
            timings["FieldProbing"] = []
        if export_vtk:
            timings["FieldExport"] = []
        timings["ComputeSolution"] = []

        if compute_impedance:
            if self.is_complex:
                self._impedances = np.ndarray(
                    shape=(len(self.signal.frequencies)), dtype=complex
                )
            else:
                self._impedances = np.ndarray(shape=(len(self.signal.frequencies)))
        for idx, frequency in enumerate(self.signal.frequencies):
            _logger.info(f"Computing at frequency: {frequency}")
            if not self.current_controlled:
                _logger.debug("Get scaled voltage values")
                voltage_values = self.get_scaled_active_contact_voltages(
                    self.signal.amplitudes[idx]
                )
                self.update_contacts(voltages=voltage_values)
            time_0 = time.time()
            self.compute_solution(frequency)
            time_1 = time.time()
            timings["ComputeSolution"].append(time_1 - time_0)
            time_0 = time_1

            if self.current_controlled:
                # TODO implement current control
                # Steps: compute impedance and rescale
                raise NotImplementedError("Current-controlled mode not yet implemented")
                continue
            if compute_impedance:
                self._impedances[idx] = self.compute_impedance()
            if export_vtk:
                self.vtk_export()
            if lattice is not None:
                potentials = self.evaluate_potential_at_points(lattice)
                time_1 = time.time()
                timings["PotentialProbing"].append(time_1 - time_0)
                time_0 = time_1

                fields = self.evaluate_field_at_points(lattice)
                field_mags = np.linalg.norm(fields, axis=1).reshape(
                    (fields.shape[0], 1)
                )
                time_1 = time.time()
                timings["FieldProbing"].append(time_1 - time_0)
                time_0 = time_1

                # Save points
                h5f_pts = h5py.File(os.path.join(self.output_path, "oss_pts.h5"), "w")
                h5f_pts.create_dataset("points", data=lattice)
                h5f_pts.close()

                # Save potential evaluation
                h5f_pot = h5py.File(
                    os.path.join(self.output_path, "oss_potentials.h5"), "w"
                )
                h5f_pot.create_dataset("points", data=lattice)
                h5f_pot.create_dataset("potentials", data=potentials)
                h5f_pot.close()
                df_pot = pd.DataFrame(
                    np.concatenate(
                        [lattice, potentials.reshape((potentials.shape[0], 1))], axis=1
                    ),
                    columns=["x-pt", "y-pt", "z-pt", "potential"],
                )
                df_pot.to_csv(
                    os.path.join(self.output_path, "oss_potentials.csv"), index=False
                )

                # Save electric field evaluation
                h5f_field = h5py.File(
                    os.path.join(self.output_path, "oss_field.h5"), "w"
                )
                h5f_field.create_dataset("points", data=lattice)
                h5f_field.create_dataset("field/field_vecs", data=fields)
                h5f_field.create_dataset("field/field_mags", data=field_mags)
                h5f_field.close()
                df_field = pd.DataFrame(
                    np.concatenate([lattice, fields, field_mags], axis=1),
                    columns=[
                        "x-pt",
                        "y-pt",
                        "z-pt",
                        "x-field",
                        "y-field",
                        "z-field",
                        "magnitude",
                    ],
                )
                if template_space:
                    df_field.to_csv(
                        os.path.join(self.output_path, "E_field_Template_space.csv"),
                        index=False,
                    )
                else:
                    df_field.to_csv(
                        os.path.join(self.output_path, "E_field_MRI_space.csv"),
                        index=False,
                    )

                if isinstance(point_model, VoxelLattice) or isinstance(
                    point_model, Lattice
                ):
                    field_mags_full = np.zeros(lattice_mask.shape[0], float)
                    # TODO set all values to -1?
                    # field_mags_full -= 1
                    # overwrite values inside mesh
                    field_mags_full[lattice_mask[:, 0]] = field_mags[:, 0]

                    if type(point_model) == VoxelLattice:
                        suffix = ""
                    else:
                        suffix = "_WA"
                    point_model.save_as_nifti(
                        field_mags_full,
                        os.path.join(self.output_path, f"E_field_solution{suffix}.nii"),
                    )
                    point_model.save_as_nifti(
                        field_mags_full,
                        os.path.join(self.output_path, f"VTA_solution{suffix}.nii"),
                        binarize=True,
                        activation_threshold=activation_threshold,
                    )

            time_1 = time.time()
            timings["FieldExport"] = time_1 - time_0
            time_0 = time_1

        # save impedance at all frequencies to file!
        if compute_impedance:
            _logger.info("Saving impedance")
            df = pd.DataFrame(
                {
                    "freq": self.signal.frequencies,
                    "real": self.impedances.real,
                    "imag": self.impedances.imag,
                }
            )
            df.to_csv(os.path.join(self.output_path, "impedance.csv"), index=False)

        return timings

    @property
    def output_path(self) -> str:
        return self._output_path

    @output_path.setter
    def output_path(self, path: str) -> None:
        """Set the path to write output.

        Notes
        -----
        Creates directory if it doesn't exist.
        """
        if not os.path.exists(path):
            os.makedirs(path)
        self._output_path = path

    @property
    def conductivity_cf(self) -> ConductivityCF:
        return self._conductivity_cf

    @property
    def signal(self) -> FrequencyDomainSignal:
        return self._signal

    @property
    def current_controlled(self) -> bool:
        return self.signal.current_controlled

    @property
    def impedances(self) -> np.ndarray:
        return self._impedances

    @property
    def conductivity(self) -> ngsolve.CoefficientFunction:
        """Return conductivity of latest solution."""
        return self._sigma

    @property
    def is_complex(self) -> bool:
        """Return the state of the data type for spaces. True if complex,
        False otherwise.

        Returns
        -------
        bool
        """
        return self._complex

    @is_complex.setter
    def is_complex(self, value: bool) -> None:
        self._complex = value

    @property
    def model_geometry(self) -> ModelGeometry:
        return self._model_geometry

    @property
    def mesh(self) -> Mesh:
        return self._mesh

    @property
    def solver(self) -> Solver:
        return self._solver

    @property
    def contacts(self) -> Contacts:
        return self._contacts

    def get_scaled_active_contact_voltages(self, factor) -> dict:
        contact_voltages = {}
        for contact in self.contacts.active:
            contact_voltages[contact.name] = (
                self._base_contacts[contact.name].voltage * factor
            )
        return contact_voltages

    def update_contacts(
        self,
        voltages: Optional[dict] = None,
        currents: Optional[dict] = None,
        surface_impedances: Optional[dict] = None,
    ) -> None:
        """TODO document."""
        if surface_impedances is None:
            surface_impedances = {}
        if currents is None:
            currents = {}
        if voltages is None:
            voltages = {}
        if self.is_complex:
            self._contacts.voltages = voltages
            self._contacts.currents = currents
            self._contacts.surface_impedances = surface_impedances
            return
        self._contacts.voltages = np.real(voltages)
        self._contacts.currents = np.real(currents)
        self._contacts.surface_impedances = np.real(surface_impedances)

    @property
    def potential(self) -> ngsolve.GridFunction:
        return self._potential

    @property
    def frequency(self) -> float:
        return self._frequency

    def evaluate_potential_at_points(self, lattice: np.ndarray) -> np.ndarray:
        """Return electric potential at specifed 3-D coordinates.

        Parameters
        ----------
        lattice : Nx3 numpy.ndarray of lattice points

        Returns
        -------
        Nx1 numpy.ndarray

        Notes
        -----
        Make sure that points outside of the computational domain
        are filtered!
        """
        mesh = self.mesh.ngsolvemesh
        x, y, z = lattice.T
        pots = self.potential(mesh(x, y, z))
        return pots

    def evaluate_field_at_points(self, lattice: np.ndarray) -> np.ndarray:
        """Return electric field components at specifed 3-D coordinates.

        Parameters
        ----------
        lattice : Nx3 numpy.ndarray of lattice points

        Returns
        -------
        Nx3 numpy.ndarray

        TODO: filter points outside of the computational domain
        """
        mesh = self.mesh.ngsolvemesh
        x, y, z = lattice.T
        fields = self.electric_field(mesh(x, y, z))
        return fields

    @property
    def current_density(self) -> ngsolve.GridFunction:
        """Return current density in A/mm^2."""
        _logger.debug(f"Compute current density at frequency {self.frequency} Hz")
        # scale to account for mm as length unit (not yet contained in conductivity)
        return 1e-3 * self.conductivity * self.electric_field

    @property
    def electric_field(self) -> ngsolve.GridFunction:
        return -ngsolve.grad(self.potential)

    def compute_impedance(self) -> complex:
        """TODO document."""
        if len(self.contacts.active) == 2:
            mesh = self._mesh.ngsolvemesh
            # do not need to account for mm because of integration
            power = ngsolve.Integrate(self.electric_field * self.current_density, mesh)
            # TODO integrate surface impedance
            voltage = 0
            for idx, contact in enumerate(self.contacts.active):
                voltage += (-1) ** idx * contact.voltage
            return voltage / power
        else:
            raise NotImplementedError(
                "Impedance for more than two active contacts not yet supported"
            )

    def vtk_export(self) -> None:
        """Export all relevant properties to VTK."""
        ngmesh = self.mesh.ngsolvemesh
        # TODO add frequency to name
        FieldSolution(self.potential, "potential", ngmesh, self.is_complex).save(
            os.path.join(self.output_path, "potential")
        )

        FieldSolution(self.electric_field, "E-field", ngmesh, self.is_complex).save(
            os.path.join(self.output_path, "E-field")
        )

        FieldSolution(self.conductivity, "conductivity", ngmesh, self.is_complex).save(
            os.path.join(self.output_path, "conductivity")
        )

        FieldSolution(
            self.conductivity_cf.material_distribution(self.mesh),
            "material",
            ngmesh,
            False,
        ).save(os.path.join(self.output_path, "material"))

    def floating_values(self) -> dict:
        floating_voltages = {}
        for contact in self.contacts.floating:
            floating_voltages[contact.name] = contact.voltage
        return floating_voltages

    def flux_space(self) -> ngsolve.comp.HDiv:
        """Return a flux space on the mesh.

        Returns
        -------
        ngsolve.HDiv

        Notes
        -----
        The HDiv space is returned with a minimum order of 1.
        It is needed for the a-posteriori error estimator
        needed for adaptive mesh refinement.

        """
        return ngsolve.HDiv(
            mesh=self.mesh, order=max(1, self._order - 1), complex=self._complex
        )

    def h1_space(self, boundaries: List[str]) -> ngsolve.comp.H1:
        """Return a h1 space on the mesh.

        Parameters
        ----------
        boundaries : list of str
            List of boundary names.

        Returns
        -------
        ngsolve.H1
        """
        dirichlet = "|".join(boundary for boundary in boundaries)
        return ngsolve.H1(
            mesh=self.mesh.ngsolvemesh,
            order=self._order,
            dirichlet=dirichlet,
            complex=self._complex,
            wb_withedges=False,
        )

    def number_space(self) -> ngsolve.comp.NumberSpace:
        """Return a number space on the mesh.

        Returns
        -------
        ngsolve.NumberSpace
            Space with only one single (global) DOF.
        """
        return ngsolve.NumberSpace(
            mesh=self.mesh.ngsolvemesh, order=0, complex=self.is_complex
        )

    def surfacel2_space(self, boundaries: List[str]) -> ngsolve.comp.SurfaceL2:
        """Return a number SurfaceL2 on the mesh.

        Returns
        -------
        ngsolve.SurfaceL2

        Notes
        -----
        The SurfaceL2 space is returned with a minimum order of 1.
        It is needed to impose floating potentials.

        """
        dirichlet = "|".join(boundary for boundary in boundaries)
        return ngsolve.SurfaceL2(
            mesh=self.mesh.ngsolvemesh,
            order=max(1, self._order - 1),
            dirichlet=dirichlet,
            complex=self.is_complex,
        )

    def get_points_in_encapsulation_layer(self, points: np.ndarray) -> np.ndarray:
        """Return mask for points in encapsulation layer."""
        encap_cf = self.mesh.ngsolvemesh.RegionCF(
            ngsolve.VOL, {"EncapsulationLayer_*": 1.0}, default=0
        )
        mesh = self.mesh.ngsolvemesh
        x, y, z = points.T
        return np.isclose(encap_cf(mesh(x, y, z)), 1.0)

    def get_points_in_csf(self, points: np.ndarray) -> np.ndarray:
        """Return mask for points in CSF."""
        material_distribution = self.conductivity_cf.material_distribution(self.mesh)
        mesh = self.mesh.ngsolvemesh
        x, y, z = points.T
        return np.isclose(
            material_distribution(mesh(x, y, z)), self.conductivity_cf.materials["CSF"]
        )
