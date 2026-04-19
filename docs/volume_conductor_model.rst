Volume conductor model
======================

The volume conductor model is the numerical core of OSS-DBSv2. It combines the
geometry, tissue properties, stimulation signal, and solver configuration to
compute electric potentials and fields in the simulated domain.

Role in the workflow
--------------------

Most standalone and Lead-DBS-driven runs eventually reach the same core steps:

1. build the geometry and mesh
2. assign conductivity information from MRI and optional DTI data
3. define boundary conditions at electrode contacts and outer surfaces
4. solve the resulting finite-element system
5. export fields, impedance estimates, and pointwise analyses

Floating and non-floating contacts
----------------------------------

OSS-DBSv2 distinguishes between several conductor variants:

- ``VolumeConductorNonFloating`` for standard active and passive contact setups
- ``VolumeConductorFloating`` when floating contacts are present
- ``VolumeConductorFloatingImpedance`` when floating contacts also carry surface
  impedance information

These variants share the same general workflow but differ in how contact
boundary conditions are enforced.

.. _stimulation-modes:

Stimulation modes
-----------------

OSS-DBSv2 supports both voltage-controlled and current-controlled stimulation.
The ``CurrentControlled`` flag in ``StimulationSignal`` selects the mode.
Depending on the mode, contacts can be **active** (carrying a Dirichlet
boundary condition), **floating** (equipotential via Lagrange multiplier), or
**floating with surface impedance** (coupled via a Robin boundary condition).
Not every combination is valid. The table below summarises the supported cases.

Voltage-controlled stimulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In voltage-controlled mode (``CurrentControlled: false``), the ``Voltage[V]``
value on each active contact is imposed directly as a Dirichlet boundary
condition.

- **Two or more active contacts**: each contact is set to its prescribed
  voltage. At least one should be at 0 V to serve as ground. This is the
  standard bipolar or multipolar voltage-controlled configuration.
- **One active contact**: a single Dirichlet condition is imposed. Floating
  contacts may carry the return current. Without a second active contact or
  floating contacts, the solution is trivially zero everywhere.
- **Floating contacts (plain)**: floating contacts are constrained to be
  equipotential surfaces via a Lagrange multiplier in the
  ``VolumeConductorFloating`` formulation. The potential on each floating
  contact is determined by current conservation.
- **Floating contacts with surface impedance**: the
  ``VolumeConductorFloatingImpedance`` formulation replaces the equipotential
  constraint with a Robin boundary condition that models the electrode–tissue
  interface impedance. The contact potential is a free variable coupled to the
  tissue potential through the admittance :math:`1/Z_s`.

Current-controlled stimulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In current-controlled mode (``CurrentControlled: true``), each contact
specifies a ``Current[A]`` value. The sum of all currents (active and floating)
must be zero (Kirchhoff's current law). Internally, OSS-DBSv2 uses an
admittance-matrix superposition to translate prescribed currents into the
potentials that the FEM solver needs.

- **Two active contacts (bipolar)**: one contact must be grounded
  (``Voltage[V]: 0``). The other is assigned a pseudo-voltage for an
  intermediate Dirichlet solve; the final solution is scaled so that the
  prescribed current flows. Floating contacts may be present in either the
  plain or surface-impedance variant.
- **One active contact (monopolar with floating)**: the single active contact
  must be grounded (``Voltage[V]: 0``). All other contacts are floating with
  prescribed currents. In the plain floating case, each floating contact is
  equipotential. In the surface-impedance case, Robin BCs are used instead.
- **No active contacts (all floating with surface impedance)**: this is only
  valid in the ``VolumeConductorFloatingImpedance`` formulation. Because there
  is no Dirichlet boundary to fix the voltage gauge, a Lagrange multiplier
  constrains the sum of floating potentials to zero
  (:math:`\sum_k u_k = 0`). Currents are imposed on the floating contacts
  directly.
- **More than two active contacts**: currently not supported in
  current-controlled mode. Use at most one grounded active contact plus
  floating contacts to model multipolar configurations with more than two
  electrodes.

Constraints and validation
^^^^^^^^^^^^^^^^^^^^^^^^^^

The code enforces several consistency checks before solving:

- In current-controlled mode, the sum of all prescribed currents must be zero.
- In multipolar current-controlled mode with plain floating contacts, exactly
  one active contact must be grounded.
- The ``FloatingImpedance`` formulation requires *every* floating contact to
  carry a surface impedance model. Mixed setups (some floating contacts with
  impedance, some without) are rejected at construction time.
- Active contacts must not carry a surface impedance model in the
  ``FloatingImpedance`` formulation.

Impedance computation
^^^^^^^^^^^^^^^^^^^^^

When ``ComputeImpedance`` is enabled, OSS-DBSv2 computes the complex impedance
at each frequency. For current-controlled setups with multiple contacts, a full
admittance matrix is assembled by superposition and then inverted. For setups
with surface impedance, the dissipation across the electrode–tissue interface is
included in the power balance so that the impedance estimate accounts for the
interface losses.

.. _surface-impedance:

Surface impedance
-----------------

Real DBS electrodes do not make perfect electrical contact with surrounding
tissue. A thin electrochemical interface forms at the metal–tissue boundary,
adding a frequency-dependent impedance that affects current flow and measured
impedance values. OSS-DBSv2 models this interface as a **surface impedance**
applied as a Robin boundary condition on the contact surface.

Physical background
^^^^^^^^^^^^^^^^^^^

The interface is often described by an equivalent circuit (Randles-type model).
Common elements include a resistor for charge-transfer resistance and a
constant-phase element (CPE) for the double-layer capacitance. The surface
impedance :math:`Z_s` relates the potential jump across the interface to the
normal current density:

.. math::

   \mathbf{j} \cdot \mathbf{n} = \frac{1}{Z_s}\,(u - u_{\text{contact}})

where :math:`u` is the tissue-side potential and :math:`u_{\text{contact}}` is
the contact potential. The surface impedance has units of
:math:`\Omega \cdot \text{mm}^2` (impedance times contact area) because the FEM
geometry uses millimetres.

JSON configuration
^^^^^^^^^^^^^^^^^^

Surface impedance is set per contact in the ``Electrodes`` section:

.. code-block:: json

   {
     "Contact_ID": 1,
     "Active": true,
     "Voltage[V]": 1.0,
     "Floating": false,
     "SurfaceImpedance": {
       "Model": "R",
       "Parameters": {"R": 500.0}
     }
   }

The ``Model`` string and ``Parameters`` dictionary are passed to the
`impedancefitter <https://impedancefitter.readthedocs.io/>`_ library, which
evaluates the equivalent circuit at each frequency. This means any model
supported by impedancefitter can be used.

Available models
^^^^^^^^^^^^^^^^

Commonly used models include:

- ``R`` — pure resistance. Parameters: ``{"R": <value in Ohm>}``.
- ``RC`` — resistor in parallel with a capacitor.
  Parameters: ``{"R": <Ohm>, "C": <Farad>}``.
- ``CPE_dl`` — constant-phase element for the double layer.
  Parameters: ``{"dl_k": <value>, "dl_alpha": <exponent 0–1>}``.
  This is the model used by Lempka et al. (2009) for DBS electrode interfaces.

At zero frequency (DC), capacitive and CPE elements have infinite impedance.
Only the resistive path contributes, so a pure CPE model yields no interface
effect at DC.

Where surface impedance can be applied
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- **On active (nonfloating) contacts**: the ``VolumeConductorNonFloating``
  formulation adds a Robin term for each active contact that carries a
  ``SurfaceImpedance`` entry. The Dirichlet voltage is replaced by a Robin
  condition coupling the tissue potential to the prescribed contact voltage.
- **On floating contacts**: the ``VolumeConductorFloatingImpedance``
  formulation couples each floating contact to the tissue via the Robin
  condition. The contact potential becomes a free variable solved for by the
  FEM system. Every floating contact must carry a surface impedance model in
  this formulation.

See :ref:`stimulation-modes` for the full list of valid contact configurations.

Effect on impedance estimates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When ``ComputeImpedance`` is enabled, the power balance includes both the
volume dissipation (tissue) and the interface dissipation (surface impedance
layer). Without accounting for the interface term, the computed impedance would
overestimate the true value because part of the dissipated power is missing.

Mesh and linear algebra
-----------------------

The volume conductor stack also includes:

- conductivity coefficient functions derived from the material model
- iterative and direct solvers
- several preconditioner strategies for larger problems

.. _solver-guidance:

Solver selection
^^^^^^^^^^^^^^^^

OSS-DBSv2 provides three solver types:

- **CG** (Conjugate Gradient) — the default. Efficient for symmetric positive
  definite systems, which covers most real-valued (non-EQS) volume conductor
  problems. Use this as the starting point.
- **GMRES** — required for non-symmetric or complex-valued systems. When
  ``EQSMode`` is enabled, the system matrix becomes complex and CG is no
  longer applicable; switch to GMRES. GMRES also handles the
  ``FloatingImpedance`` formulation where the Lagrange multiplier and Robin
  terms can break strict symmetry.
- **Direct** — factorises the system matrix with PARDISO. Robust and
  deterministic but memory-intensive. Practical for small to medium problems
  (up to ~200k DOFs) or as a reference for verifying iterative solver results.

Preconditioner selection
^^^^^^^^^^^^^^^^^^^^^^^^

The preconditioner accelerates convergence of the iterative solvers (CG and
GMRES). Available options:

- **bddc** (default) — Balancing Domain Decomposition by Constraints. The most
  efficient option for standard (nonfloating) problems, especially at higher
  FEM orders. Not compatible with the ``FloatingImpedance`` formulation — the
  code automatically falls back to ``local`` in that case.
- **local** — NGSolve's built-in Jacobi (diagonal) preconditioner. Works with
  all formulations. Sufficient for moderate problem sizes. If the solver
  produces unexpected results with ``local``, try ``customized_local`` instead.
- **customized_local** — a custom Jacobi implementation that explicitly
  constructs the inverse diagonal. More robust than the native ``local``
  preconditioner for problems with mixed spaces (e.g. FloatingImpedance with
  Lagrange multipliers). Use this if ``local`` gives warnings about non-finite
  values.
- **h1amg** — Algebraic multigrid. Can be faster than ``local`` for very large
  problems but is less tested in this codebase.
- **multigrid** — Geometric multigrid. Requires a mesh hierarchy; not commonly
  used in practice.
- **direct** — Uses a direct solve as preconditioner. Expensive but can be
  useful for debugging convergence problems.

Recommended configurations
^^^^^^^^^^^^^^^^^^^^^^^^^^

+-----------------------------+---------------+--------------------+
| Scenario                    | Solver        | Preconditioner     |
+=============================+===============+====================+
| Standard (real-valued)      | CG            | bddc               |
+-----------------------------+---------------+--------------------+
| EQS mode (complex)          | GMRES         | bddc               |
+-----------------------------+---------------+--------------------+
| FloatingImpedance           | CG or GMRES   | local              |
+-----------------------------+---------------+--------------------+
| FloatingImpedance + EQS     | GMRES         | customized_local   |
+-----------------------------+---------------+--------------------+
| Small problem / debugging   | Direct        | (not applicable)   |
+-----------------------------+---------------+--------------------+

Convergence parameters
^^^^^^^^^^^^^^^^^^^^^^

- ``MaximumSteps`` (default: 10000) — upper limit on Krylov iterations. If
  the solver hits this limit, it raises a ``RuntimeError``. Increase this
  value for very large problems or tight precision targets.
- ``Precision`` (default: 1e-12) — relative residual tolerance. For most
  applications, 1e-8 to 1e-10 is sufficient. Very tight values (1e-12 and
  below) increase iteration count without visible improvement in the solution.
- ``PrintRates`` — when ``true``, the solver prints the residual norm at each
  iteration. Useful for monitoring convergence.

Troubleshooting
^^^^^^^^^^^^^^^

- **Solver does not converge**: try increasing ``MaximumSteps`` or relaxing
  ``Precision``. If that does not help, switch from ``local`` to
  ``customized_local`` or from CG to GMRES.
- **Warning about non-finite preconditioner values**: switch from ``local`` to
  ``customized_local``. This is most common with the ``FloatingImpedance``
  formulation.
- **BDDC fails with FloatingImpedance**: this is expected — the code
  automatically switches to ``local``. No action needed.
- **Very slow convergence**: check that the mesh is not excessively refined in
  regions far from the electrode. Consider HP refinement instead of uniform
  mesh refinement for better cost-accuracy tradeoff.

.. _mesh-refinement:

Mesh refinement
---------------

Mesh quality strongly affects the accuracy of the FEM solution, especially near
electrode contacts where the field gradient is steepest. OSS-DBSv2 provides
several levels of mesh control.

Global mesh hypothesis
^^^^^^^^^^^^^^^^^^^^^^

The ``MeshingHypothesis`` block sets the baseline element size:

.. code-block:: json

   "MeshingHypothesis": {
     "Type": "Moderate"
   }

Supported types are ``VeryCoarse``, ``Coarse``, ``Moderate``, ``Fine``,
``VeryFine``, and ``Default``. A ``Custom`` type is also available, in which
case a ``CustomParameters`` dictionary is passed directly to Netgen.

Optional parameters can be added regardless of type:

- ``MaxMeshSize`` — global upper bound on element size (in mm)
- ``CurvatureSafety`` — controls refinement near curved surfaces
- ``Grading`` — mesh grading factor (higher = faster size transition)
- ``MeshSizeFilename`` — path to a Netgen mesh-size file that specifies
  target element sizes at individual points.  This is especially useful
  for refining the mesh along neuron trajectories in pathway activation
  studies.  The file can be generated from a ``PathwayPointModel`` with
  :py:meth:`~ossdbs.point_analysis.point_model.PathwayPointModel.write_netgen_meshsize_file`:

  .. code-block:: python

     from ossdbs.api import generate_point_models, load_images

     mri_image, _ = load_images(settings)
     point_models = generate_point_models(settings)
     pw = point_models[0]  # PathwayPointModel
     pw.write_netgen_meshsize_file(
         meshsize=min(mri_image.voxel_sizes),
         filename="meshsizes.txt",
     )

  Then reference the file in the JSON input:

  .. code-block:: json

     "MeshingHypothesis": {
       "Type": "Fine",
       "MeshSizeFilename": "meshsizes.txt"
     }

Local mesh sizes
^^^^^^^^^^^^^^^^

The ``MeshSize`` block allows targeted refinement on named geometry entities:

.. code-block:: json

   "MeshSize": {
     "Edges": {"E1C1": 0.05},
     "Faces": {"E1C1": 0.1, "E1C2": 0.1},
     "Volumes": {"EncapsulationLayer_1": 0.2}
   }

Keys are boundary or region names as they appear in the geometry (e.g.
``E1C1`` for electrode 1, contact 1). Contact-level ``MaxMeshSize`` and
``MaxMeshSizeEdge`` in the ``Contacts`` JSON block provide a shorthand for
the same mechanism.

Mesh element order
^^^^^^^^^^^^^^^^^^

``MeshElementOrder`` (default: 2) sets the polynomial order of the curved mesh
geometry. Higher orders better approximate the true electrode surface but
increase mesh generation cost. This is independent of the FEM polynomial order
set by ``FEMOrder``.

HP refinement
^^^^^^^^^^^^^

HP refinement combines geometric (h) and polynomial (p) refinement near
boundaries, concentrating degrees of freedom where the solution varies most:

.. code-block:: json

   "HPRefinement": {
     "Active": true,
     "Levels": 2,
     "Factor": 0.125
   }

- ``Levels`` — number of refinement layers added near boundaries
- ``Factor`` — geometric grading factor; smaller values concentrate elements
  more tightly at the boundary

HP refinement provides the best cost-accuracy tradeoff for DBS simulations,
typically reaching near-best accuracy at 60–120k degrees of freedom. It is
applied **after** any material-based bisection refinement so that the
specialised HP element types are not disrupted by subsequent bisection steps.

Material-based mesh refinement
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When tissue boundaries are poorly resolved by the initial mesh, bisection-based
refinement can be applied. This is controlled by the
``MaterialRefinementSteps`` parameter (default: 0) and uniformly bisects
elements that straddle material interfaces.

Mesh reuse
^^^^^^^^^^

For parameter studies where only the stimulation changes but the geometry stays
fixed, mesh generation can be skipped:

.. code-block:: json

   "Mesh": {
     "SaveMesh": true,
     "LoadMesh": true,
     "LoadPath": "path/to/saved_mesh.vol.gz"
   }

When loading a mesh, HP refinement is re-applied automatically if configured.

Practical guidance
^^^^^^^^^^^^^^^^^^

- Start with ``Moderate`` hypothesis and ``HPRefinement`` with 2 levels.
- Check convergence by comparing results across refinement levels before
  drawing scientific conclusions.
- Use local face/edge mesh sizes on contacts when the global hypothesis does
  not provide enough resolution near the electrode.
- For large studies, save the base mesh and reload it to avoid repeated mesh
  generation.

Related pages
-------------

- :doc:`input_settings` — JSON reference for mesh, solver, and stimulation
  settings
- :doc:`brain_geometry` — geometry construction and encapsulation layers
- :doc:`electrodes` — electrode models and custom parameters
- :doc:`materials` — dielectric models and tissue conductivity
- :doc:`stimulation_signals` — signal types and frequency-domain settings
- :doc:`point_analysis` — field evaluation on lattices and pathways
- :doc:`python_api` — scripting the full pipeline from Python

API reference
-------------

Volume conductor model
^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: ossdbs.fem.volume_conductor.volume_conductor_model
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: ossdbs.fem.volume_conductor.nonfloating
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: ossdbs.fem.volume_conductor.floating
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: ossdbs.fem.volume_conductor.floating_impedance
    :members:
    :undoc-members:
    :show-inheritance:

Material parameters
^^^^^^^^^^^^^^^^^^^

.. automodule:: ossdbs.fem.volume_conductor.conductivity
    :members:
    :undoc-members:
    :show-inheritance:


Mesh
^^^^

.. automodule:: ossdbs.fem.mesh
    :members:
    :undoc-members:
    :show-inheritance:

Solver
^^^^^^

.. automodule:: ossdbs.fem.solver
    :members:
    :undoc-members:
    :show-inheritance:


Preconditioner
^^^^^^^^^^^^^^

.. automodule:: ossdbs.fem.preconditioner
    :members:
    :undoc-members:
    :show-inheritance:
