# Runtime benchmark: best convergence-study strategy

Records a per-phase runtime breakdown of the mesh-refinement strategy that
the convergence study identified as the best cost/accuracy tradeoff, so the
same workload can be compared across machines and operating systems.

## What is measured

**Strategy** — identical for both workloads, pinned in the run scripts and
deliberately not imported from a study script, so recorded results stay
comparable as the studies evolve:

| Setting | Value |
| --- | --- |
| Meshing hypothesis | `Fine` |
| HP refinement | active, 2 levels, factor 0.125 |
| Material refinement | 1 step |
| Adaptive refinement | off |

**Two workloads**, sharing that strategy so a pair of records isolates the
effect of the workload rather than of the refinement:

| | `run_benchmark.py` (PAM) | `run_benchmark_vta.py` (VTA) |
| --- | --- | --- |
| Case | `PAM_3/BostonVerciseDirected`, monopolar | `VTA/BostonScientificVerciseDirected`, case grounding |
| Drive | current-controlled | voltage-controlled, contact 2 at 1 V |
| Signal | Rectangle, 12 octave bands | Multisine, single frequency |
| Evaluation points | 75935 pathway points | 864000 lattice points |
| Time domain | reconstructed | skipped (single frequency) |
| NEURON stage | yes | no |
| Runs on Windows | no (NEURON) | yes |
| Results directory | `results/` | `results_vta/` |

Both are small enough to run on a laptop — a benchmark is only useful if it
gets run on more than one machine.

**Phases**, from `run_report.json` (written by `main_run`) and
`VCM_report.json` (written by the volume conductor):

| Phase | Meaning |
| --- | --- |
| `load_images` | reading the MRI |
| `geometry` | electrode + brain OCC geometry construction |
| `contact_properties`, `dielectric_model`, `conductivity` | model setup |
| `meshing_and_refinement_derived` | mesh generation, material + HP refinement, FE space setup |
| `fem_solve` | summed `ComputeSolution` over all frequencies |
| `point_model_copy` | locating the evaluation points on the mesh and writing the solution to them |
| `time_reconstruction` | inverse FFT to the time domain (zero for VTA) |
| `field_export` | field output |
| `pam_total` | the NEURON pathway activation stage (absent for VTA) |
| `vta_volume_mm3` | VTA only: not a timing, but confirms two machines solved the same problem |

`meshing_and_refinement_derived` is obtained by subtracting the volume
conductor's own timings from the `VolumeConductor` phase, because mesh
generation is not instrumented separately. It is worth isolating: Netgen
meshing is largely single-threaded while the FEM solve is not, so the two
scale very differently across machines.

## Requirements

The PAM benchmark needs Linux or macOS: it times FEM **and** PAM, and PAM
needs NEURON, which is unavailable on Windows — `run_benchmark.py` refuses to
run there rather than silently reporting a partial number. The VTA benchmark
has no NEURON stage and runs anywhere, which makes it the one to use when
comparing across operating systems.

No other input is needed. The MRI input `../PAM_3/segmask.nii.gz` is tracked
in the repository despite the `*.nii.gz` ignore rule, so the benchmark runs
from a clean clone — a benchmark nobody else can run is of no use. Note it
is *not* interchangeable with `../PAM/segmask.nii.gz`: the two share a grid
but differ in 0.27 % of voxels, which changes the conductivity field and
therefore the mesh. The VTA case uses the same file: `VTA/segmask.nii.gz` is
byte-identical to it, so the VTA benchmark points at the tracked copy rather
than requiring a second 1.1 MB binary.

## Running

```bash
uv run python run_benchmark.py --dry-run          # show machine context, run nothing
uv run python run_benchmark.py                    # one PAM run
uv run python run_benchmark.py --repeats 3        # keep the fastest of three
uv run python run_benchmark.py --label hpc-node-01

uv run python run_benchmark_vta.py --repeats 3    # same flags, VTA workload
```

Each invocation writes `<results dir>/<label>-<timestamp>.json` containing the
full machine context, the pinned strategy, every run, and the fastest one —
`results/` for PAM, `results_vta/` for VTA.

Then aggregate, one directory at a time:

```bash
uv run python collect_results.py                            # PAM records
uv run python collect_results.py --results-dir results_vta  # VTA records
```

which writes `benchmark_summary.csv` and prints a Markdown table. Do not mix
the two directories: the workloads solve different problems, and comparing
their timings is meaningless.

## Contributing a result

1. Run on an otherwise idle machine — background load distorts the solve.
2. Use `--repeats 3` if you can; the fastest run is the least noisy estimate.
3. Note your threading setup. `OMP_NUM_THREADS` and friends are recorded
   automatically, and an unset variable is stored as `null`, which is *not*
   the same as `1` — it means the library chose its own default.
4. Commit the `results/*.json` or `results_vta/*.json` file. Do not edit it by
   hand.

`collect_results.py` warns when records differ in `benchmark_version`, case,
or DOF count. Those records describe different work and must not be compared —
raise `BENCHMARK_VERSION` in the corresponding run script whenever its pinned
strategy or case changes, so old records are flagged rather than silently
mixed in.

## Interpreting the numbers

Compare `dofs` and `elements` first — and `vta_volume_mm3` for VTA. If they
differ between machines, the runs solved different problems and the timings
are not comparable; this usually means a different Netgen version meshed the
geometry differently.

Note what dominates. In the PAM workload `point_model_copy` is the largest
phase, several times the linear solve, so the figures mostly reflect point
evaluation and NEURON rather than solver throughput. Read a machine
comparison as "how fast does this machine run this pipeline", not "how fast
does it do FEM".
