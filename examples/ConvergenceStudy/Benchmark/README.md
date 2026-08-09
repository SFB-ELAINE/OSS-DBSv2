# Runtime benchmark: best convergence-study strategy

Records a per-phase runtime breakdown of the mesh-refinement strategy that
the convergence study identified as the best cost/accuracy tradeoff, so the
same workload can be compared across machines and operating systems.

## What is measured

**Strategy** (pinned in `run_benchmark.py`, deliberately not imported from a
study script so that recorded results stay comparable as the studies evolve):

| Setting | Value |
| --- | --- |
| Meshing hypothesis | `Fine` |
| HP refinement | active, 2 levels, factor 0.125 |
| Material refinement | 1 step |
| Adaptive refinement | off |

**Case**: `PAM_3/BostonVerciseDirected`, single monopolar current-controlled
protocol. Chosen to be small enough to run on a laptop — a benchmark is only
useful if it gets run on more than one machine.

**Phases**, from `run_report.json` (written by `main_run`) and
`VCM_report.json` (written by the volume conductor):

| Phase | Meaning |
| --- | --- |
| `load_images` | reading the MRI |
| `geometry` | electrode + brain OCC geometry construction |
| `contact_properties`, `dielectric_model`, `conductivity` | model setup |
| `meshing_and_refinement_derived` | mesh generation, material + HP refinement, FE space setup |
| `fem_solve` | summed `ComputeSolution` over all frequencies |
| `point_model_copy` | writing the solution to the pathway points |
| `time_reconstruction` | inverse FFT to the time domain |
| `field_export` | field output |
| `pam_total` | the NEURON pathway activation stage |

`meshing_and_refinement_derived` is obtained by subtracting the volume
conductor's own timings from the `VolumeConductor` phase, because mesh
generation is not instrumented separately. It is worth isolating: Netgen
meshing is largely single-threaded while the FEM solve is not, so the two
scale very differently across machines.

## Requirements

Linux or macOS. The benchmark times FEM **and** PAM, and PAM needs NEURON,
which is unavailable on Windows — `run_benchmark.py` refuses to run there
rather than silently reporting a partial number.

The MRI input `../PAM_3/segmask.nii.gz` is required. It is not in the
repository (`*.nii.gz` is gitignored), so copy it in before the first run.

## Running

```bash
uv run python run_benchmark.py --dry-run          # show machine context, run nothing
uv run python run_benchmark.py                    # one run
uv run python run_benchmark.py --repeats 3        # keep the fastest of three
uv run python run_benchmark.py --label hpc-node-01
```

Each invocation writes `results/<label>-<timestamp>.json` containing the full
machine context, the pinned strategy, every run, and the fastest one.

Then aggregate:

```bash
uv run python collect_results.py
```

which writes `benchmark_summary.csv` and prints a Markdown table.

## Contributing a result

1. Run on an otherwise idle machine — background load distorts the solve.
2. Use `--repeats 3` if you can; the fastest run is the least noisy estimate.
3. Note your threading setup. `OMP_NUM_THREADS` and friends are recorded
   automatically, and an unset variable is stored as `null`, which is *not*
   the same as `1` — it means the library chose its own default.
4. Commit the `results/*.json` file. Do not edit it by hand.

`collect_results.py` warns when records differ in `benchmark_version` or DOF
count. Those records describe different work and must not be compared —
raise `BENCHMARK_VERSION` in `run_benchmark.py` whenever the pinned strategy
or the case changes, so old records are flagged rather than silently mixed in.

## Interpreting the numbers

Compare `dofs` and `elements` first. If they differ between machines, the
runs solved different problems and the timings are not comparable — this
usually means a different Netgen version meshed the geometry differently.
