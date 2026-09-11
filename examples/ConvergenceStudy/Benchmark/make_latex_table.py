"""Render the benchmark records as a LaTeX table.

Emits one block per machine, with a header row naming the CPU, cores and
memory, followed by one row per workload. Adding a machine means dropping
its records into the results directories and re-running this script -- the
numbers are never typed by hand.

Requires booktabs and siunitx.

Usage:
    python make_latex_table.py            # writes benchmark_table.tex
    python make_latex_table.py --stdout   # print instead
"""

import argparse
import glob
import json
import os
import re

# (results directory, label used in the table), in the order rows appear
WORKLOADS = [
    ("results", "PAM"),
    ("results_vta", "VTA (lattice)"),
    ("results_vta_ngsolve", "VTA (NGSolve only)"),
]

# (record key, column header); None renders as an em dash
COLUMNS = [
    ("dofs", "{DOFs}"),
    ("meshing_and_refinement_derived", "{Mesh}"),
    ("fem_solve", "{Solve}"),
    ("point_model_copy", "{Points}"),
    ("time_reconstruction", "{Recon.}"),
    ("field_export", "{Export}"),
    ("pam_total", "{PAM}"),
    ("wall_total", "{Total}"),
]


def load_records():
    """Return {machine key: {workload label: record}} plus machine contexts."""
    machines = {}
    contexts = {}
    for results_dir, workload in WORKLOADS:
        for path in sorted(glob.glob(os.path.join(results_dir, "*.json"))):
            with open(path) as fp:
                record = json.load(fp)
            key = record["label"]
            machines.setdefault(key, {})[workload] = record["best"]
            contexts.setdefault(key, record["machine"])
    return machines, contexts


def machine_caption(machine):
    """One-line hardware description used as the block header."""
    cpu = machine.get("cpu_model") or "unknown CPU"
    # vendor strings carry trademark noise and an inconsistently spaced clock
    for noise in ("(R)", "(TM)", "(r)", "(tm)", " CPU"):
        cpu = cpu.replace(noise, "")
    cpu = " ".join(cpu.split())
    cpu = re.sub(r"@\s*([\d.]+)\s*GHz", r"@ \\SI{\1}{\\giga\\hertz}", cpu)
    cores = machine.get("cores_physical")
    threads = machine.get("cores_logical")
    memory = machine.get("memory_gb")
    parts = [cpu]
    if cores and threads:
        parts.append(f"{cores} cores / {threads} threads")
    elif threads:
        parts.append(f"{threads} threads")
    if memory:
        parts.append(f"\\SI{{{memory:.0f}}}{{\\giga\\byte}} RAM")
    if machine.get("os"):
        parts.append(machine["os"])
    return ", ".join(parts)


def fmt(key, value):
    """Format one cell."""
    if value is None:
        return "{--}"
    if key == "dofs":
        return f"{int(value)}"
    return f"{value:.1f}"


def render(machines, contexts):
    """Return the LaTeX table as a string."""
    align = "l" + " S" * len(COLUMNS)
    lines = [
        r"% Requires \usepackage{booktabs} and \usepackage{siunitx}",
        r"\begin{table}[htbp]",
        r"  \centering",
        (
            r"  \caption{Runtime breakdown of the best convergence-study "
            r"strategy (\texttt{Fine} mesh, HP refinement with two levels and "
            r"factor $0.125$, one material refinement step). Times in "
            r"\si{\second}, fastest of three runs. \emph{Points} is the "
            r"evaluation of the solution at the pathway or lattice points, "
            r"\emph{Recon.} the reconstruction of the time domain.}"
        ),
        r"  \label{tab:benchmark}",
        rf"  \begin{{tabular}}{{{align}}}",
        r"    \toprule",
        "    Workload & " + " & ".join(header for _, header in COLUMNS) + r" \\",
        r"    \midrule",
    ]

    machine_keys = sorted(machines)
    for index, key in enumerate(machine_keys):
        if index:
            lines.append(r"    \addlinespace")
        caption = machine_caption(contexts[key])
        span = len(COLUMNS) + 1
        lines.append(rf"    \multicolumn{{{span}}}{{l}}{{\itshape {caption}}} \\")
        for _, workload in WORKLOADS:
            best = machines[key].get(workload)
            if best is None:
                continue
            cells = [fmt(col, best.get(col)) for col, _ in COLUMNS]
            lines.append(f"    {workload} & " + " & ".join(cells) + r" \\")

    lines += [r"    \bottomrule", r"  \end{tabular}", r"\end{table}"]
    return "\n".join(lines) + "\n"


def main():
    """Write or print the LaTeX table."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stdout", action="store_true", help="Print instead of writing the file."
    )
    parser.add_argument(
        "--output", default="benchmark_table.tex", help="Output file name."
    )
    args = parser.parse_args()

    machines, contexts = load_records()
    if not machines:
        print("No benchmark records found.")
        return

    table = render(machines, contexts)
    if args.stdout:
        print(table)
        return
    with open(args.output, "w") as fp:
        fp.write(table)
    print(f"Wrote {args.output} ({len(machines)} machine(s))")


if __name__ == "__main__":
    main()
