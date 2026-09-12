"""Evaluate StimSets VTA convergence.

For each mesh-refinement strategy, this reads the per-contact unit
E-field solutions (`Results_VTA_<strategy>E1C<i>/E_field_Lattice.csv`),
superimposes them for every current protocol in `Current_protocols_0.csv`
(the same superposition Lead-DBS performs in `ea_get_ossdbs_StimSets_VTA`),
and computes a VTA per protocol.

Metrics, per strategy and averaged over all protocols:
  * VTA volume (mm^3) for |E| > threshold,
  * relative volume error vs. the `best` strategy,
  * Dice overlap of the binary VTA mask vs. the `best` strategy.

The lattice point order is identical across strategies and contacts, so
superposition and Dice are simple element-wise operations.

Output: `vta_results_summary.csv` (one row per strategy) and
`vta_volumes_per_protocol.csv` (strategy x protocol volume matrix).
"""

import json
import os

import numpy as np
import pandas as pd

CONFIG_FILE = "oss-dbs_parameters_vta.json"
STIM_SETS_FILE = "OSS_sim_files_rh/Current_protocols_0.csv"
REFERENCE = "best"

# VTA activation threshold (V/m). Thresholding is done here (and in
# Lead-DBS) rather than on the FEM mesh, so the OSS-DBS config leaves
# ActivationThresholdVTA[V-per-m] = null (the on-mesh integration cannot
# handle the complex multi-frequency StimSets field).
THRESHOLD_VPM = 200.0

# Strategy result-directory prefixes (without the per-contact E1C<i> suffix).
STRATEGIES = [
    "default",
    "fine",
    "very_fine",
    "material_refinement",
    "edge_refinement",
    "fine_edge_refinement",
    "edge_voxel_refinement",
    "edge_single_material_refinement",
    "edge_double_material_refinement",
    "very_fine_edge_refinement",
    "fine_edge_single_material_refinement",
    "fine_edge_double_material_refinement",
    "hp_refinement",
    "hp_material_refinement",
    "best",
]


def load_spacing():
    """Read the lattice point spacing (mm) from the config."""
    with open(CONFIG_FILE) as fp:
        cfg = json.load(fp)
    return float(cfg["PointModel"]["Lattice"]["PointDistance[mm]"])


def load_stim_protocols():
    """Load current protocols (mA -> A, NaN -> 0). Shape (n_protocols, n_contacts)."""
    df = pd.read_csv(STIM_SETS_FILE)
    arr = df.to_numpy(dtype=float)
    arr = np.nan_to_num(arr, nan=0.0)
    return arr * 1e-3  # mA -> A


def _coord_keys(df):
    """Hashable per-point keys from lattice coordinates (rounded to um)."""
    xyz = np.round(df[["x-pt", "y-pt", "z-pt"]].to_numpy(dtype=float), 3)
    # encode as one string per row for fast set/dict alignment
    return np.array([f"{x:.3f},{y:.3f},{z:.3f}" for x, y, z in xyz])


def load_unit_fields(strategy, n_contacts):
    """Load per-contact unit E-field components for one strategy.

    Returns (keys, Ex, Ey, Ez): coordinate keys (n_points,) and field
    components each (n_points, n_contacts) in V/mm, or None if any
    contact result is missing. Different meshes keep slightly different
    boundary lattice points, so strategies are aligned on `keys`.
    """
    ex_list, ey_list, ez_list = [], [], []
    keys = None
    for contact_i in range(1, n_contacts + 1):
        csv_path = os.path.join(
            f"Results_VTA_{strategy}E1C{contact_i}", "E_field_Lattice.csv"
        )
        if not os.path.isfile(csv_path):
            print(f"  missing: {csv_path}")
            return None
        df = pd.read_csv(csv_path)
        if keys is None:
            keys = _coord_keys(df)
        ex_list.append(df["x-field"].to_numpy(dtype=float))
        ey_list.append(df["y-field"].to_numpy(dtype=float))
        ez_list.append(df["z-field"].to_numpy(dtype=float))
    return (
        keys,
        np.column_stack(ex_list),
        np.column_stack(ey_list),
        np.column_stack(ez_list),
    )


def compute_vta_masks(unit_fields, stim_protocols, threshold_vpm):
    """Binary VTA mask per protocol. Returns array (n_protocols, n_points)."""
    _keys, ex, ey, ez = unit_fields
    n_protocols = stim_protocols.shape[0]
    masks = np.empty((n_protocols, ex.shape[0]), dtype=bool)
    for p in range(n_protocols):
        s = stim_protocols[p]  # (n_contacts,) in A
        # superimpose unit solutions (V/mm)
        ex_tot = ex @ s
        ey_tot = ey @ s
        ez_tot = ez @ s
        magn_vpm = np.sqrt(ex_tot**2 + ey_tot**2 + ez_tot**2) * 1000.0  # V/mm -> V/m
        masks[p] = magn_vpm > threshold_vpm
    return masks


def read_dof(strategy):
    """Best-effort DOF/element count from the first contact's VCM report."""
    report = os.path.join(f"Results_VTA_{strategy}E1C1", "VCM_report.json")
    try:
        with open(report) as fp:
            d = json.load(fp)
        return d.get("DOF"), d.get("Elements")
    except (OSError, json.JSONDecodeError):
        return None, None


def read_time(strategy, n_contacts):
    """Total ComputeSolution time (s) summed over all contact solves.

    Mirrors the cost measure used by the PAM overview (one StimSets VTA
    estimate needs one unit solve per non-ground contact, so the cost is
    the sum across contacts). Returns None if no report is readable.
    """
    total, found = 0.0, False
    for contact_i in range(1, n_contacts + 1):
        report = os.path.join(
            f"Results_VTA_{strategy}E1C{contact_i}", "VCM_report.json"
        )
        try:
            with open(report) as fp:
                d = json.load(fp)
            total += sum(d["Timings"]["ComputeSolution"])
            found = True
        except (OSError, json.JSONDecodeError, KeyError):
            continue
    return total if found else None


def dice(mask_a, mask_b):
    """Dice overlap of two boolean masks (1.0 if both empty)."""
    a_sum = mask_a.sum()
    b_sum = mask_b.sum()
    if a_sum + b_sum == 0:
        return 1.0
    return 2.0 * np.logical_and(mask_a, mask_b).sum() / (a_sum + b_sum)


def align_to_reference(strat_keys, ref_row_of_key):
    """Indices aligning a strategy's points to the reference's points.

    Returns (strat_idx, ref_idx) so that strat_keys[strat_idx] ==
    ref_keys[ref_idx]. Meshes keep slightly different boundary lattice
    points, so Dice is computed only on points present in both.
    """
    strat_idx, ref_idx = [], []
    for j, key in enumerate(strat_keys):
        i = ref_row_of_key.get(key)
        if i is not None:
            strat_idx.append(j)
            ref_idx.append(i)
    return np.asarray(strat_idx, dtype=int), np.asarray(ref_idx, dtype=int)


def main():
    """Evaluate all strategies and write summary CSVs."""
    threshold_vpm = THRESHOLD_VPM
    spacing = load_spacing()
    voxel_volume = spacing**3  # mm^3 per lattice point
    stim_protocols = load_stim_protocols()
    n_protocols, n_contacts = stim_protocols.shape
    print(
        f"Protocols: {n_protocols}, contacts: {n_contacts}, "
        f"threshold: {threshold_vpm} V/m, spacing: {spacing} mm"
    )

    # reference (best) first
    ref_fields = load_unit_fields(REFERENCE, n_contacts)
    if ref_fields is None:
        raise FileNotFoundError(f"Reference strategy '{REFERENCE}' results missing.")
    ref_keys = ref_fields[0]
    ref_row_of_key = {key: i for i, key in enumerate(ref_keys)}
    ref_masks = compute_vta_masks(ref_fields, stim_protocols, threshold_vpm)
    ref_volumes = ref_masks.sum(axis=1) * voxel_volume

    summary_rows = []
    per_protocol = {"protocol": np.arange(n_protocols)}

    roman_idx = 0
    for strategy in STRATEGIES:
        print(f"Strategy: {strategy}")
        fields = (
            ref_fields
            if strategy == REFERENCE
            else load_unit_fields(strategy, n_contacts)
        )
        if fields is None:
            continue
        masks = (
            ref_masks
            if strategy == REFERENCE
            else compute_vta_masks(fields, stim_protocols, threshold_vpm)
        )
        # VTA volume uses the strategy's own lattice points
        volumes = masks.sum(axis=1) * voxel_volume
        per_protocol[strategy] = volumes

        # Dice needs point-aligned masks: meshes keep slightly different
        # boundary points, so align on coordinates with the reference.
        strat_idx, ref_idx = align_to_reference(fields[0], ref_row_of_key)
        n_common = strat_idx.shape[0]

        # per-protocol metrics vs. reference
        rel_err = np.where(
            ref_volumes > 0,
            np.abs(volumes - ref_volumes) / ref_volumes,
            np.nan,
        )
        dice_scores = np.array(
            [
                dice(masks[p][strat_idx], ref_masks[p][ref_idx])
                for p in range(n_protocols)
            ]
        )
        dof, elements = read_dof(strategy)
        time_s = read_time(strategy, n_contacts)
        # Roman numerals label the refinement strategies (I, II, ...) as in
        # the PAM/VTA figures; the reference run is labelled "Best".
        if strategy == REFERENCE:
            roman = "Best"
        else:
            roman_idx += 1
            roman = rf"\rom{{{roman_idx}}}"

        summary_rows.append(
            {
                "roman": roman,
                "strategy": strategy,
                "DOF": dof,
                "Elements": elements,
                "time": time_s,
                "n_common_points": int(n_common),
                "mean_VTA_volume_mm3": float(np.mean(volumes)),
                "median_VTA_volume_mm3": float(np.median(volumes)),
                "mean_rel_volume_error": float(np.nanmean(rel_err)),
                "max_rel_volume_error": float(np.nanmax(rel_err)),
                "mean_dice_vs_best": float(np.mean(dice_scores)),
                "min_dice_vs_best": float(np.min(dice_scores)),
            }
        )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv("vta_results_summary.csv", index=False)
    pd.DataFrame(per_protocol).to_csv("vta_volumes_per_protocol.csv", index=False)
    print("\nWrote vta_results_summary.csv and vta_volumes_per_protocol.csv")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
