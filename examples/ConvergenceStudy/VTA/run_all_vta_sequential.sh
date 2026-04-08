#!/usr/bin/env bash
# Run all 6 VTA convergence studies sequentially (no parallelism to avoid OOM).
# Each run executes the full run_convergence_study.py (all strategies).
set -u

cd "$(dirname "$0")"
VTA_ROOT="$(pwd)"
MASTER_LOG="${VTA_ROOT}/run_all_vta_sequential.log"

# Use venv python
VENV_PY="${VTA_ROOT}/../../../venv/bin/python"

# Limit threads to avoid ngsolve/openblas oversubscription
export OPENBLAS_NUM_THREADS=8
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8

DIRS=(
    "BostonScientificVerciseDirected/adjacent_contacts"
    "BostonScientificVerciseDirected/case_grounding"
    "BostonScientificVerciseDirected/distant_contacts"
    "Medtronic3389/adjacent_contacts"
    "Medtronic3389/case_grounding"
    "Medtronic3389/distant_contacts"
)

: > "$MASTER_LOG"
echo "=========================================" | tee -a "$MASTER_LOG"
echo "VTA sequential run started: $(date)"       | tee -a "$MASTER_LOG"
echo "Python: $VENV_PY"                          | tee -a "$MASTER_LOG"
echo "Threads: OPENBLAS/OMP/MKL = 8"             | tee -a "$MASTER_LOG"
echo "=========================================" | tee -a "$MASTER_LOG"

for d in "${DIRS[@]}"; do
    echo ""                                      | tee -a "$MASTER_LOG"
    echo "=========================================" | tee -a "$MASTER_LOG"
    echo "[$(date)] START: $d"                   | tee -a "$MASTER_LOG"
    echo "=========================================" | tee -a "$MASTER_LOG"

    cd "${VTA_ROOT}/${d}"

    # Clean up stale Results_VTA_best from the prior OOM-killed run
    # so the fresh run starts from a clean state.
    if [ -d Results_VTA_best ]; then
        rm -rf Results_VTA_best
        echo "Removed stale Results_VTA_best" | tee -a "$MASTER_LOG"
    fi

    # Capture run time and exit code
    start_ts=$(date +%s)
    if "$VENV_PY" run_convergence_study.py >> "$MASTER_LOG" 2>&1; then
        status="OK"
    else
        status="FAIL(exit=$?)"
    fi
    end_ts=$(date +%s)
    elapsed=$(( end_ts - start_ts ))

    echo "[$(date)] END:   $d  status=$status  elapsed=${elapsed}s" | tee -a "$MASTER_LOG"

    cd "$VTA_ROOT"
done

echo ""                                          | tee -a "$MASTER_LOG"
echo "=========================================" | tee -a "$MASTER_LOG"
echo "VTA sequential run finished: $(date)"      | tee -a "$MASTER_LOG"
echo "=========================================" | tee -a "$MASTER_LOG"
