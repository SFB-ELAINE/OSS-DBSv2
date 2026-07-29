#!/usr/bin/env bash
# Sequential driver for the StimSets truncation sweep.
#
# Strictly one ratio at a time (see CLAUDE.md): each invocation of
# run_truncation_study.py finishes completely before the next starts.
# Ratio 20 runs with --stage pam because its eight unit solutions were
# already produced by the reference-validity check.
#
# The untruncated reference is not re-run: it is reused from
# ../Results_PAM_hp_material_refinement, verified equivalent to the
# current code to within solver noise (3e-11 relative).
#
# Usage: setsid nohup bash run_sweep.sh > sweep.log 2>&1 < /dev/null &

set -u
cd "$(dirname "$0")"

PY="uv run python"
rm -f sweep.done

run_stage() {
  local ratio="$1" stage="$2"
  echo "=== $(date -Is) starting ratio ${ratio} (stage ${stage}) ==="
  # shellcheck disable=SC2086
  if ! $PY run_truncation_study.py "${ratio}" --stage "${stage}"; then
    echo "=== $(date -Is) ratio ${ratio} FAILED ==="
    echo "FAILED ratio=${ratio} stage=${stage}" > sweep.done
    exit 1
  fi
  echo "=== $(date -Is) finished ratio ${ratio} ==="
}

run_stage 5 all
run_stage 10 all
run_stage 20 pam
run_stage 30 all

echo "=== $(date -Is) sweep complete ==="
echo "OK" > sweep.done
