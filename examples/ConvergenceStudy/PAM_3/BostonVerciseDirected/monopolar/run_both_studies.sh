#!/bin/bash
# Run both convergence studies sequentially (floating then nonfloating)
# Usage: nohup bash run_both_studies.sh > run_both_studies.log 2>&1 &

set -e
cd "$(dirname "$0")"
source /home/julius/Documents/github/OSS-DBSv2/venv/bin/activate

echo "=== Starting floating convergence study ==="
echo "Start time: $(date)"
python run_convergence_study.py
echo "Floating study finished: $(date)"

echo ""
echo "=== Starting nonfloating convergence study ==="
echo "Start time: $(date)"
python run_convergence_study_nonfloating.py
echo "Nonfloating study finished: $(date)"

echo ""
echo "=== Both studies complete ==="
echo "End time: $(date)"
