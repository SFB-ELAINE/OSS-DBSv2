#!/usr/bin/env bash
# Run both panels strictly sequentially (CLAUDE.md rule).
# Activate the venv before sourcing this script:
#   source venv/bin/activate
#   bash examples/Electrochemistry/BCComparison/run_all.sh
set -euo pipefail
cd "$(dirname "$0")"
python run_panel_A.py
python run_panel_B.py
python summarise.py
