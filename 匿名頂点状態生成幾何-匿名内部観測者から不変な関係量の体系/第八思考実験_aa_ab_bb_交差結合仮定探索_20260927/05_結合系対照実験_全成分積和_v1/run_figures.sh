#!/bin/zsh
# 結合系の観測 row データから図を作る。 使い方: ./run_figures.sh <pattern_id>
set -e
cd "$(dirname "$0")"
export PYTHONPATH="$HOME/paper6_pydeps_py39"
export MPLBACKEND=Agg
mkdir -p logs
{
  python3 -u fig00_copy_figure_programs.py
  python3 -u fig20_analytic_reference.py "$1"
  python3 -u fig30_per_relation_figures.py "$1"
  python3 -u fig40_orbit_condition_comparison.py "$1"
  python3 -u fig50_log_radial_orbits.py "$1"
  python3 -u fig60_analytic_reference_figures.py "$1"
} 2>&1 | tee "logs/figures_$1.log"
