#!/bin/zsh
# 第七思考実験の再実行。第六と補遺のプログラムを読み込み、条件と保存先だけを変える。
# h5py と pandas は専用フォルダ ~/paper6_pydeps_py39 から読む（既存の Python 環境は変えない）。
set -e
cd "$(dirname "$0")"
export PYTHONPATH="$HOME/paper6_pydeps_py39"
export MPLBACKEND=Agg
mkdir -p logs
python3 step00_copy_programs.py            2>&1 | tee logs/step00_copy_programs.log
python3 wrap10_run_generator.py            2>&1 | tee logs/wrap10_run_generator.log
python3 wrap20_analytic_reference.py       2>&1 | tee logs/wrap20_analytic_reference.log
python3 wrap30_postprocess_figures.py      2>&1 | tee logs/wrap30_postprocess_figures.log
python3 wrap40_orbit_condition_comparison.py 2>&1 | tee logs/wrap40_orbit_condition_comparison.log
python3 wrap50_log_radial_orbits.py        2>&1 | tee logs/wrap50_log_radial_orbits.log
python3 wrap60_analytic_reference_figures.py 2>&1 | tee logs/wrap60_analytic_reference_figures.log
python3 compare_with_previous.py           2>&1 | tee logs/compare_with_previous.log
python3 inspect_Q_difference.py           2>&1 | tee logs/inspect_Q_difference.log
python3 step90_write_sha256sums.py         2>&1 | tee logs/step90_write_sha256sums.log
