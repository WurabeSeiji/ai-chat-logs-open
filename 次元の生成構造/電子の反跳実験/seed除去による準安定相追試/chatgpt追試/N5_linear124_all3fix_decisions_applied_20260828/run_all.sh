#!/bin/bash
set -e; cd "$(dirname "$0")"; export LC_ALL=en_US.UTF-8 MPLBACKEND=Agg; mkdir -p data figures results
python3 program/run_amplitude_only_fix.py > results/run.log 2> results/run_progress.log
python3 -c "
import importlib.util; spec=importlib.util.spec_from_file_location('e','program/original_engine.py'); e=importlib.util.module_from_spec(spec); spec.loader.exec_module(e); e.progress=lambda m:None
print(e.validate_against_dense(5,0,steps=300))" | tee results/validate.log
python3 program/compare_with_latest.py | tee results/compare_with_latest.log
python3 program/plot_readout_comparison.py
