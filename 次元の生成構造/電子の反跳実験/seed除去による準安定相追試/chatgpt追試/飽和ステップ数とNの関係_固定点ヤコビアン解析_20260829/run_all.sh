#!/bin/bash
cd "$(dirname "$0")"; export LC_ALL=en_US.UTF-8 MPLBACKEND=Agg; mkdir -p data figures results
python3 program/analyze_saturation_vs_N.py 2>&1 | grep -v Warning > results/run.log
python3 program/angle_dependence_of_lambda.py 2>&1 | grep -v Warning > results/angle_dependence.log
python3 program/decompose_flow_vs_discretization.py 2>&1 | grep -v Warning > results/decompose.log
python3 program/diagnose_linear_vs_run.py 6 3000 2>&1 | grep -v Warning > results/diagnose_N6.log
