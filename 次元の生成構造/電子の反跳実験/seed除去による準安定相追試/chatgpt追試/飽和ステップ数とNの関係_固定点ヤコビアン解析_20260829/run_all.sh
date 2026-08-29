#!/bin/bash
cd "$(dirname "$0")"; export LC_ALL=en_US.UTF-8 MPLBACKEND=Agg; mkdir -p data figures results
python3 program/analyze_saturation_vs_N.py 2>&1 | grep -v Warning > results/run.log
python3 program/angle_dependence_of_lambda.py 2>&1 | grep -v Warning > results/angle_dependence.log
python3 program/decompose_flow_vs_discretization.py 2>&1 | grep -v Warning > results/decompose.log
python3 program/diagnose_linear_vs_run.py 6 3000 2>&1 | grep -v Warning > results/diagnose_N6.log
( python3 program/flow_linearization_spectrum.py 4 20 0 5 2>&1 | grep -v Warning > results/flow_spectrum_N4_20_s0_5.log & python3 program/flow_linearization_spectrum.py 9 14 5 25 2>&1 | grep -v Warning > results/flow_spectrum_N9_14_s5_25.log & wait )
python3 program/hessian_inertia_energy_casimir.py 4 20 5 2>&1 | grep -v Warning > results/hessian_inertia_N4_20_s5.log
python3 program/symmetric_state_and_sigma_spread.py 2>&1 | grep -v Warning > results/symmetric_and_spread.log
