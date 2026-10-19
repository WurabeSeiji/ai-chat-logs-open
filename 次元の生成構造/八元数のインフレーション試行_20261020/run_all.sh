#!/bin/bash
set -e
cd "$(dirname "$0")"
python3 octonion_inflation_trials_20261020.py
python3 covariant_generators_on_designs_20261020.py
python3 su3_rotating_families_20261020.py
python3 analytic_kernel_types_and_stability_20261020.py
shasum -a 256 octonion_inflation_trials_20261020.py covariant_generators_on_designs_20261020.py su3_rotating_families_20261020.py analytic_kernel_types_and_stability_20261020.py results/* > SHA256SUMS.txt
