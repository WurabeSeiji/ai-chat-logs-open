#!/bin/bash
cd "$(dirname "$0")"; export LC_ALL=en_US.UTF-8 MPLBACKEND=Agg; mkdir -p data figures results; python3 program/run_mp.py | tee results/run.log
