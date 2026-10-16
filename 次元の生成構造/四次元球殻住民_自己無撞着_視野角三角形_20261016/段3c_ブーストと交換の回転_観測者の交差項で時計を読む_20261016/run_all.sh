#!/bin/bash
set -e
python3 stage3c_boost_exchange_beat.py
python3 stage3c_posthoc_per_degree_beat.py
