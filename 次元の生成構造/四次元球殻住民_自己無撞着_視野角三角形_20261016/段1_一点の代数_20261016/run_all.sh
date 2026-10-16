#!/bin/bash
set -e
python3 stage1_single_point_algebra.py
python3 stage1_posthoc_relative_closure.py
