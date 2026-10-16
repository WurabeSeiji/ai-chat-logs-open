#!/bin/bash
set -e
python3 stage2_harmonics_translation_qualification.py
python3 stage2_posthoc_boost_conditioning.py
