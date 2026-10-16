#!/bin/sh
# 全段を順に走らせる（結果は results/ に JSON で出力）
set -e
cd "$(dirname "$0")"
for s in S0_stage_gauge.py S1_generator_closure.py S2_translation_clock.py S3_stretch_triangles.py S4_self_consistency.py S5_two_body.py; do echo "== $s"; python3 "$s"; done
