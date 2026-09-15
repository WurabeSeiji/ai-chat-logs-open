#!/bin/sh
# Phase 0 監査の一括再実行（2026-09-15）
cd "$(dirname "$0")" || exit 1
python3 run_phase0_audit.py
