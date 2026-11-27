#!/bin/sh
# 英語版図再構成の一括実行（2026-09-15）
cd "$(dirname "$0")" || exit 1
python3 regenerate_inflation_figs_from_canonical_20260915.py
