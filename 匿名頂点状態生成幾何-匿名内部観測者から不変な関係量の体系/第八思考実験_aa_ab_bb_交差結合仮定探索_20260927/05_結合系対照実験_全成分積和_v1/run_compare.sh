#!/bin/zsh
# 既存データとの比較を 1 つ実行する。 使い方: ./run_compare.sh comparisons/<比較の定義>.json
set -e
cd "$(dirname "$0")"
export PYTHONPATH="$HOME/paper6_pydeps_py39"
mkdir -p logs
name="$(basename "$1" .json)"
python3 -u compare_with_existing.py "$1" 2>&1 | tee "logs/$name.log"
