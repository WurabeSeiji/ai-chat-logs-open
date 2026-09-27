#!/bin/zsh
# 実行パターンを 1 つ実行する。 使い方: ./run_pattern.sh patterns/<パターン>.json
# 表示は画面と logs/<パターン名>.log の両方に出る。裏で実行した場合は、このログを読めば進み具合が分かる。
# h5py は専用フォルダ ~/paper6_pydeps_py39 から読む（既存の Python 環境は変えない）。
set -e
cd "$(dirname "$0")"
export PYTHONPATH="$HOME/paper6_pydeps_py39"
mkdir -p logs
name="$(basename "$1" .json)"
python3 -u wrap_run_pattern.py "$1" 2>&1 | tee "logs/$name.log"
