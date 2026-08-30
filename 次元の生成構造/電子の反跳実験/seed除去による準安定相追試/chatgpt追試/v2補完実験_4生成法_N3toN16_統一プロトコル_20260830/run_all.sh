#!/bin/bash
# 全再実行：親生成＋走行前予測 → 埋め込み実験 → 走行 54 本（2 系統並列、約 1 時間）→ 集計 → 図 → SHA
set -e; set -o pipefail
cd "$(dirname "$0")"
python3 program/pass1_parents.py       > results/pass1.log
python3 program/pass2_embed_random.py  > results/pass2.log
( for N in 3 4 5 6 7 8 9 10 11 12 13 14 15 16; do for m in mp hm; do python3 program/run_dynamics.py ${m}_N$N; done; done > results/run_dynamics_A.log 2>&1 ) &
( for N in 3 4 5 6 7 8 9 10 11 12 13 14 15 16; do for m in ne rb; do [ -d data/${m}_N$N ] && python3 program/run_dynamics.py ${m}_N$N; done; done > results/run_dynamics_B.log 2>&1 ) &
wait
gzip -f data/*/treatment_linear124_amplitude_aware_timeseries.csv
python3 program/pass5_analysis.py      > results/pass5.log
python3 program/pass6_figures.py
find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 shasum -a 256 > SHA256SUMS.txt
echo done
