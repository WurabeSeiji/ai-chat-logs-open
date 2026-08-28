#!/bin/bash
# 全パッケージ再実行。原本プログラムは results/path_patches.diff の 4 箇所以外無変更。ログは results/logs/<pkg>.log
cd "$(dirname "$0")"; export LC_ALL=en_US.UTF-8 MPLBACKEND=Agg
run_job(){ d="${1%%|*}"; c="${1#*|}"; n="$(echo "$d" | tr '/' '_')"; s=$(date +%s); ( cd "rerun/$d" && eval "$c" ) > "results/logs/$n.log" 2>&1; echo "$? $(( $(date +%s)-s ))s $d" >> results/logs/_status.txt; }
export -f run_job
: > results/logs/_status.txt
# 物理走行（独立）を並列 5、その後に依存解析を直列
grep -v -E "^(N3_N16_partial|N3_N16_nontrivial|N14_N16_complete|N5_dynamics)" results/jobs.txt | xargs -P 5 -I{} bash -c 'run_job "$@"' _ {}
grep -E "^N5_dynamics" results/jobs.txt | while read -r j; do run_job "$j"; done
( cd rerun/N3_N16_partial_zero_closure_analysis_20260826 && python3 analyze_partial_zero_closures_N3_N16.py ) > results/logs/N3_N16_partial_zero_closure_analysis_20260826.log 2>&1; echo "$? partial" >> results/logs/_status.txt
grep -E "^(N3_N16_nontrivial|N14_N16_complete)" results/jobs.txt | while read -r j; do run_job "$j"; done
echo DONE >> results/logs/_status.txt
