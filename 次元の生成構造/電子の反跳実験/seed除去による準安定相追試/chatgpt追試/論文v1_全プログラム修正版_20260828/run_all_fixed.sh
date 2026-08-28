#!/bin/bash
# 修正版（FIX1-4）全パッケージ再実行。ログ results/logs/<job>.log、状態 results/logs/_status.txt（xargs -I は 255 byte 制限があるので bash ループで並列化）
cd "$(dirname "$0")"; export LC_ALL=en_US.UTF-8 MPLBACKEND=Agg
run_job(){ d="${1%%|*}"; c="${1#*|}"; n="$(echo "$d" | tr '/' '_')"; s=$(date +%s); if [ "$d" = "N16_analyze" ]; then ( cd fixed && eval "$c" ) > "results/logs/$n.log" 2>&1; else ( cd "fixed/$d" && eval "$c" ) > "results/logs/$n.log" 2>&1; fi; echo "$? $(( $(date +%s)-s ))s $d" >> results/logs/_status.txt; }
: > results/logs/_status.txt; date +%s > results/run_start_epoch.txt
# 物理走行（独立）を並列 5
while read -r j; do run_job "$j" & while [ "$(jobs -rp | wc -l)" -ge 5 ]; do sleep 1; done; done < <(grep -v -E "^(N16_analyze|N3_N16_partial|N3_N16_nontrivial|N14_N16_complete|N5_dynamics)" results/jobs.txt); wait
# 依存する解析は直列
while read -r j; do run_job "$j"; done < <(grep -E "^(N5_dynamics|N16_analyze|N3_N16_partial|N3_N16_nontrivial|N14_N16_complete)" results/jobs.txt)
echo DONE >> results/logs/_status.txt
