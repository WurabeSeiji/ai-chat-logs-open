#!/bin/bash
cd "$(dirname "$0")"; export LC_ALL=en_US.UTF-8
j=$(grep -E "^N14_N16_complete" results/jobs.txt); d="${j%%|*}"; c="${j#*|}"; s=$(date +%s)
( cd "rerun/$d" && eval "$c" ) > results/logs/N14_N16_complete_nontrivial_zero_closure_search_20260826.log 2>&1; echo "$? $(( $(date +%s)-s ))s $d (retry with compat header)" >> results/logs/_status.txt; echo DONE2 >> results/logs/_status.txt
