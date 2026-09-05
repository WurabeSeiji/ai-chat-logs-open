#!/bin/zsh
# 論文A（DOI 10.5281/zenodo.21874482）対照再現の全段実行。
# プログラムは凍結マニフェストSHA照合済みの無改変コピー。出力は本フォルダのみ（正本不可侵）。
# 各段は非致命（FAILEDを記録して続行）。進捗は timestamp つきで標準出力へ。
cd "$(dirname "$0")"
log() { echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
run() { log "--- $1 start"; shift; "$@"; local rc=$?; [ $rc -ne 0 ] && log "!!! FAILED (exit $rc)"; log "--- done (exit $rc)"; }

log "===== 対照再現開始（論文A・全段） ====="

log "=== STAGE 1: 三系 N掃引 1..20 ==="
run "neutral"  python3 -u run_nsweep_three_series_v2.py neutral  1 20
run "electron" python3 -u run_nsweep_three_series_v2.py electron 1 20
run "vacuum"   python3 -u run_nsweep_three_series_v2.py vacuum   1 20

log "=== STAGE 2: mixed 1..16 ==="
run "mixed" python3 -u run_nsweep_three_series_v2.py mixed 1 16

log "=== STAGE 3: シード型別δ一括掃引 T42000（登録グリッド全走行） ==="
run "missing_seed_sweeps" python3 -u run_missing_seed_sweeps_T42000_v1.py

log "=== STAGE 4: phase-balanced 単点 δ=0.04357 ==="
run "pb-control"  python3 -u run_phase_balanced_mixed_v1.py run --arm control --replicate ctl1 --allow-science-run
run "pb-balanced" python3 -u run_phase_balanced_mixed_v1.py run --arm phase-balanced --replicate r1 --allow-science-run

log "=== STAGE 4b: δ格子拡張 {0.01, 0.03162277660168379, 0.1} ==="
run "grid-d0.01"    zsh run_grid_chain.sh 0.01 d001
run "grid-d0.0316"  zsh run_grid_chain.sh 0.03162277660168379 d00316
run "grid-d0.1"     zsh run_grid_chain.sh 0.1 d01

log "=== STAGE 5: electron_affine16 16to64 診断（全腕） ==="
for A in A_standard B_project16 C_allowed_e1e-05 D_forbidden_e1e-05 D_forbidden_e1e-10 D_forbidden_e1e-15; do
  run "arm-$A" python3 -u run_electron_affine16_instability_diagnostic_v1.py --arm "$A"
done

log "=== STAGE 6: 第3段階 root probe（事前登録 pilot 行列 6走行） ==="
run "stage3" zsh run_stage3_chain.sh \
  mixed 0.03162277660168379 0.0316228 \
  mixed 0.04357 0.04357 \
  neutral 0.03162277660168379 0.0316228 \
  electron 0.03162277660168379 0.0316228 \
  fermion_family 0.03162277660168379 0.0316228 \
  boson_family 0.03162277660168379 0.0316228

log "=== STAGE 7: 約数類定理レジスタ位数追随検定 ==="
run "divisor" python3 -u run_divisor_class_register_order_v1.py

log "=== STAGE 8: 長時間走行オーケストレータ（T=300000 含む・最長） ==="
run "stage4-longtime" python3 -u run_stage4_longtime_orchestrator_v1.py --run-all

log "=== STAGE 9: 主張数値の集計 ==="
run "aggregate" python3 -u aggregate_paperA_claims_v1.py

log "=== STAGE 10: 正本との全数照合 ==="
run "compare" python3 -u compare_repro_vs_canonical_v1.py

log "===== 全段終了 ====="
