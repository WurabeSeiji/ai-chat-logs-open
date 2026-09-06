#!/bin/zsh
# step0/完了時 複素平面 対照テスト10ケース（読み出しのみ・新規走行なし）
cd "$(dirname "$0")"
echo "=== 対照テスト（Aのcommitted図を再現・幾何一致確認） ==="
python3 control/control_reproduce_A_figs_v1.py
echo "=== 10ケース図化（step0｜完了時） ==="
python3 plot_step0_final_10cases_v1.py
echo "=== 全38ケース図化（N=3..40） ==="
python3 plot_step0_final_allN_v1.py
