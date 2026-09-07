# N=100 make_parent 型親の生成（2026-09-07）

目的: 最新実験（make_parent 型全Nスイープ）と同型の初期データを N=100 で作る（極限挙動調査の準備）。

手順（規約: コピー→対照→最小変更）:
1. `run_n_scaling_lowrank_v1.py` = 7月正本エンジンの bit 同一コピー（SHA ba0fc19b…）
2. `make_static_parents_N3_N40_v1.py` = 正本生成器の bit 同一コピー（SHA c3230103…）
   → 無変更で N=3..40 を `parents/` へ再生成し、既存親（parents_actual_N3_N40）と **38/38 bit 一致**（対照合格）
3. `make_static_parent_N100_v1.py` = 上記との diff **2行のみ**（出力先 parents_N100/、`for N in [100]`）
   → `parents_N100/parent_static_N00100_makeparent_20260905.npz`（rng式・引数は正本同一）

監査（audit_parent_N100_v1.py → parents_N100/audit_parent_N100_v1.json）:
M=4950、親残差 9.80e-13、‖Z‖²=1、|Σz²|=3.6e-15、Z4位相残差 8.1e-13 rad、
4位置占有 [1228,1248,1231,1243]、振幅4950種・CV 2.17%、Wレイリー σ=98.867・残差 9.5e-13、
P_even=P_odd=0.500000000（N=3..40 と同一文法）。

再現: python3 make_static_parents_N3_N40_v1.py（対照）→ python3 make_static_parent_N100_v1.py → python3 audit_parent_N100_v1.py
