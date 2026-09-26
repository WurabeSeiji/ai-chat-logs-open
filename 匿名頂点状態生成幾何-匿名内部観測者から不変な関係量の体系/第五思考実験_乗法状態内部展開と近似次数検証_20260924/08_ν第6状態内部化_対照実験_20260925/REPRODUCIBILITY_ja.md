# 再現手順 — ν 第6状態内部化 対照実験

## 必要ファイル

- `run_nu_internal_state_control_v1.py`
- `plot_nu_internal_state_control_v1.py`
- `reference_audit_two_paths_rk4_snapshot.py`
- `reference_audit_method_B_strict_snapshot.py`
- `reference_verify_method_B_strict_v2_snapshot.py`
- `experiment02_sn_2p5pn_rr_C1_N12.csv`

## 実行

```bash
python run_nu_internal_state_control_v1.py
python plot_nu_internal_state_control_v1.py
```

実験プログラムは次を実行する。

1. 12ケース × 4000 macro steps の旧5状態 vs 新5+1状態比較
2. 代表3ケース × 12000 macro steps の長時間比較
3. 保存済み Paper-4 C1 CSV との照合
4. 第6状態 $\mathcal N$ の microstep 恒等保持監査
5. 新遷移関数に外部 `nu` 引数が残っていないかの静的監査

図化プログラムは保存済みCSV/JSONのみを読み、力学計算を再実行しない。

## 主要出力

- `nu_internal_state_control_timeseries.csv` — 12ケースの時系列raw/対照データ
- `nu_internal_state_control_case_summary.csv` — ケース別最大差
- `nu_internal_state_control_summary.json` — 全監査結果
- `long_run_checks.json` — 3ケース×12000 step
- `saved_C1_check.json` — 保存済みC1照合
- `nu_identity_interaction_coefficients.json` — 第6行恒等係数
- `figures/*.png`, `figures/*.svg` — 図化結果

## 期待される合格条件

```text
old/new first-five-state max abs difference = 0.0
old/new first-five-state bitwise mismatch   = 0
N vs nu0 bitwise mismatch                   = 0
N microstep hold failure                    = 0
external nu argument in new transition      = absent
saved C1 error profile old == new            = true
```

## 注意

この対照実験は ν 内部化だけを検証する。既存RK4内部の加法stageなど、別の厳密ルール監査項目は変更していない。
