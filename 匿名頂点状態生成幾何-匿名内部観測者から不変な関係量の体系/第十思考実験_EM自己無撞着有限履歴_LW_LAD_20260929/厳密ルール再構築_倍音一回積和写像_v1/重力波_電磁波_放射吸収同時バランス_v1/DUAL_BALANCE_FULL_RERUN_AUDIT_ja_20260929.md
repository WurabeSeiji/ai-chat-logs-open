# GW/EM 同時バランス 完全再実行監査 — 2026-09-29

`RUN_DUAL_BALANCE_ALL_v1.sh` を保存済みコード・固定 interaction array から完全再実行した。

## 結果

- 全4ケースの探索を最初から再実行。
- 4個の `dual_q*_v2.json` は、再実行前の保存版と **SHA-256 が完全一致**。
- したがって探索履歴を含む JSON 全体がバイト単位で一致。
- `verify_dual_balance_repro_v1.py` の独立直接再実行でも `all_exact_reproduction = true`。
- flux と delta の保存値との差は全ケースで `0.0`。

| Q/M | delta_GW | delta_EM | JSON byte exact |
|---:|---:|---:|:---:|
| 0.30 | 3.051828348031962e-09 | -9.156081764941546e-11 | yes |
| 0.60 | 7.192485214985341e-10 | -5.546214863761080e-11 | yes |
| 0.90 | 3.293698849411545e-13 | -4.487463021160879e-14 | yes |
| 0.99 | 1.445122605426921e-12 | -9.839518632864276e-13 | yes |

## 再現条件

- `strict_harmonic_rn_tensor_v1.py`
- `strict_rn_interaction_v1.npz`
- `search_dual_balance_strict_v2.py`
- `verify_dual_balance_repro_v1.py`
- `RUN_DUAL_BALANCE_ALL_v1.sh`
- `requirements_dual_balance_v1.txt`

再実行で transition / interaction array は変更していない。
