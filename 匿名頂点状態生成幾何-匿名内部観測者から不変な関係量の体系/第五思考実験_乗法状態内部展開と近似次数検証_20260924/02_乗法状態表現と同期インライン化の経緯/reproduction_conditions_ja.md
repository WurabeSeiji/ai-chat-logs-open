# 再現条件

このフォルダには複数の歴史的案を意図的に保存している。

- `audit_five_state_multiplicative.py`: exp/log 往復を使う**監査用**シェル。
- `audit_multiplicative_rk4_shell.py`: RK4の加算部だけを指数積へ持ち上げる監査。
- `audit_five_state_full_inline.py`: 同期インライン写像と狭い固定線形写像テスト。
- `audit_formula_equivalence_only.py`: 100000乱数点で逐次式とインライン式を照合。seed=20260924。
- `audit_pure_exponential_state_engine.py`: 永続5状態を乗法表現にした長時間監査。
- `validate_latest_multiplicative_state.py`: その後の状態定義確認用。

外部基準CSVは `../00_reference_from_paper4/experiment02_sn_2p5pn_rr_C1_N12.csv`。
