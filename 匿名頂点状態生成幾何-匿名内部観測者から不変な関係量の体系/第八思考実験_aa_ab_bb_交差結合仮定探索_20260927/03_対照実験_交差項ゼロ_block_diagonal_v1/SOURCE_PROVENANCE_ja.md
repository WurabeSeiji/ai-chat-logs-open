# 出典

本対照実験で比較に使用した保存済みデータは以下。

- Paper 7 self n1: `第七思考実験_自己相互作用_aa_bb_20260927/01_独立自己相互作用数値実験_v1/results/self_n1/raw_first_macro_microsteps.csv`
- Paper 6 n1 attractive: `第五思考実験_乗法状態内部展開と近似次数検証_20260924/13_charged8_integer_multiple_5case_strict_numeric_20260925/02_CASES/n1_attractive/raw/raw_first_macro_microsteps.csv`

実験コードは Paper 6 / 7 の `transition(z)` と同じ式を参照し、41×41 の状態依存行列 `S41(z)` を明示構成した。3チャネル版はその3ブロックを123×123の block diagonal 行列へ配置し、非対角ブロックは厳密に0とした。
