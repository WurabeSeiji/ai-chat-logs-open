# 再現条件

- 静的AST監査: `micro_fast` の右辺で新規辞書 `n[...]` を Load していないこと。
- 乱数500状態の1 macro-step同値性。
- one-hot相の破損、外部状態hold、11相後のcycle復帰、work reset を監査。
- 12ケース×4000 macro steps。
- 代表3ケース×12000 macro steps。
- 保存済みPaper 4 C1 CSVとの照合。
- seed=20260924。
