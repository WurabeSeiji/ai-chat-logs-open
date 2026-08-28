# N5_linear124_all3fix_decisions_applied_20260828

最新の修正プログラム（`../N5_linear124_all3fix_seedless_parentnorm_removed_40000_20260828`）に、確定した判断 1〜6（A5 未使用関数の正規化除去／A6(b) σ₁ を実際の生成子から読む＋実測位相進み／S1 seed 関数を呼ばない／R3(ii) validate を修正後の力学に書換え）だけを適用し、力学が不変であることを確認した。結果は `比較結果_判断1-6適用版vs最新版_20260828.md`。

- 力学（全観測量・全状態）は無変更版と 0.0 差で一致。変わるのは σ₁ の読出し（treatment 2.4〜3.8 → 0.36〜0.43）のみ。
- baseline で実測位相進み = (2π/L)·σ₁ が 10⁻¹⁶ で成立（読出し器の検算）。

```bash
bash run_all.sh   # 約 15 秒。対照 ../N5_linear124_all3fix_control_rerun_40000_20260828/data が必要
```
