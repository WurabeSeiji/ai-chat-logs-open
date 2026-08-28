# 論文v1_全プログラム修正版_20260828

論文 v1（DOI 10.5281/zenodo.22112009）／note 記事 n07c3e4c97e3a が引用する 15 パッケージの**全プログラム**に、振幅問題（親の正規化・初期 seed と正規化・位相のみ生成子）と回転問題（Cayley → 線形回転 exp(ANGLE·K)）の修正を機械的に適用し、全実験を再実行・図化し直して原本（旧エンジン）との差異を分析した。修正方針の決定はまだ行っていない。

- 全体分析：`全プログラム修正版_差異分析_20260828.md`
- 各パッケージ：`fixed/<pkg>/差異分析_修正版vs原本_20260828.md`
- 修正の差分：`results/fix_patches.diff`、検証：`results/verify_fixes_all.json`（14 エンジン全て OK）
- 比較図：`results/figures/cmp1..cmp6`

要点：修正版では seedless の潜伏相と指数成長が消え（onset step 1、Floquet 最大乗数 1.0018）、等分配・3+3+2+2・N14 quasi-closure も出ない。保存則（H_total = c²、ZᵀZ）と rank N−1 は保たれる。最大の論点は「親が修正後の力学の固定点でない」こと。

## 再現

```bash
cd 論文v1_全プログラム修正版_20260828 && bash run_all.sh   # 約 3 分
```
`../論文v1_全再現テスト_20260828/original/`（zip 展開物）が必要。python3 + numpy/pandas/matplotlib/scipy/Pillow、c++。
