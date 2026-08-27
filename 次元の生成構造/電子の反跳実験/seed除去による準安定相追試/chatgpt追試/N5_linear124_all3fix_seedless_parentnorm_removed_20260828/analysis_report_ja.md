# N=5 make_parent 振幅正規化除去 追試分析

## 変更点
直前の seedless/all-three-fixes 版から、`make_parent()` 内の次の1行だけを削除した。

```python
v = v / np.linalg.norm(v)
```

外部シード除去、初期 Z 正規化除去、線形指数回転、振幅保持相互作用 `Im(conj(Z_i)*Z_j)` は維持した。

## 実測
- parent residual: 0.48689998519066
- 初期/保存全ノルム H_total: 0.859649122807157
- 振幅保持版 H_perp/H_total > 0.05: step 76
- 振幅保持版 max H_perp: 0.429824543026578 (step 2188)
- 振幅保持版 max H_perp fraction: 0.499999978622691
- step 5000 PR/M: 0.489197887029825
- step 5000 amplitude range: 0.181916212445596 .. 0.587362688613061
- H_total max drift: 1.734e-13

## 重要な観察
`make_parent` 内の振幅正規化を外すと、parent のスケールが1ではなくなり、H_total は約 0.859649122807 になった。また parent residual は 0.4869 である。

したがって、この正規化は単なる表示上の正規化ではなく、現在の `make_parent` の残差評価・固定点反復とスケールの関係に影響している可能性がある。今回の追試結果は保存するが、「正規化を外した方が理論的に正しい」とはこの実験だけから判定しない。
