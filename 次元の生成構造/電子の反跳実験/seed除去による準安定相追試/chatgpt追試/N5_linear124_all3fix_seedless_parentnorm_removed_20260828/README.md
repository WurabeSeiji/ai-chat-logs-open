# N5_linear124_all3fix_seedless_parentnorm_removed_20260828

直前の `N5_linear124_all3fix_seedless_20260828` をコピーし、`make_parent()` 内の `v = v / np.linalg.norm(v)` だけを削除した追試。

維持した条件:
- N=5, M=10, L=124, 5000 step, SEED=0
- 外部シードなし: `Z0 = v.copy()`
- 初期 `Z /= norm(Z)` なし
- 線形指数回転
- 振幅保持 K: `Im(conj(Z_i) * Z_j)`
- `make_parent` の位相反復・固有対選択・収束判定等は変更なし

`program/` に実行コード、`data/` に全時系列CSVと全状態NPZ、`figures/` に図、`analysis_report_ja.md` に分析を収録。
