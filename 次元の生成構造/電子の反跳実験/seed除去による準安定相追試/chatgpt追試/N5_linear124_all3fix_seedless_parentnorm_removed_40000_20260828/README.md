# N5 linear124 all3fix seedless parentnorm removed — 40000-step sweep

直前の `N5_linear124_all3fix_seedless_parentnorm_removed_20260828` をそのまま複製し、**STEPSだけ 5000 → 40000** に変更した長時間追試。

固定条件:
- N=5, M=10, L=124, SEED=0, DELTA=0
- 外部 seed `DELTA*g` なし (`Z0=v`)
- 初期 `Z/=norm(Z)` なし
- `make_parent` 内 `v/=norm(v)` なし
- 線形指数回転 `exp((2*pi/124)K)`
- 振幅保持相互作用 `K_ij = Im(conj(Z_i) Z_j)`
- baseline phase-only K も同じ初期親から同時実行

変更点:
- `STEPS=40000`
- key-step CSV に 7500, 10000, 15000, ..., 40000 を追加

成果物:
- `program/` 実行コード
- `data/` 全時系列CSV、全状態NPZ、summary.json、key_steps.csv
- `figures/` H_perp、PR/M、振幅分散、二乗閉包残差
- `analysis_report_ja.md`
- `run_40000.log`
- `SHA256SUMS.txt`
