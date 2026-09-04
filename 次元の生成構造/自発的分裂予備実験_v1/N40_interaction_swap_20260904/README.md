# N=40 相互作用スワップ実験（2026-09-04）

## 目的

N=40 自発的分裂走行（δ=1e-15, seed=0, tol=1e-12, after=1500）について、
**step0（親＋種）を bit 同一に保ったまま、step ループの相互作用だけを入れ替え**、
インフレーション的増幅→分裂が相互作用のどの性質に依るかを単離する。

- 旧相互作用: `set_theta(np.angle(Z))`（振幅を剥ぎ取り位相のみで K_ef=sin(θ_f−θ_e) を構成）
  ＋ `K̃=K/σ_max`（毎ステップのスペクトル正規化）＋ 実直交 Cayley 回転
- 新相互作用: `one_step`＝振幅込み Hermitian H_ef=conj(z_e)z_f のユニタリ exp(−i(2π/den)H)、
  den=40（Δτ=2π/N）、正規化なし

## 手順（コピー→対照→最小差し替え）

1. `run_spontaneous_splitting_largeN_v1_savestate.py`（N40_state_readout_20260904 で
   fcurve の正本 bit 一致を確認済みの状態保存版）、`run_n_scaling_lowrank_v1.py`、
   相互作用の出典 `run_and_plot_N3_N33_legacyparent_20260903.py` を bit 同一コピー。
2. 辺順序の同一性を事前検証: `build_edges(40) == np.triu_indices(40,k=1)`（True）。
3. `run_spontaneous_splitting_largeN_v1_interactionswap.py`: savestate 版コピーに対し
   (a) `edges/adjacency/H_of/one_step` 4関数を出典から**逐語**追加（文字単位一致を機械検証、
   ALL VERBATIM）、(b) ループ内の旧3行を `Z = one_step(Z, A_int, den)` に置換、
   (c) `--den` フラグ（既定 den=N）と tag/summary への den 記録。親生成・種・正規化・
   p,q 読出し・f 測定・停止条件・rng 消費は一切不変。
4. 走行: `40 1e-15 --after=1500 --tol=1e-12`（den=40）。
5. 検証: states npz の Z0・p・q が N40_state_readout の正本走行と **bit 同一**（True）。

## 結果

| | 旧相互作用（正本） | 新相互作用（スワップ） |
|---|---|---|
| f の立ち上がり | 64 steps/decade の一定指数増幅 | 約100 step で一気に ~1e-3 |
| 0.05 交差 | τ=2011 | **なし**（cap=8000 まで未達） |
| 終値 | f=0.166（分裂・準安定） | f=0.0062（振動＋緩慢な漂い） |
| 終状態の複素平面 | リング状に位相分散 | **4束構造を保持**（束が僅かに滲む） |
| \|Z·Z\|（零閉塞） | ~1.4e-15 で保存 | 1.6e-3 まで成長（保存しない） |

**結論: 相互作用の入れ替えだけで結果が反転する。** 7月論文のインフレーション→分裂は
位相のみ・σ正規化つき相互作用に固有の現象で、振幅重み付き exp(−iΔH) 力学では同じ
初期状態（自己無撞着親＋零閉塞種）はほぼ安定に留まる。

## 実行

```bash
./run_all.sh
```

## 出力

- `largeN_splitting_result_v1/fcurve_N00040_delta1e-15_seed0_den40.csv` / `summary_*.json` / `states_*.npz`
- `fig_N40_fcurve_interaction_swap_compare.png` — f(τ) の新旧比較（step0 は bit 同一）
- `fig_N40_swap_complex_plane_final.png` — 新相互作用の終状態（τ=8000）
- `fig_N40_swap_complex_plane_final_cluster_zoom.png` — 角クラスター拡大（各5〜10本、割れ幅 ~2e-4）

## 環境

`.venv/bin/python3`（Python 3.9.6、numpy 2.0.2、macOS arm64 Accelerate）。全系列と同一。
