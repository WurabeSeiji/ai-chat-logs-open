# N=40 対照実験系（2026-09-04）

2026-07-22 の大N自発的分裂実験のうち N=40 走行を、物理・初期データを一切変更せずに
別フォルダで忠実再現し、正本図化＋状態読出し図を生成する対照実験系。

## 構成ファイル（すべてコピー元と bit 同一を diff で確認済み）

| ファイル | コピー元 | 変更 |
|---|---|---|
| run_n_scaling_lowrank_v1.py | 自発的分裂予備実験_v1/ | なし |
| run_spontaneous_splitting_largeN_v1.py | 自発的分裂予備実験_v1/ | なし |
| make_largeN_figure_v1.py | 自発的分裂予備実験_v1/ | なし |
| run_spontaneous_splitting_largeN_v1_savestate.py | 同_v1/N40_state_readout_20260904/ | 元走行との差分は状態保存の追記2箇所のみ（力学無変更） |
| plot_complex_plane_N40_v1.py | 同_v1/N40_state_readout_20260904/ | なし（読み出しのみ） |

## 実行（ラッパー）

```bash
./run_N40_only.sh
```

引数は正本 summary から同定済みの元走行条件 `40 1e-15 --after=1500 --tol=1e-12` に固定。

1. 元プログラム無変更で走行 → **GATE1**: fcurve が正本
   `../自発的分裂予備実験_v1/largeN_splitting_result_v1/fcurve_N00040_delta1e-15_seed0.csv`
   と全 3512 行 bit 一致で合格
2. 状態保存版で走行 → **GATE2**: fcurve が再び正本と bit 一致（保存追記が物理に影響しないことの検証）
3. `make_largeN_figure_v1.py` を無変更で実行 → 結果フォルダに N=40 しかないため
   元の図化のまま N=40 単独のインフレーション図になる
4. `plot_complex_plane_N40_v1.py` で追加読出し図（step0・最大step複素図・凝縮部拡大図）

## 2026-09-04 実行結果

- GATE1 PASS / GATE2 PASS（両走行とも crossing τ=2011、正本と bit 一致）
- 出力:
  - `largeN_splitting_result_v1/dormant_growth_large_n_v1.png` — **インフレーション図（N=40、元図化のまま）**
  - `largeN_splitting_result_v1/fcurve_*.csv` / `summary_*.json` / `states_*.npz`（Z0・Zfinal τ=3511・p・q）
  - `fig_N40_complex_plane_step0.png` — step0 複素図（4束構造）
  - `fig_N40_complex_plane_final.png` — 最大 step（τ=3511）複素図（リング）
  - `fig_N40_complex_plane_final_cluster_zoom.png` — 凝縮部（角クラスター）拡大図

## 静的親データ（2026-09-04 追加・make_static_parent_N40_v1.py）

本フォルダの検証済みエンジンから make_parent / zero_closure_kernel_seed を import し、
正本走行と同一手順・同一 rng 消費順で v（親）・g（零閉塞種）・Z0（正規化初期状態）を生成、
既存データを一切上書きしない別名で保存:

- `largeN_splitting_result_v1/parent_static_N40_makeparent_20260904.npz`
  （v, g, Z0, sigma, residual, n=40, seed=0, delta=1e-15, tol=1e-12, iters=1200）

**GATE PASS（2026-09-04 実行）**: 生成 Z0 が既存
`states_N00040_delta1e-15_seed0.npz` の Z0 と bit 一致（residual=6.237953231674313e-13 も一致）。
以後、新プログラム側の N=40 初期データ差し替えにはこの静的ファイルを使う
（make_parent の動的呼び出しは行わない）。

## 環境

`.venv/bin/python3`（Python 3.9.6、numpy 2.0.2、macOS arm64 Accelerate）
