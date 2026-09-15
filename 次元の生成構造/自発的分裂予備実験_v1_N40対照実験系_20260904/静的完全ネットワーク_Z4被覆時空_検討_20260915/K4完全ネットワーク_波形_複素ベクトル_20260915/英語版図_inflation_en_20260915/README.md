# 英語版図の再構成: inflation 系再掲図（図1・図2・図4）2026-09-15

## 目的

K4 論文 `time_vertex_complete_relation_closed_system_UK_K4_even_odd_harmonics_ja.md` が
再掲する inflation 系 3 図（図1: step0 十字、図2: step500 リング、図4: 3D 軌道）は
図中テキストが日本語であり、英訳版論文用の英語版が存在しなかった（図3 Mexican-hat は
元から英語のため対象外）。既存 3 図の生成スクリプトは未保存だったため、正本データから
再構成した。

## 正本（読み取りのみ）

- `…/make_parent型初期値_自己無撞着構造_20260907/full_N3_N40_sweep/states/hm_N40_den_40_states_500.npz`
  — Z 履歴 (501 step × 780 波, complex128)
- `…/make_parent型初期値_自己無撞着構造_20260907/plot3d_makeparent_selectedN_v1.py` の
  `series(40)` — 集団回転位相 θ0 の成分 (X, Y) と H⊥/H 系列 fz（式の再実装なし）

表示は正本 `plot3d_relphase_complexplane_v2.py` と同一の共回転フレーム
w(τ) = z(τ)·exp(−iθ0(τ))。θ0 の符号は同スクリプトの床静止基準
（FLOOR_STATIC_ABS=0.2, FLOOR_SIGN_RATIO=20）をコピーして自動決定。

## 対照結果（contrast_ja/ と既存図の突合）

- 図1: 4 方向帯の角度位置・振幅・タイトル数値 H⊥/H = 1.00e-31 が既存図と一致。
- 図2: リング構造・ギャップ位置・H⊥/H = 1.01e-01 が既存図と一致。
- 図4: 底 log10(H⊥/H)≈−31 からの直立柱 → 上端リング分岐の構造が既存図と一致。
- 生フレーム（θ0 除去なし）では 3D が螺旋円筒になり既存図と不一致（v1 で確認）。
  既存図が共回転フレーム表示であることの同定根拠。

数値注記: timeseries CSV（Hperp_frac）の step0 実測は 3.550e-32 だが、正本
`series()` の fz は床値 1e-31 でフロアされた系列であり、既存図・本再構成の
表示値 1.00e-31 はこちらに一致する。

## ファイル

- `regenerate_inflation_figs_from_canonical_20260915.py` — 再構成プログラム（ja/en 両方生成）
- `contrast_ja/` — 対照用日本語版 3 図（既存図との突合用）
- `en/` — 英語版 3 図（英訳版論文用）
  - `inflation_N40_phase_step0_cross_en_20260915.png`
  - `inflation_N40_phase_step500_ring_en_20260915.png`
  - `inflation_N40_relation_wave_3D_step500_en_20260915.png`
- `run_all.sh` — 一括再実行
- `SHA256SUMS.txt`

## 再実行

```sh
./run_all.sh
```
