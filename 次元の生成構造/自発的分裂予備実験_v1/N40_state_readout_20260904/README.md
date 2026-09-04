# N=40 自発的分裂走行の状態読出し（2026-09-04）

## 目的

2026-07-22 の大N自発的分裂実験（note「なぜ方向は三つで止まり、波はひとりでに増えるのか」の
`largeN_splitting_result_v1`、N=40・δ=1e-15・seed=0）は f(τ) のスカラー曲線しか保存していない。
本実験は同走行を厳密再現して step0 と最終 step（τ=3511）の全複素状態（M=780 本）を保存し、
ChatGPT_denominator_controls_N3_N33 系列と同じ複素平面読出しにかけて、
インフレーション後の準安定状態が位相凝縮構造かどうかを確認する。

## 手順（メモリ規約: コピー→対照実験→最小変更）

1. `run_spontaneous_splitting_largeN_v1.py` と依存エンジン `run_n_scaling_lowrank_v1.py` を
   `自発的分裂予備実験_v1/` から bit 同一コピー（diff で確認済み）。
2. 対照実験: コピー版を無変更で走行。元走行の引数は summary から
   `40 1e-15 --after=1500 --tol=1e-12` と同定
   （tol 既定値 1e-8 では parent_residual・crossing_tau が不一致、1e-12 で全項目 bit 一致）。
   fcurve CSV は正本 `../largeN_splitting_result_v1/fcurve_N00040_delta1e-15_seed0.csv` と全行一致。
3. 状態保存版 `run_spontaneous_splitting_largeN_v1_savestate.py`: コピーに対し
   (a) step0 状態の退避 1 行、(b) 走行後の `states_*.npz` 保存 4 行のみ追記（力学は無変更）。
   再走行後も fcurve は正本と全行一致（FCURVE_STILL_IDENTICAL）。
4. `plot_complex_plane_N40_v1.py`（読み出しのみ）で図 3 枚を生成。

## 実行

```bash
./run_all.sh
```

## 出力

- `largeN_splitting_result_v1/states_N00040_delta1e-15_seed0.npz` — Z0（step0）、Zfinal（τ=3511）、p、q
- `fig_N40_complex_plane_step0.png` — step0: 親（円偏波固有モード）＋δ種。4 本の束（対蹠 2 直径）
- `fig_N40_complex_plane_final.png` — τ=3511: 780 本がほぼ等振幅のリング状に展開
- `fig_N40_complex_plane_final_cluster_zoom.png` — 粗解像度 1/100 の最大クラスター4つの拡大。
  各 4〜5 本・内部割れ幅 ~1e-4（|z|max=3.76e-2 の ~0.5%）。12桁丸めの厳密重複は 0

## 環境

`../../ai-chat-logs-open/.venv/bin/python3`（Python 3.9.6、numpy 2.0.2、macOS arm64 Accelerate）。
ChatGPT_denominator_controls 系列と同一環境。
