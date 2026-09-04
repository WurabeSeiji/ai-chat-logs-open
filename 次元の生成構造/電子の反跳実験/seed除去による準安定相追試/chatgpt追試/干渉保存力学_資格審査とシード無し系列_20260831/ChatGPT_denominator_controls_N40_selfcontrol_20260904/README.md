# 新プログラム N=40 自己対照（2026-09-04）

## 目的

新相互作用系プログラム `run_and_plot_N3_N40_mixedseed_20260903.py`（N=40 を含む唯一の新系プログラム）
をコピーし、N=40 について正本結果（2026-09-03 走行）と同一になるかを対照する。
新→旧への一因子変形実験（N=40 ライン）の出発点。

## ファイル

- `run_and_plot_N3_N40_mixedseed_20260903.py` — 正本の bit 同一コピー（diff で確認、無変更）
- `run_and_plot_N40_only_selfcontrol.py` — 上のコピーに対し変更は **2行のみ**（diff 全文を確認済み）:
  1. `OUT` を本フォルダ `results/` に変更（無変更だと正本 results/ を上書き破壊するため必須）
  2. `for N in range(3,41)` → `range(40,41)`（N=40 のみ実行。各 N は独立）

物理・入力・dtype・演算順序は一切変更なし。N=40 の入力は正本と同じ
`hm_mp_free_N3_N40_20260901/data/hm_N40/parent_v.npz` の `v`。

## 実行と合格ゲート（2026-09-04）

`python3 run_and_plot_N40_only_selfcontrol.py` → `done N 40` / `ALL DONE`（exit 0）

- **GATE PASS**: 6分母（38,39,40,41,42,124）の `hm_N40_den_*_states_500.npz` 全配列
  （Z 501×780・N・denominator・steps）が正本と **bit 一致**
- **GATE PASS**: timeseries CSV の N=40 行 3006 行、summary CSV の N=40 行 6 行が正本と完全一致

## 記録事項

- N=40 単独実行では matmul RuntimeWarning（divide by zero / overflow / invalid value）が
  N=40 処理中に表示された。正本の全系列走行では同種警告が N=8 で先に出て以降抑制されて
  いた（numpy は同一箇所の警告を1プロセス1回しか出さない）ため、N=40 でも同じ計算が
  起きていたことが今回可視化された。結果は bit 一致しており挙動差ではない。
- 図 `fig_Hperp_denominator_controls_with_124_N3_N40.png` は N=40 のパネルのみ描画される
  （legend 警告は空パネル由来、無害）。

## 追加読出し図（plot_complex_plane_N40_selfcontrol_v1.py、読み出しのみ）

Δτ=2π/40 の `results/hm_N40_den_40_states_500.npz` から。描画部は
`自発的分裂予備実験_v1/N40_state_readout_20260904/plot_complex_plane_N40_v1.py` と同一様式。

- `fig_N40_selfcontrol_complex_plane_step0.png` — step0（hm_mp_free 解析親）
- `fig_N40_selfcontrol_complex_plane_final.png` — 最大 step（τ=500）
- `fig_N40_selfcontrol_complex_plane_final_cluster_zoom.png` — 凝縮部（角クラスター）拡大

## 次段（予定）

N=40 の初期データを、`次元の生成構造/自発的分裂予備実験_v1_N40対照実験系_20260904` の
make_parent（正本引数）で生成した**静的親ファイル**に差し替える。make_parent は動的には
呼ばず、v・g・Z0 を別々に保存し、Z0 が既存 `states_N00040_delta1e-15_seed0.npz` の Z0 と
bit 一致することを合格条件とする。
