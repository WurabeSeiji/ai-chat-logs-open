# E3-A 摂動予備実験 固定点 v3 — 解析的分析と全11図を含む完全再現手順

日付: 2026-09-21  
位置づけ: 第三思考実験へ進む前に行った **階層化二体 Kepler 系の弱結合候補の診断実験と、その固定データに対する解析的分析**。  
重要: 本パッケージは **導出された三状態則ではない**。`(a,b)` と `((ab),c)` という階層を先に与えた候補が何をするかを調べた固定点である。

## 1. 固定する候補

```text
Xin[k+1] + Xin[k-1] = tau_in  Xin[k] + epsilon Xout[k]
Xout[k+1] + Xout[k-1] = tau_out Xout[k] + epsilon Xin[k]
```

E3-A0 は数値走査、E3-A1 は固有値と摂動近似の解析、E3-A2 は同一候補・同一初期条件から軌道座標を固定、E3-A3 は **追加走査をせず** A0/A1/A2 の固定データから厳密恒等式・普遍制御量・位相頑健性を解析したものである。

## 2. 実行環境

```text
Python 3.13.5
NumPy 2.3.5
Matplotlib 3.10.8
```

依存は `requirements_E3A_fixed_point_v1.txt` に固定している。

## 3. 完全再現順序

同じフォルダで次を実行する。

```bash
python3 run_E3A0_hierarchical_two_kepler_perturbation_v1.py
python3 analyze_E3A1_perturbation_breakdown_v1.py
python3 plot_E3A_fixed_point_v1.py
python3 export_E3A2_orbit_trajectories_v1.py
python3 plot_E3A_orbits_v1.py
python3 analyze_E3A3_analytic_v1.py
python3 plot_E3A3_analysis_v1.py
python3 verify_E3A_fixed_point_v3.py
```

図化プログラムは保存済み CSV / JSON を読む。E3-A3 の時系列は、固定済み A0 の `epsilon=0.004`, phase 0 の既存条件を決定論的に再生したもので、新しいパラメータ走査ではない。

## 4. 固定データ

### E3-A0 数値実験

- `E3A0_hierarchical_two_kepler_perturbation_full_v1.json`
- `E3A0_hierarchical_two_kepler_perturbation_aggregates_v1.csv`
- `E3A0_run_stdout_v1.txt`

### E3-A1 摂動崩壊解析

- `E3A1_perturbation_breakdown_scan_v1.csv`
- `E3A1_perturbation_breakdown_analysis_v1.md`
- `E3A1_analysis_stdout_v1.txt`

### E3-A2 軌道固定データ

- `E3A2_orbit_trajectories_v1.csv`
- `export_E3A2_orbit_trajectories_v1.py`

### E3-A3 解析的分析

- `E3A3_analytic_analysis_v1.md`: 数式を含む解析本文。
- `E3A3_analytic_summary_v1.json`: 解析量、境界値、厳密式との照合結果。
- `E3A3_area_exchange_timeseries_v1.csv`: 代表固定点での局所面積交換時系列。
- `E3A3_phase_robustness_v1.csv`: 8位相の min / median / max。
- `E3A3_analysis_stdout_v1.txt`
- `analyze_E3A3_analytic_v1.py`
- `plot_E3A3_analysis_v1.py`

## 5. E3-A3 で解析的に確定した事項

### 5.1 全体面積保存は厳密

`q_in = omega(Xin_k, Xin_{k+1})`, `q_out = omega(Xout_k, Xout_{k+1})` とすると、更新則から直接

```text
q_in(k)  - q_in(k-1)  = + epsilon * omega(Xin_k, Xout_k)
q_out(k) - q_out(k-1) = - epsilon * omega(Xin_k, Xout_k)
```

したがって `q_total = q_in + q_out` は厳密保存で、局所面積変化は等量反対向きの交換である。

### 5.2 非縮退摂動の本当の制御量

```text
g = epsilon / |tau_out - tau_in|
```

である。混合角と外部正常モード中の内部重量は g の普遍関数になる。

### 5.3 二次 law-shift 近似の誤差も g だけで決まる

```text
R(g) = sqrt(1/4 + g^2) - 1/2
```

相対誤差 `R=r` の境界は `g=sqrt(r(1+r))`。

### 5.4 生の外部軌道の一次逸脱

生の外部状態は二正常モードの混合であり、二乗読み出しには `2 sin(theta) cos(theta) y_- y_+` の交差項が入る。そのため raw conic residual は非縮退小結合で一次に出る。

### 5.5 縮退は階層解釈に対して特異

`tau_in=tau_out` では任意の非零 epsilon で mixing angle = 45 deg。exact system は壊れないが、inner / outer を弱く混ざる別セクターとして扱う摂動階層は直ちに失われる。

### 5.6 三状態既約性について得た落ちる条件

この候補は固定直交変換で二つの二値則に完全分解できる。したがって真の三状態候補については、**状態に依存しない固定基底変換で二体セクターへ完全分解できる候補は、三状態既約力学としては落とす**、という判定条件を得た。

## 6. 全11図

1. `fig_E3A_01_raw_hierarchical_response_v1.png`
2. `fig_E3A_02_invariants_and_normal_modes_v1.png`
3. `fig_E3A_03_outer_law_shift_v1.png`
4. `fig_E3A_04_approximation_error_and_mixing_v1.png`
5. `fig_E3A_05_outer_mode_type_boundary_v1.png`
6. `fig_E3A_06_orbit_raw_outer_readout_v1.png`
7. `fig_E3A_07_orbit_outer_normal_mode_v1.png`
8. `fig_E3A_08_local_area_exchange_timeseries_v1.png` — 局所面積の等量交換と全体保存。
9. `fig_E3A_09_phase_robustness_v1.png` — 8位相に対する小結合一次応答の頑健性。
10. `fig_E3A_10_universal_g_control_v1.png` — g による近似誤差と混合重量の普遍曲線。
11. `fig_E3A_11_degenerate_vs_nondegenerate_mixing_v1.png` — 縮退と非縮退の質的差。

## 7. 完全性・再現性

- `verify_E3A_fixed_point_v3.py` は一時ディレクトリで A0/A1/A2/A3 と全11図を再生成する。
- `verification_v3.txt` に固定時照合結果を保存する。
- `SHA256SUMS_v3.txt` に v3 の主要ファイルの SHA-256 を保存する。

固定時、E3-A0/A1/A2/A3 の全固定データおよび図1〜11はすべて `MATCH` した。
