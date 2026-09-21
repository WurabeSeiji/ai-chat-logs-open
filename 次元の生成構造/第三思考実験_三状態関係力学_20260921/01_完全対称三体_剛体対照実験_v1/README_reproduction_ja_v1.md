# E3-B 完全対称三体・剛体対照実験 固定点 v1

日付: 2026-09-21
位置づけ: 第三思考実験の最初の N=3 対照実験。完全対称な三状態から個体差が自発生成されるかを検査する。

## 1. 候補則

```text
X_i[k+1] + X_i[k-1] = tau X_i[k] + epsilon sum_{j != i} X_j[k]
```

三状態 a,b,c は完全グラフ K3 で同じ結合 epsilon を持つ。初期状態は原点を重心とする正三角形で、全置換に対して対称。

period n=31 を固定し、相対モードが

```text
lambda_relative = tau - epsilon = 2 cos(2 pi / 31)
```

となるように tau を設定した。

この候補は最終三体則ではなく、完全対称 N=3 のゼロ点を確定する対照実験である。

## 2. 事前宣言した落ちる条件

- E3B1: 重心閉包が 1e-10 を超えて崩れる。
- E3B2: 三辺長の等値性または剛体性が 1e-10 を超えて崩れる。
- E3B3: 三角形面積の保存が 1e-10 を超えて崩れる。
- E3B4: 二乗読み出し後の三角形の等値性・剛体性が 1e-10 を超えて崩れる。
- E3B5: eps>0 で回帰した6方向結合が入力 epsilon から 1e-10 相対以上ずれる。
- E3B6: 方向付き逆向き比 C_{i<-j}/C_{j<-i} が 1 から 1e-10 以上ずれる。
- E3B7: 状態ごとの半径または一歩面積が 1e-10 相対以上に分裂する。

すべて PASS。

## 3. 実験範囲

```text
epsilon = 0, 1e-5, 1e-4, 1e-3, 5e-3, 1e-2
steps   = 3100 per epsilon
period  = 31
```

## 4. 再現順序

同じフォルダで実行する。

```bash
python3 run_E3B0_symmetric_three_body_rigid_v1.py
python3 analyze_E3B1_symmetric_three_body_v1.py
python3 plot_E3B_symmetric_three_body_v1.py
python3 verify_E3B_fixed_point_v1.py
```

図化プログラムは保存済み CSV のみを読み、数値実験を再実行しない。

## 5. 固定データ

- `E3B0_symmetric_three_body_results_v1.json`: 実験条件、全判定、6方向結合を含む結果。
- `E3B0_symmetric_three_body_summary_v1.csv`: epsilon ごとの集約診断。
- `E3B0_symmetric_three_body_trajectories_v1.csv`: 全 epsilon・全 step の三状態座標、二乗読み出し座標、三辺長、面積、一歩面積。
- `E3B0_run_stdout_v1.txt`: 実行時標準出力。
- `E3B1_symmetric_three_body_analysis_v1.md`: 解析的縮約と物理解釈。

## 6. 解析的結論

重心ゼロなので

```text
X_a + X_b + X_c = 0
sum_{j != i} X_j = -X_i
```

ゆえに候補則は完全対称部分空間で厳密に

```text
X_i[k+1] + X_i[k-1] = (tau - epsilon) X_i[k]
```

へ縮約する。

したがって完全対称三体系は、同じ二値二階則を 120 deg 位相だけずらして走らせる剛体運動となり、個体差は生成しない。

6方向結合は全て同じで、質量比的読みは

```text
m_a : m_b : m_c = 1 : 1 : 1
```

以外を与えない。電荷比的読みも三者の等値以上は区別できない。

## 7. 図

- `fig_E3B_01_rigid_triangle_snapshots_v1.png`: 状態空間での正三角形剛体回転。
- `fig_E3B_02_quadratic_readout_orbits_v1.png`: 二乗読み出し後の3軌道。位相差だけで同じ軌道族。
- `fig_E3B_03_pair_distances_v1.png`: 三辺長の相対変化。数値精度内で一定。
- `fig_E3B_04_symmetry_residuals_v1.png`: 重心閉包・辺等値・面積・個体差の残差。
- `fig_E3B_05_masslike_ratios_v1.png`: 方向付き結合比の 1 からの偏差。

## 8. 固定結論

**完全対称な三体系では、三体になっただけでは質量比・電荷比の非自明な個体差は読めない。比率 1:1:1 がそのまま結果である。**

次段階では、外部から粒子名を入れず、関係または保存量に最小の非対称を入れ、その非対称が方向付き応答比へ写るかを検査する。
