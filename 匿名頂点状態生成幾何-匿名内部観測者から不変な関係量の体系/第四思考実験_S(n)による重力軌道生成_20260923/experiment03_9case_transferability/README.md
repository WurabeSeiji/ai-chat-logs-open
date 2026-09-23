# Paper 4 Experiment 03 — 9-case transferability of one S(n) rule

目的は、個々の軌道へ係数をフィットすることではなく、**同一の S(n) 生成則を変更せず**、異なる初期軌道条件へ適用したときに Paper 3 の PN 軌道族を追従できるかを確認することです。

## 条件

- 質量比: `q = 1`（`nu = 0.25`）固定
- `p = 24.3, 60, 120`
- `e = 0.20, 0.45, 0.55`
- 合計 9 ケース
- S(n) 展開次数: `N = 12`
- 4000 steps / radial orbit
- 3 radial-orbit 相当を計算
- ケースごとの再フィットなし。変更するのは初期 `p0,e0` のみ。

## 比較対象

基準: Paper 3 と同じ Newtonian + 1PN + 2PN 保存項 + 2.5PN 放射反作用。

S(n): Experiment 02 と同じ規則。

1. Schwarzschild 型角度ゲージの N 次展開
2. Paper 3 と同じ局所 2.5PN 放射反作用項
3. Newtonian osculating element `(p,e)` の RK4 漸化
4. `det S_k = 1` の内部写像
5. 二次読出し `Phi(X)=(a^2-b^2,2ab)`

## raw データ

`raw/` は**図化より前に保存**されます。

- `pn_integrator_<case>.csv`: PN RK4 の生積分 row
- `pn_plot_<case>.csv`: PN の図化に実際に使用する一様時刻 row
- `sn_<case>.csv`: S(n) の全 step row。`p,e,X0,X1` に加えて 2.5PN RR の局所診断量も保存

`plot_9case_suite.py` は dynamics を再計算せず、この raw CSV だけを読みます。

## 集計

- `summary/case_metrics.csv`
- `summary/case_metrics.json`
- `summary/cases.json`

主指標:

- physical azimuth `phi` を揃えた radial RMS / a
- 近点移動誤差 [rad/orbit]
- 近点半径の縮小率
- `max |det S - 1|`
- internal state map error

## 図

- `experiment03_9case_orbit_grid`: 9条件の軌道比較
- `experiment03_radial_rms_over_a`: RMS誤差
- `experiment03_precession_error`: 近点移動誤差
- `experiment03_peri_shrink_error`: 近点縮小率誤差

各図は PNG と SVG の両方を出力します。

## 実行

```bash
python code/run_all_9case.py
```


### Grid note

当初候補の `p=24.3, e=0.65` は Paper 3 PN 基準生成器の近点速度較正が発散したため、比較基準を同一アルゴリズムで保つ目的で高離心率を `e=0.55` とした。これは S(n) 側の失敗ではなく、基準 PN 数値生成器側の安定範囲による。
