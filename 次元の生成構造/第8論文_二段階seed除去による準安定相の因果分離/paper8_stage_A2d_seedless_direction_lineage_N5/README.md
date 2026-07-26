# Stage A2d: 保存済み完全無seed N=5方向系譜

Stage A2aで取得済みの完全無seed軌道

```text
Z0 = v
initial kernel seed = OFF
metastable transverse seed = OFF
Benettin reinjection = OFF
```

を再実行せず、保存済み `dominant_plane_steps.npy` と
`dominant_plane_values.npy` から、既存論文7定義

```text
D34(t) = s4_new_dirs(B0, Bdom(t))
P34(t) = D34(t) D34(t)^T
```

を再構成する後処理ラッパである。

## 実行

```bash
python3 analyze_saved_seedless_lineage.py
```

## 入力

- `../paper8_stage_A2a_seedless_N5/raw/A2a_N5_seedless_f64_e1/`
- `../paper8_stage_A2a_seedless_N5/raw/A2a_N5_seedless_f64_e2/`
- 比較参照のみ:
  `../paper8_stage_A2c_direction_lineage_N5/processed/early_vs_late_direction_overlap.csv`

A2aの2実行は全入力ファイルをバイト比較し、方向平面配列もbitwise比較する。
`B0` は、A2a実装で同じ親状態から同じ関数により再計算・保存された
step 0の `Bdom` を使用する。

## 標本上の制限

A2aの方向平面は原則5 step間隔で、f指定水準の初回通過stepが追加保存されている。
したがって、早期対後期の主角・overlap・射影距離は直接計算できるが、
A2cの全1-step連続回転最大値とは同一統計ではない。

このパッケージでは、標本間主角と `angle / delta_step` を記述量として保存する。
後者は未保存stepの最大回転を復元する量ではない。

## 出力

- `processed/seedless_sampled_lineage_timeseries.csv`
- `processed/seedless_sampled_continuity.csv`
- `processed/seedless_lineage_by_f_level.csv`
- `processed/seedless_lineage_by_q_band.csv`
- `processed/seedless_lineage_summary.json`
- `figures/figure01_seedless_early_vs_late_lineage.{png,svg}`
- `figures/figure02_seedless_sampled_rotation.{png,svg}`
- `reports/stage_A2d_seedless_direction_lineage_N5_report.md`
- `logs/input_manifest.json`

既存A2cのseedあり `D34/P34` が利用可能な場合は、同じ保存stepで
無seed対seedあり部分空間も直接比較する。これは、それぞれの後期代表方向への
overlap曲線を比較する量とは別である。
