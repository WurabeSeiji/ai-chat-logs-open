# 荷電二体系・解析近似基準解 v1
## 再現手順

## 1. raw data

```bash
python solve_charged_binary_analytic_reference_v1.py
```

明示条件:

```bash
python solve_charged_binary_analytic_reference_v1.py \
  --rows 40001 \
  --r0 50 \
  --rf 20 \
  --mass-ratio 1 \
  --lambda-A 0.30 \
  --lambda-B -0.30
```

生成:

- `charged_binary_analytic_reference_v1_raw.csv`
- `charged_binary_analytic_reference_v1_metadata.json`

基準 raw CSV SHA-256:

```text
f4f49789e0555903f77cd76b34a3a617d6bf0970a5a077ef0e8ccf38cc111689
```

## 2. validation

```bash
python validate_charged_binary_analytic_reference_v1.py
```

生成: `validation_results.json`

基準:

```text
max_relative_dt_dr_error               9.800091824565288e-10
max_relative_dphi_dr_error             9.681061897354175e-11
max_relative_dr_dt_from_energy_balance 5.545165674397686e-16
max_abs_kepler_identity_error          4.440892098500626e-16
max_abs_fgw_over_2forb_minus1          0
max_abs_fem_over_forb_minus1           0
max_abs_energy_balance_residual         9.138831184862806e-12
```

## 3. 図化

plotter は raw CSV のみを読み、軌道を再計算しない。

```bash
python plot_charged_binary_analytic_reference_v1.py
```

生成:

- figure01 relative orbit
- figure02 COM orbits
- figure03 radius vs time
- figure04 radiation power
- figure05 cumulative radiated energy
- figure06 EM/GW power ratio
- figure07 charged 1PN Kepler diagnostic
- figure08 radiation frequencies

各図は SVG/PNG を生成するが、Google Drive には PNG を保存しない。

## 4. 監査規則

1. solver → raw CSV の順で生成する。
2. metadata の raw CSV SHA-256 を確認する。
3. validation を独立に実行する。
4. plotter は raw CSV だけを入力する。
5. plotter から solver を呼ばない。
6. 条件を変えた run で v1 基準ファイルを無断上書きしない。
7. reference branch の式・誤差・読み出しを anonymous generator の状態更新へ戻さない。
