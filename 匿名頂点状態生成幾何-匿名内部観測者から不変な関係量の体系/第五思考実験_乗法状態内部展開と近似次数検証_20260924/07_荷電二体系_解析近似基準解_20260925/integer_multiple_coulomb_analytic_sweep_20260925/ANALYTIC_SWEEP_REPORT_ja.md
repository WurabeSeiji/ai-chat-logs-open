# 整数倍クーロン解析系列：現行解析法そのままの事前実験

## 固定した実験条件

- 解析プログラム: `solve_charged_binary_analytic_reference_v1.py` を無改変で使用。
- `G=M=c=1`, `mass_ratio=1`, `r0=50`, `rf=20`, `rows=40001` は現行実験のまま固定。
- 変更したのは `lambda_A`, `lambda_B` のみ。
- 基準 `|lambda|=0.30` に対し、強度変更は整数倍 `n=1,3,4` のみ。
- 引力: `(lambda_A,lambda_B)=(+0.30n,-0.30n)`。
- 斥力: `(lambda_A,lambda_B)=(+0.30n,+0.30n)`。
- 途中で現行解析法が失敗してもバッチ全体は停止せず、失敗をそのまま記録。

- solver SHA-256: `c3871bdcda6d1701437c73b6f8aea1bc5184f6e8318d528c6e3f7c1a84f23a4f`

## ケースと実行結果

| case | n | lambda_A | lambda_B | |Fc/Fg| | Q | Z=1-Q | current method |
|---|---:|---:|---:|---:|---:|---:|---|
| n1_attractive | 1 | 0.30 | -0.30 | 0.09 | -0.09 | 1.09 | success |
| n1_repulsive | 1 | 0.30 | 0.30 | 0.09 | 0.09 | 0.91 | success |
| n3_attractive | 3 | 0.90 | -0.90 | 0.81 | -0.81 | 1.81 | success |
| n3_repulsive | 3 | 0.90 | 0.90 | 0.81 | 0.81 | 0.19 | success |
| n4_attractive | 4 | 1.20 | -1.20 | 1.44 | -1.44 | 2.44 | success |
| n4_repulsive | 4 | 1.20 | 1.20 | 1.44 | 1.44 | -0.44 | failed_current_method |

## 成功ケースの現行解析法出力

### n1_attractive

- `t_final = 169032.317607064`
- `phi_final = 767.088999401591`
- `orbital_cycles = 122.086006046179`
- `P_EM/P_GW initial = 1.72018348623853`
- `P_EM/P_GW final = 0.688073394495413`
- `max_abs_energy_balance_residual = 9.13883e-12`

### n1_repulsive

- `t_final = 574545.646661031`
- `phi_final = 2287.91263636999`
- `orbital_cycles = 364.132605440694`
- `P_EM/P_GW initial = 0`
- `P_EM/P_GW final = 0`
- `max_abs_energy_balance_residual = 9.3577e-12`

### n3_attractive

- `t_final = 17450.0318499065`
- `phi_final = 104.752681063412`
- `orbital_cycles = 16.6719069933707`
- `P_EM/P_GW initial = 9.3232044198895`
- `P_EM/P_GW final = 3.7292817679558`
- `max_abs_energy_balance_residual = 1.27326e-11`

### n3_repulsive

- `t_final = 13179536.0110803`
- `phi_final = 23981.2117594408`
- `orbital_cycles = 3816.72839284849`
- `P_EM/P_GW initial = 0`
- `P_EM/P_GW final = 0`
- `max_abs_energy_balance_residual = 1.95381e-12`

### n4_attractive

- `t_final = 7507.62015701966`
- `phi_final = 52.4575417577687`
- `orbital_cycles = 8.34887707319841`
- `P_EM/P_GW initial = 12.2950819672131`
- `P_EM/P_GW final = 4.91803278688525`
- `max_abs_energy_balance_residual = 1.68392e-11`

## 現行解析法が停止したケース

### n4_repulsive

- `lambda_A=lambda_B=+1.20`
- `Q=+1.44`
- `Z=1-Q=-0.44`
- 現行 solver は次の既存条件チェックで停止した:
  `Z=1-lambda_A*lambda_B must be >0 for an attractive circular benchmark`
- このチェックは削除・回避していない。現行方式での不成立結果として、そのまま保存した。

## 注意

この系列は、新しい一般解法へ切り替えたものではない。現行の leading-order adiabatic quasi-circular analytic benchmark を、パラメータだけ変更して機械的に適用した結果である。
