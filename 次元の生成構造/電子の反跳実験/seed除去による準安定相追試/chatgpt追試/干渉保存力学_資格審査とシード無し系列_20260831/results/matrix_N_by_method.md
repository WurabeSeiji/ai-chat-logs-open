# 干渉保存力学・シード無し系列：予測 vs 実測（機械判定）

予測一致 38/40

| tag | 予測 | 実測 | 一致 | res_new/r²(親) | ρ−1 | 閉塞max | H⊥率max | 重なり欠損max | PR/M末 | 時計比 | 着地 | μ比 | 有理近似 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| hm_N10 | eq_inflating | inflating | OK | 1.4e-15 | 3.0e-03 | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.200 | 1.000000 | no | -20.0000 | -20/1 |
| hm_N11 | eq_inflating | inflating | OK | 9.0e-16 | 3.8e-03 | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.182 | 1.000000 | YES | -24.7500 | -99/4 |
| hm_N12 | eq_inflating | inflating | OK | 1.4e-15 | 4.7e-03 | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.167 | 1.000000 | YES | -30.0000 | -30/1 |
| hm_N13 | eq_inflating | inflating | OK | 2.8e-15 | 5.7e-03 | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.154 | 1.000000 | YES | -35.7500 | -143/4 |
| hm_N14 | eq_inflating | inflating | OK | 1.9e-15 | 6.8e-03 | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.143 | 1.000000 | YES | -42.0000 | -42/1 |
| hm_N15 | eq_inflating | inflating | OK | 1.9e-15 | 7.9e-03 | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.133 | 1.000000 | YES | -48.7500 | -195/4 |
| hm_N16 | eq_inflating | inflating | OK | 1.4e-15 | 9.2e-03 | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.125 | 1.000000 | YES | -56.0000 | -56/1 |
| hm_N3 | eq_neutral | held | OK | 5.9e-16 | 4.0e-05 | 1.4e-13 | 5.7e-27 | 6.7e-16 | 1.000 | 1.000000 | YES | 1.0000 | 1/1 |
| hm_N4 | eq_neutral | held | OK | 1.1e-15 | 1.6e-04 | 4.4e-11 | 1.5e-22 | 6.7e-16 | 1.000 | 1.000000 | YES | 1.0000 | 1/1 |
| hm_N5 | eq_neutral | drifting | NG | 6.3e-16 | 4.6e-04 | 1.7e-07 | 1.7e-21 | 8.1e-15 | 1.000 | 1.000000 | no | 1.0000 | 1/1 |
| hm_N6 | eq_neutral | drifting | NG | 8.2e-16 | 7.7e-04 | 1.6e-01 | 1.1e-04 | 6.7e-03 | 0.983 | 1.000000 | no | 0.8852 | 8/9 |
| hm_N7 | eq_inflating | inflating | OK | 1.8e-15 | 1.2e-03 | 1.0e+00 | 9.9e-01 | 1.0e+00 | 0.286 | 1.000000 | no | -8.7453 | -35/4 |
| hm_N8 | eq_inflating | inflating | OK | 1.0e-15 | 1.7e-03 | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.250 | 1.000000 | no | -12.0000 | -12/1 |
| hm_N9 | eq_inflating | inflating | OK | 6.7e-16 | 2.3e-03 | 1.0e+00 | 9.9e-01 | 1.0e+00 | 0.222 | 1.000000 | no | -15.7500 | -63/4 |
| ne_N10 | non_equilibrium | inflating | OK | 7.7e-01 | — | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.200 | 1.000001 | YES | -16.9492 | -203/12 |
| ne_N11 | non_equilibrium | inflating | OK | 7.7e-01 | — | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.182 | 0.999309 | YES | -20.9746 | -21/1 |
| ne_N12 | non_equilibrium | inflating | OK | 7.7e-01 | — | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.167 | 1.000000 | YES | -25.4237 | -178/7 |
| ne_N13 | non_equilibrium | inflating | OK | 9.0e-01 | — | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.154 | 1.000000 | YES | -30.2966 | -303/10 |
| ne_N14 | non_equilibrium | inflating | OK | 7.7e-01 | — | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.143 | 1.000000 | YES | -35.5932 | -178/5 |
| ne_N15 | non_equilibrium | inflating | OK | 7.7e-01 | — | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.133 | 1.000028 | YES | -41.3136 | -413/10 |
| ne_N16 | non_equilibrium | inflating | OK | 7.7e-01 | — | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.125 | 1.000000 | YES | -47.4576 | -522/11 |
| ne_N3 | non_equilibrium | inflating | OK | 3.2e-01 | — | 7.7e-01 | 7.7e-01 | 9.9e-01 | 0.909 | 0.999759 | no | 0.4915 | 1/2 |
| ne_N4 | non_equilibrium | inflating | OK | 4.5e-01 | — | 1.0e+00 | 9.6e-01 | 1.0e+00 | 0.488 | 0.999911 | no | -1.2317 | -11/9 |
| ne_N5 | non_equilibrium | inflating | OK | 2.9e-01 | — | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.400 | 0.999945 | no | -3.6715 | -11/3 |
| ne_N6 | non_equilibrium | inflating | OK | 7.7e-01 | — | 1.0e+00 | 9.9e-01 | 1.0e+00 | 0.333 | 0.999976 | no | -5.0847 | -61/12 |
| ne_N7 | non_equilibrium | inflating | OK | 5.6e-01 | — | 1.0e+00 | 9.9e-01 | 1.0e+00 | 0.286 | 0.999974 | no | -8.0290 | -8/1 |
| ne_N8 | non_equilibrium | inflating | OK | 7.7e-01 | — | 1.0e+00 | 9.8e-01 | 1.0e+00 | 0.250 | 1.000030 | no | -10.1695 | -61/6 |
| ne_N9 | non_equilibrium | inflating | OK | 9.6e-01 | — | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.222 | 1.000001 | no | -11.5809 | -139/12 |
| rb_N10 | non_equilibrium | inflating | OK | 1.5e+00 | — | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.200 | 0.997071 | YES | -12.9934 | -13/1 |
| rb_N11 | non_equilibrium | inflating | OK | 1.7e+00 | — | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.182 | 0.999032 | YES | -15.5324 | -171/11 |
| rb_N12 | non_equilibrium | inflating | OK | 1.2e+00 | — | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.167 | 0.999091 | YES | -20.9825 | -21/1 |
| rb_N13 | non_equilibrium | inflating | OK | 1.3e+00 | — | 1.0e+00 | 9.9e-01 | 1.0e+00 | 0.154 | 0.999076 | YES | -24.0492 | -289/12 |
| rb_N14 | non_equilibrium | inflating | OK | 1.3e+00 | — | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.143 | 0.996928 | YES | -28.9408 | -347/12 |
| rb_N15 | non_equilibrium | inflating | OK | 1.5e+00 | — | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.133 | 0.995446 | YES | -31.8968 | -319/10 |
| rb_N16 | non_equilibrium | inflating | OK | 1.8e+00 | — | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.125 | 0.991017 | YES | -34.2032 | -171/5 |
| rb_N5 | non_equilibrium | inflating | OK | 4.9e-01 | — | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.400 | 0.999920 | no | -3.4806 | -7/2 |
| rb_N6 | non_equilibrium | inflating | OK | 1.5e+00 | — | 1.0e+00 | 9.9e-01 | 1.0e+00 | 0.333 | 0.994032 | no | -3.8037 | -19/5 |
| rb_N7 | non_equilibrium | inflating | OK | 9.0e-01 | — | 1.0e+00 | 9.9e-01 | 1.0e+00 | 0.286 | 0.999606 | no | -6.5072 | -13/2 |
| rb_N8 | non_equilibrium | inflating | OK | 1.2e+00 | — | 1.0e+00 | 9.9e-01 | 1.0e+00 | 0.250 | 0.997222 | no | -8.9404 | -107/12 |
| rb_N9 | non_equilibrium | inflating | OK | 1.3e+00 | — | 1.0e+00 | 1.0e+00 | 1.0e+00 | 0.222 | 0.999231 | YES | -11.0315 | -11/1 |
