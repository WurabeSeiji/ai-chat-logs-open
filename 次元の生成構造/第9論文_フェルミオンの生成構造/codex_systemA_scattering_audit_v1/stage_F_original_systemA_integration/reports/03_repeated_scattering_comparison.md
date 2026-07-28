# 03 代表条件の反復比較

`N_A=1,N_B=63,R0=0.55` では両入力が奇数倍音カーネルであり、全128衝突を通じて `c_A,c_B,c_mean=-1` が機械精度内で維持された。したがって反転Candidate 1の `R_eff` は時間変動せず、κごとの定数となった。

| mode | kappa | R_eff | min L gap (collision) | min N_eff gap (collision) | classification |
|---|---:|---:|---:|---:|---|
| C0 | 0.01 | 0.55 | 0.000114791 (110) | 0.708042 (110) | baseline_equivalent |
| reversed_C1 | 0.01 | 0.553889848459 | 4.9746e-05 (29) | 0.306838 (29) | period_and_amplitude_shift |
| reversed_C1 | 0.1 | 0.588721624283 | 2.47058e-05 (44) | 0.152388 (44) | period_and_amplitude_shift |
| reversed_C1 | 1 | 0.886123958197 | 4.38476e-05 (89) | 0.270456 (89) | period_and_amplitude_shift |

これは「型依存応答の実装」が既存の交換運動を別の一定反射率の交換運動へ写したという数値観察である。パリティ純度が交換に伴って変化した、という観察ではない。
