# Nall_linear1000000_steps5_mpmath50_equimodular_selfconsistent_directHperp_20260828

刻み 2π/1000000・5 step を **多倍長（mpmath 50 桁）**で走らせ、倍精度の読出し床（10⁻³¹）の下で H⊥ がどうなるかを見る（N=5,8,10,16,20）。親は各 N の倍精度パッケージ `…linear100000_steps50…/data/states_treatment.npz` の Z[0]（倍精度で求めた自己無撞着親、残差 r は倍精度水準のまま）。相互作用・回転・読出しは同じ定義を mp で再実装（`program/run_mp.py`、exp(ΔK) は Taylor 10 項）。結果は `results/summary.json`、`data/N*_mp_timeseries.csv`、`figures/compare_N_L1000000_5_mp50.png`、`実験結果_…md`。

```bash
python3 program/run_mp.py   # 数分
```
