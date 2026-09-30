# 厳密ルール再構築 実験結果 v1

## 結論

2026-09-29 の監査で、旧 EM-harmonic, RWZ, RN-coupled, Kastor--Traschen 実装のうち、本研究固有の `STRICT_EXPERIMENT_RULES_ja` に照らして不適合な部分を確認したため、状態生成系を全面再構築した。

新実装では、永続物理状態は有限履歴を表す **複素倍音係数だけ**である。全次状態は、固定 sparse quadratic interaction array による

\[
Z_{n+1}=\mathcal T[Z_n,Z_n]
\]

の一回積和で同期生成される。

実行 transition の実コードは RN, KT, spin の全てで本質的に

```python
mon = z[pi] * z[pj]
z_next = M @ mon
```

だけである。

## 再実行した実験

### A. RN coupled anonymous two-channel radiation

同一 interaction array、同一コードで初期状態の Q だけを変更した。

- Q/M = 0, 0.3, 0.6, 0.9, 0.99
- anonymous channel 0 initial / anonymous channel 1 initial の両方
- Nr=61, harmonic history K=4, 700 steps

Q=0 では other-channel flux は完全に0。
Q>0 では連続的に他チャネルへの変換が発生した。

channel 0 initial:

- 0.3: 0.0016139064
- 0.6: 0.0062901234
- 0.9: 0.0137090482
- 0.99: 0.0164396584

channel 1 initial:

- 0.3: 0.0014597175
- 0.6: 0.0057071048
- 0.9: 0.0125027843
- 0.99: 0.0150217465

### B. Q=0 放射・吸収バランス

同一 strict map で初期波束位置を探索した。

\[
x_0\approx 2.97890625
\]

で

\[
(F_{out}-F_{in})/(F_{out}+F_{in})
\approx -6.44\times10^{-5}.
\]

### C. harmonic compilation 独立照合

strict harmonic map と、外部 validation 用 direct sample recurrence を120 step 比較。

\[
\max|\Delta|=3.0654\times10^{-17}.
\]

### D. Kastor--Traschen full nonlinear exact sector

外部 `t -> exact state` 注入を廃止し、scale state 自身を

\[
a' = a h
\]

として同一 quadratic map で更新した。

32 step まで外部 exact benchmarkとの差は 0--1.1e-16。
Hamiltonian/Gauss residual は N=33,49,65 の格子細分化で減少。
aa/ab/bb 再構成誤差は0。

### E. spin hierarchy

同じ一回 quadratic map で half-angle state と ordinary-angle state を更新。

- 2pi sign flip error = 1.16e-14
- 4pi return error = 2.40e-14
- weight-2 covariance error = 4.48e-16

## 旧実装の扱い

以下は削除しないが、本研究の strict universal-interaction implementation としては superseded とする。

- `em_retarded_harmonic_lad_v3.py`: persistent state は倍音だが内部 RK4 が R4 不適合。
- `rwz_spin2_flux_v1.py`: leapfrog previous/current array が strict state rule 不適合。
- `rn_coupled_gem_v1.py`: 同上。
- `kastor_traschen_two_center_v1.py`: exact-state sampling であり state-generated dynamics ではない。

これらは標準理論の reference/benchmark としてのみ保持する。

## 未実行範囲

一般の asymptotically-flat full nonlinear charged binary 3D numerical relativity は、現在の計算資源では未実行である。本再構築はそれを実行したと主張しない。
