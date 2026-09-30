# 一般 double-Reissner–Nordström 厳密解による質量・電荷比スイープと aa/ab/bb 読出し

Version 2 / 2026-09-29

## 目的

同じ一つの Einstein–Maxwell 厳密解族を用いて、

- 電荷なし
- 質量優勢
- 質量と電荷が近い領域
- 極値領域 |Q|=M
- 電荷優勢

を別ルールへ切り替えずに連続的に扱えるかを確認する。また、一体場の中から自己寄与と二体非線形交差寄与を読み出し、内部複素状態 a,b の二次関係 aa,2ab,bb の spin-weight 2 型変換を同時検証する。

## 使用した厳密解

Manko の物理パラメータ表示による一般 double-Reissner–Nordström 解を実装した。任意パラメータは

(M,m,Q,q,R)

である。

mu = (m Q - M q)/(R+M+m)

Sigma^2 = M^2-Q^2+2 mu Q

sigma^2 = m^2-q^2-2 mu q

nu = R^2-Sigma^2-sigma^2+2 mu^2

kappa = M m-(Q-mu)(q+mu)

Ernst potentials と Weyl metric functions は論文の閉形式 A,B,C から直接評価した。

|Q|=M の等質量・等電荷ケースでは Sigma=sigma=0 となり一般表示が退化するため、同じ解族の厳密極限である二中心 Majumdar–Papapetrou 解を用いた。数値的に 0.999999 等へずらして回避してはいない。

## 実行条件

等質量スイープ:

M=m=1, R=8, Q=q=lambda

lambda = 0, 0.1, 0.5, 0.9, 1.0, 1.1, 2.0, 10.0

評価点:

(rho,z)=(1.3,0.0),(2.0,1.0),(3.5,-0.7)

加えて一般5パラメータ性確認として

M=1.2, m=0.8, Q=0.35, q=-0.15, R=6.5

を評価した。この場合 mu != 0 である。

## 結果

代表点 (rho,z)=(1.3,0) における結果:

| Q/M | Sigma | f | Phi | 分類 |
|---:|---:|---:|---:|---|
| 0.0 | 1.000000 | 0.3803236 | 0 | 電荷なし |
| 0.1 | 0.994987 | 0.3810452 | 0.0309958 | 強い質量優勢 |
| 0.5 | 0.866025 | 0.3987597 | 0.1564275 | 質量優勢 |
| 0.9 | 0.435890 | 0.4431925 | 0.2878491 | 質量・電荷近接 |
| 1.0 | 0 | 0.4593160 | 0.3222715 | extremal / MP exact limit |
| 1.1 | 0.458258 i | 0.4777924 | 0.3575125 | over-extreme |
| 2.0 | 1.732051 i | 0.8042261 | 0.7328650 | 電荷優勢 |
| 10.0 | 9.949874 i | 2.3796967 | -1.0788540 | 強い電荷優勢 |

Q/M>1 では Sigma,sigma は虚数になるが、物理 metric functions f, Phi, exp(2gamma) は実数へ戻る。これは double-RN 解の hyperextreme analytic continuation である。

最大虚部は lambda=10 の強い over-extreme ケースでも約 6e-15 であった。

## Einstein–Maxwell 方程式の数値残差

実装確認のため、静電 Einstein–Maxwell 方程式

 div(f^-1 grad A_t)=0

 f Delta f = |grad f|^2 + 2 f |grad A_t|^2

を有限差分で評価した。

lambda <= 2 では代表点群の相対残差は概ね 1e-6 以下、非対称5パラメータ例では Maxwell 残差 2.7e-8、Einstein scalar residual 3.3e-7 程度であった。

lambda=10 では強い hyperextreme 場のため固定差分幅で 1e-5 程度まで悪化したが、物理量の虚部は機械精度であり、閉形式評価自体は実数性を維持した。

## aa / ab / bb の扱い

一般 double-RN 解は MP のように単純な

U^2 = aa + ab + bb + ...

という多項式ではない。したがって一般領域で「ab=2ab」と先に決めることはしない。

代わりに任意の厳密 observable O について

O_cross = O_full - O_a - O_b + O_vac

という inclusion-exclusion を用いる。

これは

- O_a: a 単体の exact Reissner–Nordström field
- O_b: b 単体の exact Reissner–Nordström field
- O_cross: full nonlinear two-body interaction remainder

を厳密場から読み出す方法である。

したがって一般領域では aa,bb は単体 exact sectors、ab は「二体を同時に置いたことで初めて生じる非線形 remainder」と読むのが自然である。

MP 極限ではこの remainder が明示的な二次交差項へ還元される。

## spin-2 型読出し

内部複素状態

a=A exp(i phi_a), b=B exp(i phi_b)

に共通横断面回転

a -> exp(i psi)a, b -> exp(i psi)b

を施すと

aa=a^2 -> exp(2 i psi) aa

ab=2ab -> exp(2 i psi) ab

bb=b^2 -> exp(2 i psi) bb

となる。

257個の回転角で数値確認した最大誤差:

- aa: 7.0e-16
- ab: 7.0e-16
- bb: 2.8e-16
- total=(a+b)^2: 1.37e-15

従って aa と bb からも spin-weight 2 型読出しは可能であり、ab だけが spin-2 を担うわけではない。

ただしこれは internal symmetric-square の kinematic weight-2 test である。今回の double-RN 解は静的なので、physical radiative gravitational-wave spin-2 を証明したことにはならない。次段階では radiative Weyl curvature / h_+,h_x と直接比較する必要がある。

## 現時点の結論

1. 同じ一般 double-RN 厳密解族で、電荷なしから強い電荷優勢まで parameter change のみで扱える。
2. |Q|=M は数値的回避ではなく exact Majumdar–Papapetrou limit として接続できる。
3. 一般5パラメータ非対称ケース mu != 0 も同じ実装で動作する。
4. 一般 exact field では aa/ab/bb を単純多項式へ固定せず、single-body exact sectors + nonlinear cross remainder として一つの解から読むべきである。
5. aa,ab,bb の symmetric-square は全て同じ weight-2 変換則を持つ。
6. 未解決なのは「放射する重力波の physical spin-2」であり、静的 exact anchor と混同しない。

## 再現ファイル

- double_rn_exact_sweep_v2.py
- double_rn_exact_sweep_v2.json
