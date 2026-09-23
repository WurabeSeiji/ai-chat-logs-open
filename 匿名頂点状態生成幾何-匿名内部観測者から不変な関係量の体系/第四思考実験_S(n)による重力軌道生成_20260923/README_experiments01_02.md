# 論文4 数値実験01・02：S(n)による重力軌道生成

作成日: 2026-09-23

## 目的

論文3の C1 条件を基準として、同じ初期形状パラメータから S(n) 系列で軌道を生成し、
保存的 S(n) と 2.5PN 放射反作用を追加した S(n) を比較する。

基準条件:

- G = c = M = 1
- a = 30
- b = 27
- q = m_A/m_B = 1
- nu = 0.25
- e0 = sqrt(1-b^2/a^2) = 0.4358898943540673
- p0 = a(1-e0^2) = 24.3
- S(n) 展開次数 N = 12
- 4000 step / radial cycle

論文3基準軌道は `generate_paper3_reference_orbit.py` の同じ 2PN 保存項 + 2.5PN 放射反作用を用いる。

## 実験01：保存的 S(n)

内部状態:

    X_{k+1} = S_k X_k

位置読出し:

    Phi(a,b) = (a^2-b^2, 2ab)

角度ゲージ:

    g_N(chi) = sum_{m=0}^N c_m ((3+e cos chi)/p)^m

    c_0 = 1
    c_{m+1} = ((2m+1)/(m+1)) c_m

半径:

    r = p/(1+e cos chi)

各 S_k は det(S_k)=1 を満たす。

## 実験02：保存的 S(n) + 2.5PN 放射反作用漸化

S(n) の構造は実験01と同一。p,e を固定せず、論文3と同じ局所 2.5PN 放射反作用項を
Newtonian osculating elements に作用させて RK4 で chi に対して更新する。

Paper 3 と同じ放射反作用項:

    C = (8/5) nu / r^3
    an = rdot (18 v^2 + 2/(3r) - 25 rdot^2)
    bv = -(6 v^2 - 2/r - 15 rdot^2)
    a_RR = C (an n + bv v)

Newtonian osculating relations:

    p = h^2
    e^2 = 1 + 2 E p
    dE/dt = v . a_RR
    dp/dt = 2 C bv p
    de/dt = (p dE/dt + E dp/dt)/e

これを dt/dchi = r^2/sqrt(p) と組み合わせ、p(chi), e(chi) を RK4 更新する。
2.5PN の散逸項に Newtonian osculating relations を使うのは最小実験であり、
保存側の高次補正を RR 評価に入れると 3.5PN 以上の補正を含むため、今回は入れていない。

## 実行結果（C1, N=12）

### 論文3基準

- 近点移動: 0.926885 rad/orbit（比較スクリプトでの抽出値）
- 近点半径:
  - start: 16.923303
  - 1回後: 16.825431
  - 2回後: 16.725070

### 実験01

- 近点移動: 0.960214 rad/orbit
- radial RMS / a: 0.0251045
- 近点半径は保存的なので 16.923303 のまま
- max |det S - 1| = 5.55e-16

### 実験02

- p: 24.3 -> 23.704787
- e: 0.4358899 -> 0.4203470
- 近点移動: 0.969957 rad/orbit
- radial RMS / a: 0.0195011
- radial RMS は実験01より 22.32% 改善
- 近点半径:
  - start: 16.923303
  - 1回後: 16.846739
  - 2回後: 16.768800
- max |det S - 1| = 5.55e-16

## 初回評価

2.5PN 放射反作用漸化を追加すると、論文3に見られる軌道縮小を明確に再現し、
半径方向の一致は改善した。一方、近点移動角は実験01よりわずかに過大になった。
これは散逸項の導入自体は有効だが、論文3が有限質量比 nu=0.25 の 1PN+2PN 保存力学であるのに対し、
S(n) 保存部が Schwarzschild 型である差が残っていることと整合的である。

## 実行

    python run_all_experiments01_02.py

## 主な出力

- `reference_pn_C1.csv`
- `experiment01_sn_conservative_C1_N12.csv`
- `experiment02_sn_2p5pn_rr_C1_N12.csv`
- `comparison_experiments_01_02_C1_N12.json`
- `compare_C1_PN_vs_Exp01_Exp02_N12.svg`
- PNG はローカル確認用のみで Google Drive には保存しない。
