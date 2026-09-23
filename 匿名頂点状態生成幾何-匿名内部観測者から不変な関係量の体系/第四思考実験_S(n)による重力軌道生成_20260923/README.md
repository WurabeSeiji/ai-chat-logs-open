# 第四思考実験 数値実験01 — S(n)による重力軌道生成と論文3 PN軌道の比較

日付: 2026-09-23

## 目的

論文3の強い場二体軌道生成器を基準とし、同じ初期軌道形状 C1 (a=30, b=27, m_A/m_B=1) に対して、

- 基準: 論文3と同じ Newtonian + 1PN + 2PN 保存項 + 2.5PN 放射反作用
- 比較: `X_{k+1}=S_k X_k`, `det S_k=1` と論文2の二次読出し `Phi(X)=(a^2-b^2,2ab)`

から生成した軌道を重ねて比較する。

## S(n) モデル

Schwarzschild の Darwin パラメータを用いる。

```text
r(chi) = p / (1 + e cos chi)
```

角度ゲージは

```text
dphi/dchi = (1 - 2(3+e cos chi)/p)^(-1/2)
```

を N 次まで二項展開する。

```text
c_0 = 1
c_{m+1} = ((2m+1)/(m+1)) c_m

g_N(chi) = sum_{m=0}^N c_m ((3+e cos chi)/p)^m
```

内部状態は

```text
X_k = sqrt(r_k) [cos(phi_k/2), sin(phi_k/2)]^T
```

とし、一ステップを

```text
S_k = R(phi_{k+1}/2) diag(rho_k, 1/rho_k) R(-phi_k/2)
rho_k = sqrt(r_{k+1}/r_k)
```

で更新する。各ステップで `det(S_k)=1` である。位置は論文2と同じ二次写像で読む。

## 初回実行条件

- C1: a=30, b=27, mass ratio=1, nu=0.25
- 3 radial orbits
- S(n) expansion order N=12
- 4000 steps / radial orbit

## 初回結果

- 論文3 PN 近点移動: 約 0.92689 rad/orbit
- S(n), N=12 近点移動: 約 0.96021 rad/orbit
- 差: 約 0.03333 rad/orbit (約 1.91 deg/orbit)
- 半径差 RMS / a: 約 2.51e-2
- max |det(S)-1|: 約 5.55e-16
- S による内部状態更新の最大数値誤差: 約 1.56e-13

基準PNは2.5PN放射反作用を含むため軌道が徐々に縮小する。一方、今回の最小 S(n) モデルは保存的 Schwarzschild 型であり、その差が後半の軌道差として現れる。したがって本結果は完全一致ではなく、S(n) による保存的重力軌道が論文3の強い場軌道へどの程度近づくかを見る最初の比較である。

## ファイル

- `generate_paper3_reference_orbit.py` — 論文3の生成器のコピー
- `run_reference_pn_orbit.py` — C1基準軌道だけを生成
- `run_sn_gravity_orbit.py` — 任意次数 N の S(n) 軌道を生成
- `compare_orbits.py` — 二軌道を比較し SVG/PNG を生成
- `run_all.py` — 上記を順に実行
- `results/reference_pn_C1.csv`
- `results/reference_pn_C1_meta.json`
- `results/sn_gravity_C1_N12.csv`
- `results/sn_gravity_C1_N12_meta.json`
- `results/comparison_C1_N12.json`
- `figures/compare_C1_PN_vs_Sn_N12.svg`
- `figures/compare_C1_PN_vs_Sn_N12.png` — ローカル生成のみ。Google Drive には保存しない。

## 再実行

```bash
python run_all.py
```
