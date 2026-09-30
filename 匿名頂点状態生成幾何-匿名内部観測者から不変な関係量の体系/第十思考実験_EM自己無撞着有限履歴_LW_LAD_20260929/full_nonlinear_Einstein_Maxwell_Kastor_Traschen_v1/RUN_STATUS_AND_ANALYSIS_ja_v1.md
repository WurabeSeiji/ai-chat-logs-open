# full nonlinear Einstein-Maxwell 二体系：実行状況と結果

## 1. 今回実際に実行したもの

今回この計算環境で実際に実行したのは、二中心 Kastor-Traschen (KT) 厳密解である。

\[
 ds^2=-\Omega^{-2}dt^2+a(t)^2\Omega^2(dx^2+dy^2+dz^2),
\]

\[
 \Omega=1+\sum_A \frac{M_A}{a(t)r_A},\qquad a(t)=e^{Ht},\qquad \Lambda=3H^2.
\]

各中心は \(|Q_A|=M_A\) を持つ。今回は

\[
M_1=M_2=1,\quad Q_1=Q_2=1,\quad d=8,\quad H=-0.05
\]

とした。\(H<0\) は収縮枝であり、二中心は背景 de Sitter スケールに対して接近する。

KT は Einstein-Maxwell-\(\Lambda\) 方程式の full nonlinear exact solution であり、線形摂動ではない。

## 2. 3+1 変数

\[
W=a\Omega=a+\sum_A\frac{M_A}{r_A}
\]

と置くと

\[
\gamma_{ij}=W^2\delta_{ij},\qquad \alpha=\frac{a}{W},\qquad \beta^i=0,
\]

\[
K_{ij}=-H\gamma_{ij},\qquad K=-3H.
\]

電場は符号規約を除いて

\[
E_i=-\partial_i\ln W,
\qquad
E^i=-\frac{\partial_i W}{W^3}.
\]

## 3. 制約式の数値検証

格子上で Hamiltonian constraint と Maxwell Gauss constraint を独立に有限差分評価した。

Hamiltonian constraint:

\[
{}^{(3)}R+K^2-K_{ij}K^{ij}-16\pi\rho_{EM}-2\Lambda=0.
\]

Gauss constraint:

\[
D_iE^i=0
\]

（点電荷位置を除く）。

格子サイズは N=33,49,65。中心から半径 1.2 未満と外側境界1セルを評価から除外した。

### t=0 の RMS residual

| N | h | Hamiltonian RMS | Gauss RMS |
|---:|---:|---:|---:|
| 33 | 0.75 | 5.5740e-3 | 1.3935e-3 |
| 49 | 0.50 | 1.9527e-3 | 4.8819e-4 |
| 65 | 0.375 | 1.1803e-3 | 2.9507e-4 |

推定収束次数は successive grids で約 2.59, 1.75。点源近傍を固定物理半径で除外した単純2次差分として、概ね2次収束と整合する。

Momentum constraint は

\[
K_{ij}=-H\gamma_{ij}
\]

かつ Poynting flux がゼロなので解析的に 0 である。

## 4. 動的二中心

背景物理距離 \(a(t)d\) は

| t | a(t) | a(t)d |
|---:|---:|---:|
| 0 | 1.000000 | 8.000000 |
| 4 | 0.818731 | 6.549846 |
| 8 | 0.670320 | 5.362560 |

へ減少した。

したがって今回の実行は、静的 double-RN ではなく full nonlinear かつ動的な二中心 Einstein-Maxwell 解を実際に評価している。

## 5. aa, ab, bb の非線形同時存在

中点で

\[
W^2=(a+u_a+u_b)^2
\]

を展開すると

\[
W^2=a^2+2au_a+2au_b+u_a^2+2u_au_b+u_b^2.
\]

したがって

\[
aa=u_a^2,\qquad ab=2u_au_b,\qquad bb=u_b^2
\]

は一つの full nonlinear metric factor 内に同時に存在する。

直接展開値と \(W^2\) の差は全時刻で最大

\[
2.22\times10^{-16}
\]

であった。

## 6. 重要な境界：今回まだ実行していないもの

今回 **実行していない** のは、asymptotically-flat な一般 charged binary merger を 3D numerical relativity で時間発展し、\(\Psi_4\) と Maxwell Newman-Penrose scalar を null infinity / extraction sphere で抽出する計算である。

理由は物理ではなく現在の実行環境である。このセッションには Einstein Toolkit / GRChombo / SpECTRE と MPI runtime がインストールされておらず、利用可能資源は 5 CPU core、約 5.8 GiB RAM である。production NR の charged binary evolution をここで実行済みと主張することはできない。

Einstein Toolkit の現行 Hypatia release は 2026-07-10 公開で、BSSN/Z4c spacetime evolution と matter support を含む。公式 gallery の比較的大きな merger 例では 256 cores / 140 GB 程度が記載されており、一般3D merger はこの環境の規模を大きく超える。

## 7. 本研究に対する現在の到達点

実行済み：

1. static exact Einstein-Maxwell: RN / MP / double-RN
2. dynamic exact full nonlinear Einstein-Maxwell-Lambda: two-center Kastor-Traschen（今回）
3. spin-2 radiation on exact Schwarzschild background: RW/Zerilli
4. coupled spin-2/spin-1 radiation on exact RN background: RN coupled perturbation
5. retarded Maxwell + Lorentz-Dirac finite-history evolution

未実行：

6. asymptotically-flat full nonlinear charged two-body merger in 3D numerical relativity

したがって、full nonlinear Einstein-Maxwell 自体への対応は今回実行した。しかし「一般 charged binary の放射まで含む full nonlinear 3D evolution」は依然として別計算資源を要する。

## 8. 参考

- D. Kastor and J. Traschen, Cosmological Multi-Black Hole Solutions, Phys. Rev. D 47, 5370 (1993), arXiv:hep-th/9212035.
- Einstein Toolkit, current release Hypatia (2026-07-10), https://einsteintoolkit.org/
