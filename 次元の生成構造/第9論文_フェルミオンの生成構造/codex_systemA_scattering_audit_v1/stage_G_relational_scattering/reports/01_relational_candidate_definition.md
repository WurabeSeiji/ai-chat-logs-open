# relational_C1の数学定義

## 由来別復調と一次元関係波

全状態チャネル \(\Psi_X\)（\(X=A,B\)）から、由来 \(\lambda\in\{A_0,B_0\}\) ごとに正しいηモードへ射影し、対応搬送波を除去したカーネルを

\[
k_X^{(\lambda)}(u)
=D_\lambda^{-1}\Pi_\lambda\Psi_X
\]

とする。

関係量に使う一次元カーネル波は

\[
a(u)=k_A^{(A_0)}(u)+k_A^{(B_0)}(u),
\qquad
b(u)=k_B^{(A_0)}(u)+k_B^{(B_0)}(u)
\]

と定義する。すなわち、由来別に正しく復調した後でコヒーレントに合成する。由来ラベルを直和のまま内積する定義は採用しない。直和内積は全状態のA/Bチャネル直交性をそのまま保存し、既存System Aでは関係量が恒等的にゼロとなるため、本Stageが定義する「復調後の波形関係」と異なる量になる。

## 複素重なりと関係強度

\[
z_{AB}
=
\frac{\langle a,b\rangle}{\|a\|\|b\|},
\qquad
\Gamma_{AB}=|z_{AB}|^2。
\]

Cauchy–Schwarz不等式

\[
|\langle a,b\rangle|\le\|a\|\|b\|
\]

から

\[
0\le\Gamma_{AB}\le1
\]

を得る。

コヒーレント合成波のノルム二乗が `1e-24` 未満の場合、\(\Gamma_{AB}\) は0と置かず数値不成立とする。

## 応答と散乱角

\[
\bar c=\frac{c_A+c_B}{2},
\qquad
F_{\mathrm{rel}}=-\bar c\,\Gamma_{AB},
\]

\[
\Delta\theta_{\mathrm{rel}}
=\kappa\rho(\theta_0)F_{\mathrm{rel}},
\qquad
\theta_{\mathrm{eff}}
=\theta_0+\Delta\theta_{\mathrm{rel}}。
\]

したがって

\[
\boxed{
\theta_{\mathrm{eff}}
=
\theta_0-\kappa\rho(\theta_0)
\frac{c_A+c_B}{2}\Gamma_{AB}
}
\]

である。

## 純状態での符号

- 純奇数同型：\(\bar c=-1\) より \(F_{\mathrm{rel}}=+\Gamma_{AB}\)
- 純偶数同型：\(\bar c=+1\) より \(F_{\mathrm{rel}}=-\Gamma_{AB}\)
- 純異型：\(\bar c=0\) より \(F_{\mathrm{rel}}=0\)

## 分類

- \(\Gamma_{AB}\) とrelational_C1：モデル定義
- \(0\le\Gamma_{AB}\le1\)：数学的帰結
- パリティ符号と関係強度が散乱角を制御する：作業仮説
- \(\Gamma_{AB}\) を自然界の相互作用強度と同一視すること：未導出
