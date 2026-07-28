# 反転Candidate 1の一定パリティ区間における退化

## 対象

反転Candidate 1を

\[
\theta_{\mathrm{eff}}(n)
=
\theta_0-\kappa\rho(\theta_0)
\frac{c_A(n)+c_B(n)}{2}
\]

とする。

## 証明

反復区間 \(I\) の全ての \(n\) について

\[
c_A(n)=c_A^\ast,\qquad c_B(n)=c_B^\ast
\]

と仮定する。このとき

\[
\bar c(n)=\frac{c_A^\ast+c_B^\ast}{2}
\]

は \(n\) に依存しない。したがって

\[
\theta_{\mathrm{eff}}(n)
=
\theta_0-\kappa\rho(\theta_0)
\frac{c_A^\ast+c_B^\ast}{2}
\equiv\theta^\ast
\]

も区間 \(I\) で一定である。散乱係数

\[
t^\ast=e^{i\theta^\ast}\cos\theta^\ast,\qquad
r^\ast=-ie^{i\theta^\ast}\sin\theta^\ast
\]

も一定になるため、この区間の更新はC0を一定角 \(\theta^\ast\)、すなわち一定反射率

\[
R^\ast=\sin^2\theta^\ast
\]

で実行することと等価である。

純奇数対では \(c_A^\ast=c_B^\ast=-1\) なので

\[
\theta^\ast=\theta_0+\kappa\rho(\theta_0)。
\]

純偶数対では \(c_A^\ast=c_B^\ast=+1\) なので

\[
\theta^\ast=\theta_0-\kappa\rho(\theta_0)。
\]

以上により、反転Candidate 1は一定パリティ区間では動的状態依存散乱ではなく、C0の一定反射率への再パラメータ化に退化する。

## 分類

- 退化命題：数学的帰結
- 純奇数・純偶数の符号方向：モデル定義からの導出済み帰結
- これを自然界のフェルミオン・ボゾン散乱と同一視すること：未導出

この証明後、Stage Gでは反転Candidate 1の追加R掃引を行わない。
