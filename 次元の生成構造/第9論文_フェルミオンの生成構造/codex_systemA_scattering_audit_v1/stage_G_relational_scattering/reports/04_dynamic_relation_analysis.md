# 関係量の動的性

## 機械判定

- 動的閾値: `1.0e-10`
- relational_C1の最大 `Delta Gamma`: `3.7123082385903672e-15`
- relational_C1の最大 `Delta R_eff`: `4.4408920985006262e-16`
- 判定: `gamma_constant`, `R_eff_constant`
- 相関54件のうち定数系列として未定義: `54`

未定義相関を0へ置換していない。

## なぜ一定になったか

復調後の初期関係波を \(a_0,b_0\) とし、両者が単位ノルムで実重なり \(s=\langle a_0,b_0\rangle\) を持つと、Gram行列は

\[
G_0=
\begin{pmatrix}1&s\\s&1\end{pmatrix}
=I+s\sigma_x。
\]

System Aの各衝突行列は

\[
U_n=
\begin{pmatrix}r_n&t_n\\t_n&r_n\end{pmatrix}
=r_n I+t_n\sigma_x
\]

である。よって \([G_0,U_n]=0\)。さらに \(U_n\) はユニタリなので

\[
G_{n+1}=U_nG_nU_n^\dagger=G_n。
\]

この帰結は \(U_n\) の角度が状態依存でも成立する。したがって本Stageの二条件では、\(\Gamma=|s|^2\) は保存量となった。代表条件では `Gamma=1/32`、31系列custom packetでは `Gamma=1/2` である。

## 中心判定

relational_C1は今回の既存System A条件では `constant_relation_reparameterization` に退化した。関係量を追加したというモデル定義は成立するが、「純パリティ区間で散乱率を動的にする」という作業仮説はこの対称更新・初期Gram条件では実現しなかった。新候補は生成しない。
