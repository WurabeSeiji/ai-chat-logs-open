# E3-A4 固定結論 — 摂動近似では三状態既約力学は得られない

日付: 2026-09-21

## 位置づけ

第三思考実験へ進む前に、二体 Kepler 系を二階層に置き弱く線形結合する摂動候補を検証した。
この候補は、三状態の既約力学を導出する目的に対して成立するかを判定するための予備実験である。

## 候補則

$$
X^{in}_{k+1}+X^{in}_{k-1}=\tau_{in}X^{in}_k+\varepsilon X^{out}_k,
$$

$$
X^{out}_{k+1}+X^{out}_{k-1}=\tau_{out}X^{out}_k+\varepsilon X^{in}_k.
$$

## 解析結果

階層空間の行列

$$
H=\begin{pmatrix}
\tau_{in} & \varepsilon\\
\varepsilon & \tau_{out}
\end{pmatrix}
$$

は状態に依存しない固定直交変換で厳密に対角化できる。
したがって正常モード $Y_\pm$ は

$$
Y_{\pm,k+1}+Y_{\pm,k-1}=\lambda_\pm Y_{\pm,k}
$$

という二つの独立な二値則に戻る。

従って、この摂動候補は三状態固有の既約力学ではなく、二つの二体近似が弱く結合されただけの系である。

$$
\boxed{
\text{perturbative hierarchical model}
\;\Longrightarrow\;
\text{two coupled two-body sectors}
\;\Longrightarrow\;
\text{two exact normal-mode two-body laws}
}
$$

## 同時に得られた副次結果

局所面積

$$
q_{in,k}=\omega(X^{in}_k,X^{in}_{k+1}),\qquad
q_{out,k}=\omega(X^{out}_k,X^{out}_{k+1})
$$

は

$$
q_{in,k}-q_{in,k-1}=+\varepsilon\,\omega(X^{in}_k,X^{out}_k),
$$

$$
q_{out,k}-q_{out,k-1}=-\varepsilon\,\omega(X^{in}_k,X^{out}_k)
$$

を満たすので、

$$
\boxed{q_{in}+q_{out}=\mathrm{const}}
$$

は厳密保存である。

非縮退では本当の摂動パラメータは

$$
\boxed{g=\frac{\varepsilon}{|\tau_{out}-\tau_{in}|}}
$$

である。

生の外部読み出しは正常モード交差項により一次で単一円錐から外れるが、正常モード自体はそれぞれ二体 Kepler 型の円錐軌道を厳密に保つ。
従って生の軌道変形は、新しい三体力や post-Keplerian 力の導出ではなく、二正常モードの混合として説明できる。

## 本来の目的に対する判定

本来欲しいのは、無名な三状態 $a,b,c$ から直接導出され、固定基底変換で二体セクターへ完全分解できない三状態既約則である。

今回の摂動候補はこの条件を満たさない。
従って、摂動近似を延長しても本来の答えには到達しない。

第三思考実験本体へ持ち帰る新しい落ちる条件は次である。

$$
\boxed{
\text{状態に依存しない固定基底変換で二体セクターへ完全分解できる候補は、三状態既約力学として採用しない。}
}
$$

## 論文に記載する短い結論案

> 階層化した二つの二体 Kepler 系を弱く線形結合する摂動模型を検証した。全体面積保存量は厳密に保存され、局所保存量は二セクター間で交換される。しかし系全体は状態に依存しない固定直交変換によって二つの独立な二値正常モードへ厳密分解できる。したがって、この摂動模型は三状態固有の既約力学を生成せず、二体近似が二つ結合された系に留まる。よって第三思考実験では、固定基底変換で二体セクターへ完全分解できない三状態則を探索する必要がある。

## 再現物

同フォルダに以下を固定保存している。

- 数値実験プログラム
- 解析プログラム
- 完全 JSON データ
- 集約 CSV データ
- 摂動崩壊解析 CSV
- 軌道座標 CSV
- 全11図
- 図化プログラム
- 再現検証プログラム
- requirements
- SHA-256 manifest
- verification logs

詳細は `README_reproduction_ja_v3.md` を参照。
