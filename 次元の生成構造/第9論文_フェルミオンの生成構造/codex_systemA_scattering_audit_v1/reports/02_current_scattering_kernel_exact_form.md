# 02 現行散乱核の厳密形

## 1. 係数生成

**[コード上の事実]** 外部入力 `R_input` から、

\[
\delta=2\arcsin\sqrt{R_{\mathrm{input}}}
\]

を計算する。

```text
散乱源:140-143
System A:718-720
```

散乱係数は、

\[
t=e^{i\delta/2}\cos\frac{\delta}{2},
\qquad
r=-i e^{i\delta/2}\sin\frac{\delta}{2}
\]

である。

```text
散乱源:134-137
```

よって、

\[
T=|t|^2=1-R_{\mathrm{input}},
\qquad
R=|r|^2=R_{\mathrm{input}}.
\]

さらに \(c=\sqrt{1-R}\), \(s=\sqrt R\) とおけば、

\[
t=(1-R)+i\sqrt{R(1-R)},
\]

\[
r=R-i\sqrt{R(1-R)}.
\]

## 2. 生の二チャネル更新

**[コード上の事実]**

\[
\widetilde a=ra+tb,
\qquad
\widetilde b=ta+rb.
\]

すなわち、

\[
\begin{pmatrix}
\widetilde a\\
\widetilde b
\end{pmatrix}
=
U(R)
\begin{pmatrix}
a\\
b
\end{pmatrix},
\qquad
U(R)=
\begin{pmatrix}
r&t\\
t&r
\end{pmatrix}.
\]

```text
System A:718-736
```

## 3. ユニタリ性

**[数学的帰結]**

\[
|r|^2+|t|^2=1,
\qquad
r^*t+t^*r=0
\]

なので、

\[
U^\dagger U=I.
\]

したがって正規化前の二チャネル結合ノルムは保存される。

\[
\|\widetilde a\|^2+\|\widetilde b\|^2
=
\|a\|^2+\|b\|^2.
\]

**[数値観測]** 独立診断では \(R=0,0.25,0.5,0.6971778791282474,1\) の全点で、ユニタリ誤差は最大 \(2.26\times10^{-16}\) だった。

```text
logs/current_scattering_diagnostic.json
```

## 4. 実コードの最終更新

**[コード上の事実]** 実コードは生の出力をそのまま使わず、各チャネルを別々に正規化する。

\[
a'=\frac{\widetilde a}{\|\widetilde a\|},
\qquad
b'=\frac{\widetilde b}{\|\widetilde b\|}.
\]

したがって実装された写像は、

\[
\mathcal S_R(a,b)=
\left(
\frac{ra+tb}{\|ra+tb\|},
\frac{ta+rb}{\|ta+rb\|}
\right).
\]

これは一般には非線形であり、ユニタリな線形写像 \(U(R)\) そのものではない。

## 5. 現在の既定入力での簡約

**[コード上の事実]** System Aは `hair_enabled=True` を固定し、A側に \(m_A=1\)、B側に \(m_B=2\) の \(\eta\) 位相を与える。

```text
散乱源:42-47
散乱源:88-101
System A:718-723
```

16点の周期 \(\eta\) 格子上でこの2モードは直交する。初期状態がそれぞれ単位ノルムなら、

\[
\langle a,b\rangle=0,
\quad
\|\widetilde a\|=\|\widetilde b\|=1.
\]

また \(U(R)\) がユニタリなので、この直交性と各チャネルノルムは反復後も保存される。

**[数学的帰結]** 現行System Aの既定経路では、各衝突後の `normalize` は丸め誤差を除いて恒等作用である。

## 6. パリティとの可換性

二チャネルへ同じ半周期移動を

\[
\mathcal P=\operatorname{diag}(P,P)
\]

として作用させる。

**[数学的帰結]** \(r,t\) は状態に依存しない複素スカラーなので、

\[
U(R)\mathcal P=\mathcal P U(R).
\]

さらにノルムは \(P\) で不変だから、チャネル別正規化を含む実装でも、

\[
\mathcal S_R(\mathcal P\Psi)=\mathcal P\mathcal S_R(\Psi)
\]

である。

よって現行散乱核は偶数・奇数セクターを相互変換しない。ただし、異なるセクターを持つA/Bを線形混合するため、各出力チャネルのセクター重量は \(R,T\) に従って混ざり得る。

## 7. 判定

**[コード上の事実]** `r,t` は `R_input` だけで決まり、次を参照しない。

- FFT係数
- 倍音番号
- 半周期移動 \(P\)
- 半周期相関
- 偶数・奇数セクター重量
- `HarmonicCase.mode`

**[数学的帰結]** 現行散乱核は型を保存し得るが、型ごとの異なる散乱係数を生成しない。

**[作業仮説]** 「偶数倍音ボゾン型・奇数倍音フェルミオン型」は、この現行散乱核から導出されたものではない。
