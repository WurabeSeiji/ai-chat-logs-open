# 厳密ルール再構築：倍音主状態・一回積和写像 実装監査

## 0. 目的

本書は、2026-09-29 の再監査で確認されたルール逸脱を修正し、状態生成系を以下の絶対条件へ戻した実装の監査記録である。

- 永続状態は **倍音係数だけ**。
- 全次状態は、更新前の同一状態だけを読む。
- 次状態全体は、同一の固定 interaction array による **一回の積和写像**で生成する。
- 内部に RK stage、`state + increment`、exact-state injection、読み出し値の帰還を残さない。
- ケース依存の電荷量は transition の外部引数にせず状態に含める。
- GW / EM / spin / aa / ab / bb は生成側のラベルではなく readout とする。

基準文書：`STRICT_EXPERIMENT_RULES_ja` の R1--R14。

---

# 1. 旧実装で確認された違反

## 1.1 `em_retarded_harmonic_lad_v3.py`

倍音係数は persistent state だったが、内部更新に RK4 を使用していた。

```python
k1,d1=deriv(...)
k2,_=deriv(...)
k3,_=deriv(...)
k4,_=deriv(...)
out=tuple(s + dt*(a+2*b+2*c+d)/6 ...)
```

したがって R4 の「内部を含めた乗法形式」に不適合である。

## 1.2 `rwz_spin2_flux_v1.py` / `rn_coupled_gem_v1.py`

`Gprev, G, Eprev, E` などが persistent state として保持され、

```python
Gn = 2*G - Gprev + ...
En = 2*E - Eprev + ...
```

という通常の加法 leapfrog で更新されていた。標準方程式の数値解としては有効だが、R2/R4 の本研究固有の生成規則には不適合である。

## 1.3 `kastor_traschen_two_center_v1.py`

```python
for t in times:
    a=exp(H*t)
    W=a+Ustatic
```

として外部時刻から exact state を直接生成していた。これは exact solution evaluator としては正しいが、状態から次状態を生成する universal-interaction implementation ではない。

以上の旧コードは参照計算として保存するが、本研究の「厳密ルール準拠生成器」としては **superseded** とする。

---

# 2. 新実装の共通構造

状態を

\[
Z_n=\{C_{v,m}^{(n)}\}
\]

とする。`v` は匿名状態チャネル、`m` は有限履歴の倍音番号である。

時間サンプル履歴を一切 persistent にせず、

\[
s_{v,h}=\sum_m C_{v,m}e^{2\pi i mh/K}
\]

で表現される有限履歴全体を倍音係数で保持する。

sample-space の合法な積項

\[
s'_{o,h_o}\;{+}{=}\;c\,s_{i,h_i}s_{j,h_j}
\]

を、初期化時に Fourier 基底へ完全展開する。

すると、実行時の状態更新は固定 sparse interaction array \(\mathcal T\) による

\[
\boxed{
Z'_{\alpha}
=
\sum_{\beta,\gamma}
\mathcal T_{\alpha\beta\gamma}
Z_{\beta}Z_{\gamma}
}
\]

だけになる。

実コードでは interaction tensor の非零 pair を monomial vector に圧縮し、

```python
mon = z[pi] * z[pj]
return M @ mon
```

とする。

これは transition の全コードであり、`transition()` 内には

- `if`
- `for`
- RK stage
- `copy`
- `+=`, `-=`
- `log`, `arg`

が存在しないことを AST 監査で確認した。

---

# 3. RN coupled spin-sector 実装

対象：`strict_harmonic_rn_tensor_v1.py`

## 3.1 初期状態

匿名波チャネルを \(\psi_0,\psi_1\) とし、物理的な spin-2 / spin-1 の名称は transition では使用しない。

背景質量は自然単位

\[
M=1
\]

とする。

ケース依存量 \(Q\) は外部 transition parameter ではなく状態へ入れる。その多項式係数を一回積和で扱えるよう、

\[
q_0=1,
\quad q_1=Q,
\quad q_2=Q^2,
\quad q_3=Q^3,
\quad q_4=Q^4
\]

を状態チャネルとして初期化し、すべて同じ interaction array で恒等保持する。

初期化後、transition は `q` を関数引数として受け取らない。

## 3.2 根拠方程式

Reissner--Nordström 背景の odd-parity gravito-electromagnetic coupled master system を

\[
\partial_t^2
\begin{pmatrix}G\\E\end{pmatrix}
=
\partial_{r_*}^2
\begin{pmatrix}G\\E\end{pmatrix}
-
\begin{pmatrix}
V_{GG}&V_{GE}\\
V_{GE}&V_{EE}
\end{pmatrix}
\begin{pmatrix}G\\E\end{pmatrix}
\]

とする。

\[
f(r)=1-\frac{2}{r}+\frac{Q^2}{r^2},
\qquad
\partial_{r_*}=f\partial_r
\]

より

\[
\partial_{r_*}^2\psi
=f^2\partial_r^2\psi
+ff'\partial_r\psi.
\]

\(\ell=2\) について実装した potential は

\[
V_{GG}
=f\left(\frac{\ell(\ell+1)}{r^2}
-\frac6{r^3}+\frac{4Q^2}{r^4}\right),
\]

\[
V_{EE}
=f\left(\frac{\ell(\ell+1)}{r^2}
+\frac{4Q^2}{r^4}\right),
\]

\[
V_{GE}
=f\frac{2Q\sqrt{(\ell-1)(\ell+2)}}{r^3}.
\]

これらはすべて \(Q^0\) から \(Q^4\) の有限多項式になるため、固定係数配列と状態 \(q_0,\ldots,q_4\) の積和だけへ展開できる。

## 3.3 「前時刻」も別状態にしない

通常の二階差分

\[
\psi^{n+1}
=2\psi^n-\psi^{n-1}+\Delta t^2\mathcal L\psi^n
\]

に必要な \(\psi^{n-1}\) は外部配列 `prev` として保持しない。

有限履歴 \(K=4\) 全体が倍音係数 \(C_m\) の中に存在し、履歴 shift と newest-sample generation をまとめて Fourier 基底の一つの interaction tensor にコンパイルしている。

従って loop 境界を越える物理状態は `z` のみである。

## 3.4 実行ループ

```python
for n in range(steps):
    old = z
    z = transition(old, M, pi, pj)
```

`old` は同一 step 内の read-only 参照であり、次 step の物理状態として別に保持されない。

実際の transition は

```python
def transition(z, M, pi, pj):
    mon = z[pi] * z[pj]
    return M @ mon
```

だけである。

## 3.5 読出し

DFT decode、波形、flux、spin-2 / spin-1 の意味付けは観測専用で、`z=transition(...)` の右辺へ戻らない。

\[
Z_n\rightarrow \text{observation}
\]

のみである。

---

# 4. 同一 interaction array での電荷スイープ

最終実験条件：

\[
Q/M=0,\;0.3,\;0.6,\;0.9,\;0.99,
\]

\[
N_r=61,
\quad K=4,
\quad \Delta t=0.12,
\quad 700\ \text{steps}.
\]

全ケースで同じ interaction array を使用し、変更したのは初期状態中の \(Q\) だけである。

pure channel-0 初期状態から channel-1 へ出た boundary-flux 比：

| Q/M | other-channel fraction |
|---:|---:|
| 0.00 | 0 |
| 0.30 | 0.0016139064 |
| 0.60 | 0.0062901234 |
| 0.90 | 0.0137090482 |
| 0.99 | 0.0164396584 |

pure channel-1 初期状態から channel-0 への比：

| Q/M | other-channel fraction |
|---:|---:|
| 0.00 | 0 |
| 0.30 | 0.0014597175 |
| 0.60 | 0.0057071048 |
| 0.90 | 0.0125027843 |
| 0.99 | 0.0150217465 |

\(Q=0\) では交差変換は数値上完全に0であり、\(Q\neq0\) で連続的に立ち上がる。

状態内の \(Q^p\) consistency error は最大でも約

\[
1.1\times10^{-16}
\]

である。

---

# 5. spin-2 波の放射・吸収バランス再実験

同じ strict harmonic interaction を \(Q=0\) のまま使用し、初期波束位置だけを初期状態として変化させた。

\[
x_0\approx2.97890625
\]

で

\[
\frac{F_{\rm out}-F_{\rm in}}
{F_{\rm out}+F_{\rm in}}
\approx-6.44\times10^{-5}
\]

を得た。

これは旧 RWZ leapfrog の結果を転記したものではなく、今回の倍音主状態・一回積和写像で再探索した値である。

---

# 6. harmonic map 自体の独立照合

同じ離散 RN recurrence を sample-history 配列で外部参照計算し、strict harmonic map と120 step 比較した。

全履歴・全チャネルの最大差：

\[
\boxed{3.0654\times10^{-17}}
\]

したがって Fourier 基底への interaction-array compilation 自体が物理更新を変更していないことを確認した。

この sample-space 計算は **validation only** であり、strict experiment の transition には使用しない。

---

# 7. full nonlinear Kastor--Traschen sector の再構築

対象：`strict_harmonic_kt_tensor_v1.py`

旧コードの

\[
t\mapsto a(t)=e^{Ht}
\]

という外部 exact-state injection を廃止した。

新状態は匿名 scalar channels

\[
(1,\ a,\ h,\ m_1,\ m_2,\ x_1,\ x_2)
\]

の有限履歴を倍音係数として保持する。ここで

\[
h=e^{H\Delta t}
\]

自身も状態であり、transition は外部 \(H\) を読まない。

一回写像は sample-space で

\[
a' = ah,
\]

\[
1'=1\cdot1,
\quad
h'=h\cdot1,
\quad
m_A'=m_A\cdot1,
\quad
x_A'=x_A\cdot1
\]

だけから構成し、これを Fourier 基底へ固定二次 interaction としてコンパイルする。

実コード：

```python
def transition(z,M,pi,pj):
    return M @ (z[pi] * z[pj])
```

生成後の状態からのみ

\[
W=a+\frac{m_1}{r_1}+\frac{m_2}{r_2}
\]

を観測し、Kastor--Traschen exact geometry と比較する。

32 step 後まで、状態生成された scale factor と外部 exact benchmark の差は

\[
0\ \text{to}\ 1.1\times10^{-16}
\]

である。

また

\[
W^2
=a^2+2au_a+2au_b+u_a^2+2u_au_b+u_b^2
\]

から

\[
aa=u_a^2,
\quad ab=2u_au_b,
\quad bb=u_b^2
\]

を readout し、再構成誤差は全チェックポイントで0であった。

重要：`aa`, `ab`, `bb` は persistent state でも別 interaction でもない。

---

# 8. spin-1/2 型 / spin-2 型読出しの再構築

対象：`strict_harmonic_spin_tensor_v1.py`

匿名複素状態 \(V,a,b\) と固定回転状態 \(R_{1/2},R_1\) をすべて倍音 persistent state に入れ、

\[
V'=V R_{1/2},
\qquad
a'=aR_1,
\qquad
b'=bR_1
\]

を同一固定二次 interaction で更新する。

\[
R_{1/2}=e^{i\pi/N},
\qquad
R_1=e^{i2\pi/N},
\qquad
R_{1/2}^2=R_1.
\]

readout では

\[
V(N)=-V(0),
\qquad
V(2N)=V(0)
\]

および

\[
aa=a^2,
\quad ab=2ab,
\quad bb=b^2
\]

の weight-2 covariance を確認する。

\(N=64\) の結果：

- 2\(\pi\) sign-flip error: \(1.16\times10^{-14}\)
- 4\(\pi\) return error: \(2.40\times10^{-14}\)
- weight-2 one-step covariance error: \(4.48\times10^{-16}\)
- \(R_{1/2}^2-R_1\): 0

spin の名称は readout interpretation であり、transition 内に spin 分岐は存在しない。

---

# 9. R1--R14 監査

| Rule | 新実装 |
|---|---|
| R1 状態閉包 | ケース依存 \(Q\), その有限履歴、必要な生成状態を倍音状態へ内包 |
| R2 一回写像 | `M @ (z[pi]*z[pj])` の一回のみ |
| R3 恒等写像 | unit, q-powers, source strengths 等も同じ map の積で保持 |
| R4 乗法形式 | transition 内に RK / `state+increment` / log-exp 往復なし |
| R5 同期更新 | 全 monomial は旧 `z` のみから生成 |
| R6 持越し禁止 | loop 境界を越える物理状態は `z` のみ |
| R7 読出し非帰還 | decode / flux / geometry / spin readout は transition に戻さない |
| R8 無名性 | 実行 map は匿名 channel index。grav/EM/spin 名称は readout のみ |
| R9 固定係数 | grid/discretization と標準方程式係数のみ。ケース依存 Q は状態 |
| R10 自然単位 | \(G=M=c=1\) |
| R11 原則と数値一致を分離 | 原則監査を先に実施し、その後で数値照合 |
| R12 厳密一致 | 旧不適合コードを superseded と明記、新コード・data・auditを同時保存 |
| R13 再現性 | code, interaction arrays, raw JSON, logs, hashes を保存 |
| R14 無断例外禁止 | 今回新しい例外は導入していない |

---

# 10. 重要な限界

この再構築によって確認したのは、

1. finite-history harmonic state を唯一の persistent state とできること、
2. RN coupled spin-sector の離散方程式を一つの固定二次積和 interaction へコンパイルできること、
3. Q=0 の分離、Q≠0 のチャネル変換、放射/吸収バランス点を同じ map で再現できること、
4. KT exact sector の時間状態を外部時刻なしの乗算状態生成へ置換できること、
5. aa/ab/bb および spin-like 量を readout として保持できること、

である。

一般の asymptotically-flat full nonlinear charged-binary merger の3D numerical relativity を、この計算資源で実行したことを意味しない。その検証は HPC が必要であり、本稿では主張しない。
