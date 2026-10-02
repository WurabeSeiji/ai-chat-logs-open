# 万能相互作用 Einstein–Maxwell 乗法更新 完全仕様 v3.2
## — 背景時空も粒子も先に置かず、厳密な Einstein–Maxwell 方程式を状態依存乗法写像 S(Ψ) で発展させ、系の唯一の時計と a,b の内部時計を読み出す —

**著者:** 木原範昭 (Noriaki Kihara)  
**仕様整理:** ChatGPT（v3.1 まで）、Claude（v3.2 改定）  
**版:** v3.2  
**日付:** 2026-10-01  
**置換対象:** `universal_interaction_einstein_maxwell_multiplicative_spec_ja_v3.1_20260930.md`  
**重要:** 本仕様はコード実装ではない。実装開始前の正本仕様である。

---

# 0A. v3.2 改定（2026-10-01、木原範昭の指示による）

## 0A.1 本実験の目的（v3.1 までの誤った問いを置き換える）

v3.1 までは本実験の問いを「既存の加法更新と状態依存乗法更新が一致するか」としていた。これは誤りである。§23 が示すとおり、任意の加法更新は affine lift によって乗法形式に書けるので、その一致は恒等的に成り立ち、実験の問いにならない。

本実験の目的は次である。

1. 本研究の理論系はスケール不変を原則とする。したがって、天体スケールで成り立つ厳密な方程式（一般相対論＋Maxwell 方程式）が、原子スケール（Planck 長）でも成り立つという仮説に立つ。
2. その仮説のもとで、これまでの理論系で未解決だった**動力学**を、厳密な Einstein–Maxwell 方程式から導く。
3. このとき次の三原則を守る。
   - **背景時空を先に設定しない。** 時空は、重力による歪みが生じて初めて読み出せるようになる。座標格子はゲージのラベルであり、背景ではない。
   - **粒子すら仮定しない。** 限られた状態の変化だけから、時空と粒子様の動力学が読み出せるかを見る。外から素粒子状態を与えてはならない。一方、結果として素粒子状態に類似したものが読み出された場合に、そう解釈することは禁止しない。
   - **状態は波であり、読み出せるのは同一位相内の現在値だけである。** 読出し・演算・書込みは同時に行われる。これが積和（乗法更新）だけを前提にする理由である。

乗法更新を要求する理由は、位相と振幅の解釈が極めて容易になるからである。複素数の積では、振幅は掛け算で、位相は足し算で、一回の積の中で同時に更新される。状態を積で進めれば、各刻みの位相の進みと振幅の変化をそのまま読み出せる。

加法参照系（§6）は、実装が意図どおりの計算をしているかを確かめる**実装検算**としてのみ残す。実験の問いではない。

## 0A.2 方程式は変更しない

Einstein–Maxwell 方程式は変更しない。背景の空間も時間もない系では「点」に意味がない。したがって、点の発散を避けるために導入された仮定（regular Kerr–Newman ゲート、excision の必須化、地平面での境界条件）は外す。この系で Einstein–Maxwell 方程式が成り立たないと、実装者が勝手に仮定してはならない。

## 0A.3 単位と電荷の規約

\[
\boxed{G=q_0=c=1}
\]

とする。\(q_0\) は電荷素量の 1/3 である。電子・陽子の電荷は \(Q=\mp3\) である。v3.1 までの文書に残っていた \(G=M=c=1\)（質量を単位とする規約）は誤りである。

電荷の規約は、物理ではなく選択の問題である。本仕様は、使用する厳密解（Reissner–Nordström・Kerr–Newman 計量、\(f=1-2M/r+Q^2/r^2\)）と標準理論の式をそのまま用いるため、Einstein–Maxwell の標準文献が用いる **Gauss 幾何単位**（\(G=c=1,\ 4\pi\epsilon_0=1\)）を採用する。微細構造定数は標準的な定義

\[
\alpha=\frac{e^2}{\hbar c}
\]

で扱う。Heaviside–Lorentz 規約で書かれた先行論文（反復交換散乱における有限位数共鳴、透過振幅 \(1-R=e=\sqrt{4\pi\alpha}\)）と比較する場合に限り、\(e_{\rm HL}=\sqrt{4\pi}\,e_{\rm Gauss}\) と換算する。

## 0A.4 系の唯一の時計 δτ と、1 刻み内の 4 段

\(\delta\tau\) は数値計算の刻み幅ではない。系の内部の唯一の時計の 1 刻みである。放射があれば系の時計は変わるので、\(\delta\tau\) も変わってよい。

本仕様では、系の時計の 1 周期がちょうど

\[
\boxed{N_{\rm tick}=31}
\]

刻みになるように \(\delta\tau\) を置く。1 刻みの中は 4 段（§5）であり、段で数えると \(31\times4=124\) で 1 周期が閉じる。

\(\delta\tau\) を細かくして連続極限への収束を見ることはしない。与えた \(\delta\tau\) で得た結果がそのまま結果であり、連続解の近似ではない。

## 0A.5 注意（将来の誤読防止）

> **1 刻み内を 4 段にしていること、および系の時計の周期を 31 刻みとなるよう \(\delta\tau\) を選んでいることは、微細構造定数 \(\alpha\) の値に対応する有限位数根 \(m/n=23/124\)（\(124=31\times4\)、反対称固有値 \(\lambda_a^{31}=i\)）から予測値として与えたものである。**
>
> **本研究は、最小の公理系が自己無撞着になることを条件と考えている。そのため循環のように見えるが、これは隠れたパラメータの導入ではない。任意に与えたパラメータが自己無撞着になるには特定の値でなければならないことを示そうとするものである。**
>
> 4 段は、スピン 1/2 の位相の一巡（\(1,i,-1,-i\)。\(2\pi\) で \(-1\)、\(4\pi\) で \(+1\)）と構造上対応する。これはスピン 1/2 の条件にたまたま合っていた可能性があり、スピン 2 など他のスピンでは段数が異なる可能性がある。
>
> 同じ理由で、a,b の初期スピン \(J=\hbar/2=9/(2\alpha)\)（Gauss 規約、\(Q=3\)）も \(\alpha\) から与える。標準理論は \(\alpha\) や質量比を外から与える任意パラメータとして扱うが、本系ではそれらの値が自己無撞着の条件になる。なぜ物理がその値を選ぶのかは不明である。その値を選べば系が閉じる、ということだけを確かめる。導出は必ず循環するので、「入力だから導出ではない」という批判は本研究の枠組みに当たらない。

## 0A.6 予想は成功条件にしない

ロックの有無、ロック比、閉曲面の選び方への依存、\(\delta\tau\) の変化、\(\alpha\)・\(\hbar\) の読み出し値などについての予想は、外れても系の動作に影響しない。これらは読み出しとして保存し、結果が出てから解釈する。実験の成功条件は、方程式を変えずに適用していること、隠れた仮定がないこと、規則どおりに状態を更新していること、同じ条件で同じ結果が再現されること、全読み出しを保存していることに限る。

---

# 0. 改定理由（v3.1 までの記録）

v2.0 は、Einstein–Maxwell–GH の連続 RHS を `K(Ψ)` へ分解し、CF4 と matrix exponential で

\[
Ψ_{n+1}=e^B e^A Ψ_n
\]

と更新する方針を採った。

しかし今回の研究目的は、既存の数値相対論がすでに実行している

\[
Ψ_{n+1}=Ψ_n+\Delta Ψ_n
\]

型の差分更新を、物理方程式・空間差分・境界処理・時間積分を変えずに、

\[
\boxed{Ψ_{n+1}=S(Ψ_n)Ψ_n}
\]

という**乗法的な万能相互作用**へ写像することである。

したがって v3.0 では、

1. Einstein–Maxwell の数値 RHS は既存の加法数値系と同一とする。
2. 正式な状態更新は乗法系 `S(Ψ)Ψ` のみとする。
3. 同一条件で加法参照系を並走させ、両者の一致を逐step検算する。
4. `T`, `δτ`, `Dτ`, `κτ`, `ONE` をすべて状態に含める。
5. 加法参照系のために乗法系の状態ルールを緩めない。
6. v2.0 で唯一の時間更新法としていた CF4/expmv は本仕様では採用しない。

と改定する。

## 0.1 v3.1 監査修正

v3.1 は v3.0 の中心相互作用

\[
\boxed{
Psi_{n+1}^{(\times)}=S(Psi_n^{(\times)})Psi_n^{(\times)}
}
\]

および加法検算系との二経路比較を変更しない。今回の修正は、前チャットで確定した許容範囲に合わせ、実装時に誤読され得る契約だけを明確化する。

1. 正式物理状態 `Psi^(x)` は初期化後、stageを含めて `Psi += DeltaPsi` / `u += Deltau` で更新してはならない。
2. 加法参照系 `Psi^(+)` は独立した検算専用状態として保持してよく、加法更新も許される。ただし正式系へ一切帰還させない。
3. `F_h` は入力された正式状態を非破壊で評価し、RHS を別出力として返す。
4. case 情報の利用自体は許されるが、case に応じて物理方程式、RHS、interaction、正式 time-evolution route を切り替えてはならない。
5. solver workspace、ghost buffer、FFT plan、preconditioner、allocation/cache 等の計算補助情報は許される。ただし hidden physical/dynamical state として次状態を左右してはならない。
6. `S[q,d]Psi` は dense 行列を明示生成しなくてもよい。sparse/structured operator application は許されるが、正式状態を in-place の `+=` で作る shortcut は禁止する。

---

# 1. 状態更新の形式

本実験の目的は §0A.1 に定める。本節は、その目的のための状態更新の形式を定める。

物理・数値系が与える一step写像を

\[
\boxed{
\mathcal G_h: Ψ_n \mapsto Ψ_{n+1}
}
\]

とする。

通常の数値実装では、動的場 `u` に対して

\[
u_{n+1}=u_n+\Delta u_n
\]

と加法的に更新する。

万能相互作用では同じ \(\Delta u_n\) を用いるが、状態更新そのものは

\[
\boxed{
Ψ_{n+1}=S_n Ψ_n,
\qquad
S_n=S(Ψ_n)
}
\]

という行列作用だけで行う。

実装検算として、同じ初期状態、同じ空間離散化、同じ境界処理、同じ時間刻み、同じ時間積分公式に対して

\[
Ψ_{n+1}^{(\times)}=Ψ_{n+1}^{(+)}
\]

が数値丸め誤差の範囲で毎step成立することを確認する。これは実装が意図どおりの計算をしているかの確認であり、実験の問いではない（§0A.1、§23）。

---

# 2. 絶対ルール

1. **正式系は乗法系である。** 次状態は必ず `S(Ψ)Ψ` で生成する。
2. 正式系で `Ψ += ΔΨ`、`u += Δu` を次状態更新として使用してはならない。
3. 加法系は**検算専用**であり、正式系の次状態生成へ帰還させない。
4. 乗法系と加法系は同じ Einstein–Maxwell RHS evaluator を使用する。
5. 乗法系と加法系で物理方程式、空間差分、境界条件、ゲージ、constraint処理、時間積分係数を変えてはならない。
6. `case_id` 等のケース情報によって、物理方程式、RHS、interaction、正式 time-evolution route を切り替えることは禁止する。初期条件、topology、validation種別、出力、readout、停止条件など、物理法則そのものを変更しない用途へのケース情報利用は許される。
7. 質量・電荷・スピン等の初期ラベルを時間発展係数へ直接再注入しない。
8. `T_τ`, `δτ`, `D_τ`, `κ_τ`, `ONE` は状態成分である。
9. 外部 `dt`、外部現在時刻、外部step counterを物理更新へ使用してはならない。
10. 数値解法に必要な固定係数・差分係数・境界係数・時間積分係数は、仕様で固定された恒等状態として保持する。
11. readoutは一方向であり dynamics へ帰還しない。
12. 加法系との一致を得るために乗法系の定義をstep途中で変更してはならない。
13. `F_h` は入力状態を read-only として扱い、正式状態・stage状態を in-place 更新してはならない。
14. `S[q,d]Psi` の実装は operator application として行う。fresh output bufferへの structured/sparse matvec は許されるが、正式状態に対する `+=` shortcut は禁止する。
15. 単位は \(G=q_0=c=1\)、電荷は Gauss 幾何単位とする（§0A.3）。質量を単位とする \(G=M=c=1\) を使ってはならない。
16. Einstein–Maxwell 方程式を変更しない。点の発散を理由にした正則性ゲート・excision・地平面境界条件を導入しない（§0A.2）。
17. `δτ` は恒等状態ではなく、系の時計の 1 周期が \(N_{\rm tick}=31\) 刻みになるよう各 step の開始時に状態から定める（§3.2）。

---

# 3. 完全状態ベクトル Ψ

状態を

\[
\boxed{
Ψ=
(Ψ_{unit},Ψ_{time},Ψ_{phys0},Ψ_{numerics},Ψ_{grid},Ψ_{field})^T
}
\]

とする。

## 3.1 unit state

\[
\boxed{ONE=1}
\]

`ONE` は乗法的 affine lift の基準成分である。

## 3.2 time state

時間関連状態として少なくとも

\[
\boxed{
T_τ,
\delta τ,
D_τ,
\kappa_τ
}
\]

を持つ。

定義：

\[
\boxed{T_τ=e^{\kappa_τ τ}}
\]

\[
\boxed{D_τ=e^{\kappa_τ\delta τ}}
\]

したがって監査式は

\[
\boxed{
\delta τ=\frac{\log D_τ}{\kappa_τ}
}
\]

である。

v3.2 では `κ_τ≠0`, `T_τ>0`, `D_τ>0`, `δτ>0` を必須とする。

`δτ` は状態から読む。外部引数 `dt` を dynamics に渡してはならない。

### 3.2.1 δτ は系の唯一の時計の刻み（v3.2）

v3.1 では `δτ,Dτ` を恒等状態としていた。これは誤りである。`δτ` を固定することは、放射によって系が変化しない（定常である）ことを先に仮定するのと同じだからである。

v3.2 では、`δτ` を系の内部の唯一の時計の 1 刻みとし、各 step の開始時に状態から次のように定める。

\[
\boxed{
\delta\tau_n=\frac{P_{\rm sys}(Ψ_n)}{N_{\rm tick}},
\qquad
N_{\rm tick}=31
}
\]

\[
\boxed{
D_{τ,n}=e^{\kappa_τ\delta\tau_n}
}
\]

ここで \(P_{\rm sys}=2\pi/\Omega_{\rm sys}\) は系の時計（a,b の相対回転）の 1 周期である。束縛した二体系の時計は一つであり、自由度は相対回転の角振動数 \(\Omega_{\rm sys}\) 一個である（各天体の固有時計は、その時計に対する固定比＝赤方偏移で決まる）。

- \(N_{\rm tick}=31\) と 1 刻み内の 4 段の由来と位置づけは §0A.5 の注意書きのとおりである。
- `δτ,D_τ` は状態成分であり、`S(Ψ)` の同じ一回の作用で次 step の値が生成される。
- `δτ` は時間状態そのものであり、readout ではない。したがって §2 の「readout を dynamics へ帰還させない」規則には当たらない。
- \(\Omega_{\rm sys}(Ψ)\) の推定式は現在状態だけから評価し、履歴を使わない。採用する推定式は §24 の実装開始前監査で一意に確定する。
- 与えた \(N_{\rm tick}\) で得た結果がそのまま結果である。`δτ` を細かくして連続極限への収束を見ることはしない。

`κ_τ` は恒等状態である。

時間更新は

\[
\boxed{T_{τ,n+1}=D_{τ,n} T_{τ,n}}
\]

とする。したがって読み出し時刻は刻みの積算

\[
τ_n=\frac{\log T_{τ,n}}{\kappa_τ}=τ_0+\sum_{k<n}\delta\tau_k
\]

であり、step 番号の定数倍ではない。

## 3.3 dynamic field state

物理発展する場をまとめて

\[
\boxed{u}
\]

と書く。

GH–Maxwell 表現を採る場合は従来どおり

- \(g_{ab}\)
- \(\Pi_{ab}\)
- \(\Phi_{iab}\)
- \(H_a\)
- \(\Theta_a\)
- \(E^i\)
- \(B^i\)

を含む。

v3.1 の乗法化は、この物理 field representation 自体を変更しない。

---

# 4. 数値 Einstein–Maxwell RHS の唯一性

空間離散化後の Einstein–Maxwell 系を

\[
\boxed{
\frac{du}{dτ}=F_h(Ψ)
}
\]

と定義する。

ここで `h` は空間離散化全体を表し、

- spatial derivative
- domain/interface coupling
- boundary treatment
- gauge terms
- constraint damping / cleaning（採用する場合）
- matter backreaction

を含む。

重要なのは、**乗法系と加法系の両方が同じ `F_h` 関数を呼ぶこと**である。

`F_h^(×)` と `F_h^(+)` の二種類を作ってはならない。

実装監査では、RHS evaluator のコードハッシュが両経路で同一であることを記録する。

## 4.1 RHS evaluator の非破壊契約

`F_h` は次の契約を満たす。

\[
\boxed{
F_h:\;Psi_{in}^{RO}\longmapsto k
}
\]

ここで `RO` は read-only を意味する。`F_h` の評価中に、入力された正式状態またはstage状態を in-place に変更してはならない。

禁止例：

```text
Psi[field] += correction
Psi[:] = filtered_Psi
apply_boundary_in_place(Psi)
```

許容例：

- ghost / boundary 用の別workspace生成
- derivative / flux / constraint の一時配列
- filter結果をRHS側の別配列へ反映
- Newton / Krylov / FFT等のsolver workspace
- FFT plan、preconditioner、allocation cache等の計算補助情報

ただし、これらの補助情報が `Psi` に存在しない hidden physical/dynamical state として、同じ `Psi`・同じ恒等設定に対する次状態を履歴依存で変更してはならない。キャッシュは、状態・固定設定から定まる演算の再利用・高速化に限る。

---

# 5. 時間積分公式

v3.2 baseline は classical explicit RK4（4 段）とする。

**4 段の位置づけ（v3.2）**：4 段は `δτ` より細かい途中の時刻を作るものではなく、1 刻みの中で位相の一巡を表すものと読む。1 刻みの発展 \(e^{z}\) を \(z^4\) の項まで合わせると、回転 \(z=i\theta\) では各項の位相が \(1,i,-1,-i\) と 4 分の 1 周ずつ進み、4 項で一巡する（\(i^4=1\)）。これはスピン 1/2 の二重被覆（\(2\pi\) で \(-1\)、\(4\pi\) で \(+1\)）と構造上対応する。段数の由来と、スピン 2 など他のスピンで段数が異なる可能性については §0A.5 の注意書きのとおりである。

Butcher tableau：

\[
\begin{array}{c|cccc}
0   & 0   & 0   & 0 & 0\\
1/2 & 1/2 & 0   & 0 & 0\\
1/2 & 0   & 1/2 & 0 & 0\\
1   & 0   & 0   & 1 & 0\\
\hline
    & 1/6 & 1/3 & 1/3 & 1/6
\end{array}
\]

係数

\[
a_{rs},\quad b_r,\quad c_r
\]

は `Ψ_numerics` の恒等状態として保持する。

`time_integrator_id = RK4_CLASSICAL` も恒等状態として保持する。

別の時間積分法を使用する場合は、別runで state を変更し、加法系と乗法系の両方を同時に同じtableauへ変更する。

---

# 6. 加法参照系 — 検算専用

加法系は数値相対論で通常行う更新を、そのまま参照計算として実行する。

step開始状態を

\[
u_n^{(+)}=u_n
\]

とする。

時間 stage は状態から生成する。

\[
T_r=D_τ^{c_r}T_{τ,n}
\]

または等価に

\[
T_r=e^{\kappa_τc_r\delta τ}T_{τ,n}
\]

とする。

RK4 stage：

\[
k_1^{(+)}=F_h(Ψ_n^{(+)})
\]

\[
u_2^{(+)}=u_n+\frac{\delta τ}{2}k_1^{(+)}
\]

\[
k_2^{(+)}=F_h(Ψ_2^{(+)})
\]

\[
u_3^{(+)}=u_n+\frac{\delta τ}{2}k_2^{(+)}
\]

\[
k_3^{(+)}=F_h(Ψ_3^{(+)})
\]

\[
u_4^{(+)}=u_n+\delta τ k_3^{(+)}
\]

\[
k_4^{(+)}=F_h(Ψ_4^{(+)})
\]

final increment：

\[
\boxed{
\Delta u_{RK4}^{(+)}
=\frac{\delta τ}{6}
(k_1^{(+)}+2k_2^{(+)}+2k_3^{(+)}+k_4^{(+)})
}
\]

final field：

\[
\boxed{
u_{n+1}^{(+)}=u_n+\Delta u_{RK4}^{(+)}
}
\]

時間：

\[
\boxed{T_{τ,n+1}^{(+)}=D_{τ,n}T_{τ,n}}
\]

`δτ,Dτ` は §3.2.1 と同じ規則で次 step の値を定める。`ONE,κτ` とその他恒等状態は変更しない。

**この加法系は検算結果を作るためだけに存在し、正式な状態列を生成しない。**

---

# 7. 乗法正式系 — stage更新も全て行列作用

乗法系では、RK stage state を作るときも `u += ...` を実行しない。

任意の field increment vector \(q\) と time factor \(d\) に対し、乗法 affine operatorを

\[
\boxed{
S[q,d]
}
\]

と定義する。

state ordering を

\[
Ψ=(ONE,T_τ,δτ,D_τ,κ_τ,I_1,\ldots,I_m,u_1,\ldots,u_N)^T
\]

とすると、`S[q,d]` は次だけを非恒等にする。

### ONE

\[
S_{ONE,ONE}=1
\]

### Tτ

\[
\boxed{S_{T_τ,T_τ}=d}
\]

### 恒等状態

\[
S_{I_j,I_j}=1
\]

### dynamic fields

\[
\boxed{S_{u_A,u_A}=1}
\]

\[
\boxed{S_{u_A,ONE}=q_A}
\]

その他は0。

よって

\[
\boxed{
S[q,d]Ψ
=
(1,dT_τ,δτ,D_τ,κ_τ,I_1,\ldots,I_m,u+q)^T
}
\]

となる。

これは**加算更新をコードで実行したものではなく、拡張状態に対する一回の行列乗算**である。

---

# 8. RK4 stage の乗法生成

step開始状態を

\[
Ψ_n^{(×)}=Ψ_n
\]

とする。

## stage 1

\[
\boxed{k_1^{(×)}=F_h(Ψ_n^{(×)})}
\]

stage 2 operator：

\[
q_2=\frac{\delta τ}{2}k_1^{(×)}
\]

\[
d_2=D_τ^{1/2}
\]

\[
\boxed{
Ψ_2^{(×)}=S[q_2,d_2]Ψ_n^{(×)}
}
\]

\[
\boxed{k_2^{(×)}=F_h(Ψ_2^{(×)})}
\]

## stage 3

\[
q_3=\frac{\delta τ}{2}k_2^{(×)}
\]

\[
d_3=D_τ^{1/2}
\]

\[
\boxed{
Ψ_3^{(×)}=S[q_3,d_3]Ψ_n^{(×)}
}
\]

\[
\boxed{k_3^{(×)}=F_h(Ψ_3^{(×)})}
\]

## stage 4

\[
q_4=\delta τ k_3^{(×)}
\]

\[
d_4=D_τ
\]

\[
\boxed{
Ψ_4^{(×)}=S[q_4,d_4]Ψ_n^{(×)}
}
\]

\[
\boxed{k_4^{(×)}=F_h(Ψ_4^{(×)})}
\]

---

# 9. final universal interaction S(Ψn)

RK4 increment を

\[
\boxed{
\Delta u_n^{(×)}
=\frac{\delta τ}{6}
(k_1^{(×)}+2k_2^{(×)}+2k_3^{(×)}+k_4^{(×)})
}
\]

とする。

正式な万能相互作用配列は

\[
\boxed{
S(Ψ_n)=S[\Delta u_n^{(×)},D_τ]
}
\]

である。

v3.2 では、同じ `S(Ψ_n)` の作用で時間状態の次の値も生成する。すなわち `S(Ψ_n)` は §7 の非恒等成分に加えて、

\[
S_{\delta\tau,\delta\tau}=0,\qquad S_{\delta\tau,ONE}=\delta\tau_{n+1},
\qquad
S_{D_τ,D_τ}=0,\qquad S_{D_τ,ONE}=e^{\kappa_τ\delta\tau_{n+1}}
\]

を持つ。\(\delta\tau_{n+1}\) は §3.2.1 の規則で状態から定める値であり、外部から与えない。

従って正式な次状態は

\[
\boxed{
Ψ_{n+1}^{(×)}=S(Ψ_n)Ψ_n^{(×)}
}
\]

のみで生成する。

この式が本仕様の**唯一の正式 time evolution route**である。

---

# 10. 「乗法系」と呼ぶための禁止事項

次はすべて禁止する。

1. `u_new = u + Delta_u` を正式系で実行してから `S` にコピーする。
2. `S` をログ用にだけ生成し、実際の更新は加法で行う。
3. 加法系の最終 `u_new` を乗法系へコピーする。
4. 加法系の stage state を乗法系へコピーする。
5. 加法系の `k1..k4` を乗法系へコピーする。
6. 乗法系の `k1..k4` を加法系へコピーする。
7. 二経路で異なる RHS、異なる境界条件、異なる差分係数を用いる。
8. 検算誤差を減らすため片側だけ補正する。
9. `Tτ` だけ外部時刻で更新する。
10. `δτ` を外部 `dt` から与える。

両経路は**同じ関数定義を独立に評価する**。

## 10.1 operator application の実装契約

正式系での `S[q,d]Psi` は、数学上の operator application を実装する。dense matrix の明示生成は要求しない。

例えば structured operator kernel により fresh output buffer `Psi_out` へ

\[
(Psi_{out})_A=\sum_B S_{AB}(Psi_{in})_B
\]

を直接評価してよい。`S` が疎であることを利用した component-wise 実装も許される。

ただし正式状態に対し

```text
Psi_next = Psi_current.copy()
Psi_next[field] += q
```

のように、状態加算を time-evolution route として実行することは禁止する。fresh output へ operator の各成分を評価することと、正式状態を `+=` で更新することを区別する。

---

# 11. 二経路の独立性

加法系と乗法系は同じ source code の `F_h` を呼ぶが、stage state と stage RHS の保存領域は分離する。

必要な一時量：

加法側：

\[
Ψ_r^{(+)},\quad k_r^{(+)}
\]

乗法側：

\[
Ψ_r^{(×)},\quad k_r^{(×)}
\]

である。

双方を同じメモリ参照にaliasしてはならない。

比較前に片方を他方へ上書きしてはならない。

---

# 12. stage-by-stage equivalence audit

最終値だけでなく各stageを比較する。

\[
\boxed{
E_{stage,r}
=\|u_r^{(×)}-u_r^{(+)}\|_\infty
}
\]

\[
\boxed{
E_{rhs,r}
=\|k_r^{(×)}-k_r^{(+)}\|_\infty
}
\]

時間状態も

\[
\boxed{
E_{T,r}=|T_r^{(×)}-T_r^{(+)}|
}
\]

を保存する。

最終step：

\[
\boxed{
E_{step}
=\|Ψ_{n+1}^{(×)}-Ψ_{n+1}^{(+)}\|_\infty
}
\]

relative error：

\[
\boxed{
E_{rel}
=\frac{\|Ψ_{n+1}^{(×)}-Ψ_{n+1}^{(+)}\|_\infty}
{\max(1,\|Ψ_{n+1}^{(+)}\|_\infty)}
}
\]

基準：

\[
\boxed{E_{rel}<\epsilon_{equiv}}
\]

を必須とする。

不一致時に一方を補正せず、runを `MULTIPLICATIVE_ADDITIVE_MISMATCH` として停止する。

---

# 13. identity-state audit

両経路で毎step、

\[
ONE=1
\]

\[
\kappa_{τ,n+1}=\kappa_{τ,n}
\]

\[
N_{{\rm tick},n+1}=N_{{\rm tick},n}=31
\]

および全固定数値設定、初期ラベル、grid設定が不変であることを監査する。

v3.2 では `δτ,D_τ` は恒等状態ではない（§3.2.1）。代わりに、各 step で \(\delta\tau_{n+1}\) が §3.2.1 の規則から再計算した値と一致していること、および乗法系と加法系で同じ値になっていることを監査する。

乗法系では該当行が単位行でなければFAIL。

加法系では該当成分のincrementが0でなければFAIL。

---

# 14. time-state consistency audit

毎step開始時と終了時に

\[
\boxed{
C_{δτ}
=\delta τ-\frac{\log D_τ}{\kappa_τ}
}
\]

を監査する。

\[
|C_{δτ}|<\epsilon_{time-state}
\]

を必須とする。

また、step番号 `n` は物理時刻の正本ではない。

物理時刻 readout は

\[
\boxed{
τ=\frac{\log T_τ}{\kappa_τ}
}
\]

だけから行う。これは刻み \(\delta\tau_k\) の積算であり、系の時計で数えた時刻である。a,b それぞれの固有時は、読出し仕様で lapse（赤方偏移）を用いて別に読み出す。

---

# 15. Einstein–Maxwell RHS の物理内容

v3.1 は v2.0 で定めた Einstein–Maxwell–GH field RHS の物理式を、時間更新表現の変更だけを理由に変更しない。

すなわち `F_h(Ψ)` は、採用した GH–Maxwell 表現において

- `g_ab`
- `Pi_ab`
- `Phi_iab`
- `H_a`
- `Theta_a`
- `E^i`
- `B^i`

の全 RHS を同時に返す。

ただし v3.1 では、その RHS を `K` へ無理に因数分解してから exponentiationする必要はない。

\[
\boxed{
F_h(Ψ)\to \Delta u\to S[\Delta u,D_τ]
}
\]

という順序を正本とする。

物理 RHS の変更は、加法参照系・乗法正式系の双方に同時に適用されなければならない。

---

# 16. 空間離散化と境界処理

空間離散化・domain/interface・境界処理は `F_h` の定義の一部である。

これらを乗法化のために新設計しない。

正式実験では、採用した数値 Einstein–Maxwell solver の

- derivative stencil / spectral derivative
- interpolation / prolongation
- interface treatment
- outer boundary
- numerical dissipation

をそのまま共通 `F_h` に使用する。

v3.2 では、点の発散を避けるための excision（計算領域からの切り取り）や地平面境界条件を導入しない（§0A.2）。a,b の種の厳密式が持つ特異的な部分は、初期化仕様 v3.1 に従い、厳密式のまま扱う。

二経路比較で重要なのは、それらが**完全に同一**であることである。

再現性に影響する固定係数は `Ψ_numerics` または `Ψ_grid` の恒等状態へ記録する。

---

# 17. 初期化との接続契約

初期化器は、一つの初期状態

\[
\boxed{Ψ_0}
\]

だけを生成する。

加法系用と乗法系用の異なる初期化を作ってはならない。

実験開始時に `Ψ_0` をbitwise copyして、

\[
Ψ_0^{(+)}=Ψ_0^{(×)}
\]

とする。

以後は両経路を独立進行させる。

---

# 18. readoutとの分離

readout は両経路について同じ関数

\[
\mathcal R(Ψ)
\]

を適用する。

比較する主要物理量：

- a,b を囲む閉曲面上の閉じた位相として読む質量・電荷・角運動量（地平面は使わない）
- separation / orbital radius
- phase / angular velocity（系の時計 \(\Omega_{\rm sys}\)）
- `δτ` の系列
- a,b の内部時計 \(\omega_{\rm int}=M/J\) と、系の時計とのロック比
- GW flux
- EM flux
- constraint norms

詳細は読出し仕様 v3.1 に定める。

である。

まずstate-level equivalenceを検証し、その後readout-level equivalenceを検証する。

readout差だけでstateを補正してはならない。

---

# 19. 正式実験の比較ログ

最低限、各stepについて次を保存する。

```text
step
T_tau_mul
T_tau_add
tau_mul
tau_add
delta_tau
D_tau
kappa_tau
N_tick
Omega_sys
P_sys
stage2_state_error_inf
stage3_state_error_inf
stage4_state_error_inf
k1_error_inf
k2_error_inf
k3_error_inf
k4_error_inf
final_state_error_inf
final_state_error_rel
identity_state_error_inf
time_state_consistency_error
```

さらに物理 readout 差を別表に保存する。

---

# 20. FAIL code

最低限：

```text
STATE_ONE_INVALID
TIME_STATE_INVALID
TIME_STATE_INCONSISTENT
IDENTITY_STATE_DRIFT
RHS_NONFINITE
MULTIPLICATIVE_STAGE_MISMATCH
MULTIPLICATIVE_RHS_MISMATCH
MULTIPLICATIVE_ADDITIVE_MISMATCH
INITIAL_STATE_MISMATCH
NUMERICAL_CONFIG_MISMATCH
```

FAIL時に自動補正して続行してはならない。

---

# 21. v2.0 から撤回する項目

v3.1 では、v2.0 の以下を**万能相互作用の必須条件としては撤回**する。

1. `K(Ψ)` の21 blockへRHSを分解すること。
2. `KΨ=F_direct` を万能相互作用成立の中心監査とすること。
3. CF4を唯一のtime updateとすること。
4. `S=e^Be^A` を唯一のS生成法とすること。
5. matrix exponential action `AL_MOHY_HIGHAM_2011_EXPMV` を必須とすること。

これらは将来、乗法写像の別factorizationとして比較対象にしてよいが、v3.1 baselineでは使用しない。

---

# 22. v2.0 から維持する原則

以下は維持する。

1. 状態以外の hidden mutable physical/dynamical state を次状態決定へ使用しない。solver workspace、ghost buffer、FFT plan、preconditioner、allocation/cache 等の計算補助情報は、物理状態を追加せず同じ状態・同じ固定設定に対する演算を再現する範囲で許される。
2. case-specific な物理方程式・RHS・interaction・正式 time-evolution route を作らない。case情報を初期条件、topology、validation、readout、出力、停止条件に使うことは許される。
3. 初期物理ラベルを現在の dynamics係数へ戻さない。
4. readoutを dynamicsへ戻さない。
5. `ONE` を明示状態として持つ。
6. `T_τ,D_τ,κ_τ` を状態として持つ。
7. 数値再現性に必要な固定設定を状態へ保存する。
8. Einstein–Maxwell は全ケースで同一方程式を用いる。

v3.1 はこれに `δτ` の明示状態化と、RHS非破壊契約・許容実装範囲の明確化を追加した。v3.2 は、目的の書き換え（§0A.1）、単位 \(G=q_0=c=1\)（§0A.3）、`δτ` を系の時計から定める規則（§3.2.1）、4 段の位置づけ（§5）、excision・地平面境界条件の除去（§16）を追加する。

---

# 23. 万能相互作用の数学的意味

一般の加法一step map

\[
u' = u+\Delta u(u)
\]

は、`ONE=1` を追加した拡張状態

\[
\tilde u=(ONE,u)^T
\]

に対して

\[
\boxed{
\tilde u'
=
\begin{pmatrix}
1 & 0\\
\Delta u(u) & I
\end{pmatrix}
\tilde u
}
\]

と書ける。

従って任意の既存差分時間発展は、同じ一step incrementを使って状態依存乗法写像へ持ち上げられる。

**v3.2 訂正**：v3.1 はここで「既存加法更新 ≡ 状態依存乗法更新が full Einstein–Maxwell 数値発展の全状態について成立するか」を本研究で検証すべき内容としていた。これは誤りである。上の式が示すとおり、この一致は恒等的に成り立つので、検証の対象にならない。加法系との一致は実装検算にとどめる。

本研究で乗法形式を採る理由は、位相と振幅の解釈が容易になることである（§0A.1）。本実験の目的は §0A.1 に定める。

---

# 24. 実装開始前監査

正式実験コードを実装・実行する前に、以下が**本仕様または明示参照された上位・下位仕様／採用solver契約のいずれかで一意に確定**していることを確認する。すべてを本ファイルへ重複記載する必要はない。

- [ ] 共通 `F_h` の採用実装が確定
- [ ] `Ψ` state schema が確定
- [ ] `ONE,Tτ,δτ,Dτ,κτ` がstateに存在
- [ ] RK4 tableauがstateに存在
- [ ] 加法参照系の式が確定
- [ ] 乗法stage operator `S[q,d]` が確定
- [ ] final `S(Ψ_n)` が確定
- [ ] stage-level equivalence audit が確定
- [ ] identity-state audit が確定
- [ ] time-state consistency audit が確定
- [ ] 初期化器が単一 `Ψ0` を返す
- [ ] 両経路で同じ `F_h` を使用することが静的に監査可能
- [ ] 単位 \(G=q_0=c=1\)、Gauss 幾何単位が全コードで統一されている
- [ ] 系の時計 \(\Omega_{\rm sys}(Ψ)\) の推定式（現在状態だけから評価、履歴不使用）が一意に確定している
- [ ] `δτ` の規則 \(\delta\tau=P_{\rm sys}/N_{\rm tick}\)、\(N_{\rm tick}=31\) が確定している
- [ ] excision・地平面境界条件・正則性ゲートが使われていない

一つでも実装者の裁量で物理結果を変え得る未定義が残る場合は、その項目を確定してから正式実験へ進む。一方、通常のsolver workspace、memory layout、並列化、cache等、物理状態遷移を変えない実装詳細まで本仕様で固定する必要はない。

---

# 25. 最終定義

v3.2 における万能相互作用は、現在状態から同じ Einstein–Maxwell 数値RHSを評価し、4 段の increment を構成し、系の時計から次の `δτ` を定めた後、

\[
\boxed{
S(Ψ_n)=S[\Delta u_{RK4}(Ψ_n),D_{τ,n};\ \delta\tau_{n+1}]
}
\]

として、

\[
\boxed{
Ψ_{n+1}=S(Ψ_n)Ψ_n
}
\]

と更新するものである。

加法系は

\[
Ψ_{n+1}^{(+)}=Ψ_n+\Delta Ψ_n
\]

を別メモリで計算し、**実装検算にのみ使用する**。

実装検算の条件は

\[
\boxed{
Ψ_{n+1}^{(×)}\simeq Ψ_{n+1}^{(+)}
}
\]

である。

**正式な状態列は乗法系だけを採用する。**
