# 万能相互作用 Einstein–Maxwell 乗法更新・加法検算 完全仕様 v3.0
## — 既存の数値 Einstein–Maxwell 差分発展を、状態依存乗法写像 S(Ψ) として厳密に表現し、加法参照系と逐step検算する —

**著者:** 木原範昭 (Noriaki Kihara)  
**仕様整理:** ChatGPT  
**版:** v3.0  
**日付:** 2026-09-30  
**置換対象:** `universal_interaction_operator_complete_spec_ja_v2.0_20260930.md` の時間更新・相互作用表現部  
**重要:** 本仕様はコード実装ではない。実装開始前の正本仕様である。

---

# 0. 改定理由

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

---

# 1. 最上位目的

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

完成条件は、同じ初期状態、同じ空間離散化、同じ境界処理、同じ時間刻み、同じ時間積分公式に対して

\[
\boxed{
Ψ_{n+1}^{(\times)}=Ψ_{n+1}^{(+)}
}
\]

が数値丸め誤差の範囲で毎step成立することである。

---

# 2. 絶対ルール

1. **正式系は乗法系である。** 次状態は必ず `S(Ψ)Ψ` で生成する。
2. 正式系で `Ψ += ΔΨ`、`u += Δu` を次状態更新として使用してはならない。
3. 加法系は**検算専用**であり、正式系の次状態生成へ帰還させない。
4. 乗法系と加法系は同じ Einstein–Maxwell RHS evaluator を使用する。
5. 乗法系と加法系で物理方程式、空間差分、境界条件、ゲージ、constraint処理、時間積分係数を変えてはならない。
6. `case_id` による dynamics 分岐は禁止する。
7. 質量・電荷・スピン等の初期ラベルを時間発展係数へ直接再注入しない。
8. `T_τ`, `δτ`, `D_τ`, `κ_τ`, `ONE` は状態成分である。
9. 外部 `dt`、外部現在時刻、外部step counterを物理更新へ使用してはならない。
10. 数値解法に必要な固定係数・差分係数・境界係数・時間積分係数は、仕様で固定された恒等状態として保持する。
11. readoutは一方向であり dynamics へ帰還しない。
12. 加法系との一致を得るために乗法系の定義をstep途中で変更してはならない。

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

v3.0 では `κ_τ≠0`, `T_τ>0`, `D_τ>0`, `δτ>0` を必須とする。

`δτ` は状態から読む。外部引数 `dt` を dynamics に渡してはならない。

時間更新は

\[
\boxed{T_{τ,n+1}=D_τ T_{τ,n}}
\]

とする。

`δτ,Dτ,κτ` は v3.0 では恒等状態である。

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

v3.0 の乗法化は、この物理 field representation 自体を変更しない。

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

---

# 5. 時間積分公式

v3.0 baseline は classical explicit RK4 とする。

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
\boxed{T_{τ,n+1}^{(+)}=D_τT_{τ,n}}
\]

`ONE,δτ,Dτ,κτ` とその他恒等状態は変更しない。

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
\delta τ_{n+1}=\delta τ_n
\]

\[
D_{τ,n+1}=D_{τ,n}
\]

\[
\kappa_{τ,n+1}=\kappa_{τ,n}
\]

および全固定数値設定、初期ラベル、grid設定が不変であることを監査する。

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

だけから行う。

---

# 15. Einstein–Maxwell RHS の物理内容

v3.0 は v2.0 で定めた Einstein–Maxwell–GH field RHS の物理式を、時間更新表現の変更だけを理由に変更しない。

すなわち `F_h(Ψ)` は、採用した GH–Maxwell 表現において

- `g_ab`
- `Pi_ab`
- `Phi_iab`
- `H_a`
- `Theta_a`
- `E^i`
- `B^i`

の全 RHS を同時に返す。

ただし v3.0 では、その RHS を `K` へ無理に因数分解してから exponentiationする必要はない。

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
- puncture / excision treatment
- numerical dissipation

をそのまま共通 `F_h` に使用する。

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

- horizon位置
- separation / orbital radius
- phase / angular velocity
- horizon mass / charge / spin
- GW flux
- EM flux
- constraint norms

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

v3.0 では、v2.0 の以下を**万能相互作用の必須条件としては撤回**する。

1. `K(Ψ)` の21 blockへRHSを分解すること。
2. `KΨ=F_direct` を万能相互作用成立の中心監査とすること。
3. CF4を唯一のtime updateとすること。
4. `S=e^Be^A` を唯一のS生成法とすること。
5. matrix exponential action `AL_MOHY_HIGHAM_2011_EXPMV` を必須とすること。

これらは将来、乗法写像の別factorizationとして比較対象にしてよいが、v3.0 baselineでは使用しない。

---

# 22. v2.0 から維持する原則

以下は維持する。

1. 状態以外のmutable persistent情報を読まない。
2. case-specific dynamicsを作らない。
3. 初期物理ラベルを現在の dynamics係数へ戻さない。
4. readoutを dynamicsへ戻さない。
5. `ONE` を明示状態として持つ。
6. `T_τ,D_τ,κ_τ` を状態として持つ。
7. 数値再現性に必要な固定設定を状態へ保存する。
8. Einstein–Maxwell は全ケースで同一方程式を用いる。

v3.0 はこれに `δτ` の明示状態化を追加する。

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

本研究で検証すべき内容は、

\[
\boxed{
\text{既存加法更新}
\equiv
\text{状態依存乗法更新}
}
\]

が full Einstein–Maxwell 数値発展の全状態について成立するかである。

---

# 24. 実装開始前監査

コード実装を開始してよいのは以下が全て明示された後のみ。

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

一つでも未定義なら正式実験コードを書いてはならない。

---

# 25. 最終定義

v3.0 における万能相互作用は、現在状態から同じ Einstein–Maxwell 数値RHSを評価し、RK4 incrementを構成した後、

\[
\boxed{
S(Ψ_n)=S[\Delta u_{RK4}(Ψ_n),D_τ]
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

を別メモリで計算し、**検算にのみ使用する**。

最終監査条件は

\[
\boxed{
Ψ_{n+1}^{(×)}\simeq Ψ_{n+1}^{(+)}
}
\]

である。

**正式な状態列は乗法系だけを採用する。**
