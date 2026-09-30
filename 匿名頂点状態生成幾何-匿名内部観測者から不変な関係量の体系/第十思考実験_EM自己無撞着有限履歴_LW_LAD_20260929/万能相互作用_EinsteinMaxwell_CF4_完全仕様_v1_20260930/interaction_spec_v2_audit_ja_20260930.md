# 相互作用仕様 v2.0 静的監査報告

**対象:** `universal_interaction_operator_complete_spec_ja_v2.0_20260930.md`  
**日付:** 2026-09-30

## 結論

相互作用部について、`admissible ψ_n -> K(ψ_n) -> CF4 -> S(ψ_n) -> ψ_{n+1}` を一意に定めるために必要な式と演算契約を v2.0 へ明示した。

この監査は **初期化仕様の完成を意味しない**。また `D_i,V_*,P_*` の生成法そのものは grid/operator construction 仕様へ分離した。ただし相互作用から見てこれら完成済み operator は `ψ_grid` の必須恒等状態であり、外部隠れ自由度を禁止した。

## v1.0 から修正した主項目

1. `∂_0 g_ab = β^i Φ_iab - α Π_ab` を明示し、4D Christoffel の生成を一意化。
2. `K_ij = 1/2 Π_ij + 1/2 n^a(Φ_ija+Φ_jia)` を明示。
3. gauge-driver から `∂_0 H_a` を一意に作り、`∇_(a H_b)` を完全定義。
4. Maxwell の contravariant state `E^i,B^i` に対して curved-space metric derivative contribution `∂_j(αγ_kl)` を明示。
5. Lie derivative は partial derivative 形式に固定し、covariant derivative との混在を禁止。
6. Maxwell `F_ab` 再構成の4D volume-form orientation と spatial vector 4D extension を固定。
7. Einstein matter source を必ず trace-reversed `T_ab - 1/2 g_ab T` から形成。
8. canonical affine lift を「spatial derivative が直接 field に作用する項だけ dynamic block、それ以外は全部 ONE 列」に固定。
9. nonlinear multiplication を `P_* diag(f_*) V_*` に固定し、retained coefficient 直接積を禁止。
10. `D_i,V_*,P_*` を state contract に昇格し、外部 boundary/ghost/penalty 依存を禁止。
11. CF4 の `K_i` と `k_i=hK_i` を区別し、積順序 `e^B e^A` を固定。
12. exponential action を `AL_MOHY_HIGHAM_2011_EXPMV` に固定。
13. 単位規約 `c=1, 4π ε0=1` を interaction 内で固定し、silent unit conversion を禁止。
14. v1.0 の `rac12` 転記誤りを `1/2` に修正。

## 外部標準との照合

### GH first-order equations

Lindblom–Scheel–Kidder–Owen–Rinne 系および SpECTRE `gh::TimeDerivative` の `g_ab, Π_ab, Φ_iab` RHS と照合。

### Gauge driver

Lindblom & Szilagyi (2009) の first-order gauge-driver

- `∂t H_a - β^k ∂k H_a = -μ(H_a-F_a)+θ_a`
- `∂t θ_a + η β^k ∂k H_a = -η θ_a`

および damped-wave target と照合。

### Maxwell

3+1 Maxwell electrovac の

- `(∂t-L_β)E^i = ε^{ijk}D_j(αB_k)+α K E^i`
- `(∂t-L_β)B^i = -ε^{ijk}D_j(αE_k)+α K B^i`

と照合し、contravariant field representation に必要な metric derivative term を明示展開した。

### CF4

Celledoni–Marthinsen–Owren の fourth-order commutator-free scheme と stage/final exponential coefficients を照合。

### expmv

Al-Mohy & Higham (2011) の scaling + truncated Taylor matrix-exponential action を唯一の exponential action method とした。

## 自動静的チェック

以下をすべて PASS：

- GH `g` RHS present
- GH `Pi` RHS present
- GH `Phi` RHS present
- gauge `H` RHS present
- gauge `Theta` RHS present
- Maxwell `E` RHS present
- Maxwell `B` RHS present
- `∇_a H_b` construction present
- Maxwell metric derivative term present
- trace-reversed matter source present
- operator-state contract present
- canonical K audit present
- CF4 present
- fixed expmv method present

## 明確に対象外のもの

本仕様の対象外は以下のみ。

1. 初期化器の具体実装
2. `D_i,V_*,P_*` を生成する grid/operator construction
3. horizon/GW/EM readout

これらは `S(ψ)` の数式を勝手に補う余地ではなく、それぞれ独立した仕様として完成させる必要がある。

## 実装禁止条件

上記対象外の仕様が未完成であっても、相互作用コードの式を別方式へ変更して補ってはならない。
