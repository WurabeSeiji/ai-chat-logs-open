# 論文0 v1.4 挿入確定テキスト(claude.ai 査読確定版 → Claude Code 反映用)

**性格**: 先行研究調査(Claude Code)を claude.ai が査読し、挿入位置・文言を確定したもの。Claude Code はこれを §1/§4.6/参考文献に機械的に反映し、付録A数式修正(済)・related_identifiers 拡張と併せて v1.4 として再ビルド → newversion 投稿。二者検算プロトコル。

---

## 挿入1:§1 末尾(存在閾値 (2) の直後)に一文追加

**日本語**:
> なお、辺長一定の正則多胞体が球面上に存在するための条件は、Schläfli 行列・角欠損による古典的制約として知られる [Coxeter 1973]。本稿の存在閾値 $R^{*}_d$ はその系譜に属するが、本稿の新規性は閾値そのものではなく、これを次元の天井(§4.6)として読む点にある。

**English**:
> The condition for a regular polytope of fixed edge length to exist on a sphere is classically known as a constraint via the Schläfli matrix and angular deficit [Coxeter 1973]. The existence threshold $R^{*}_d$ here belongs to that lineage; the novelty of this paper is not the threshold itself but reading it as a ceiling on dimension (§4.6).

---

## 挿入2:§4.6「シリーズへの含意」の直前に新小節「先行研究との関係」を追加

**日本語**:
### 先行研究との関係

「次元が一意でなくスケールに依存しうる」という着想自体は新しくない。拡散で測るスペクトル次元は量子重力(CDT・漸近安全性・LQG)で短距離 $\sim 2$・長距離 $\sim 4$ へ流れることが知られ [Carlip 2019; Ambjørn–Jurkiewicz–Loll 2005]、これは連続的・拡散的な実効次元である。次元を連続パラメータとして扱う技法(次元正則化)や、測度から定まる非整数次元(フラクタル次元)も、次元を整数に限らない点で関連する。定曲率空間で曲率半径が幾何的限界を半径依存で与える点では、Böröczky のパッキング密度上界 [Böröczky 1978] が最も近い先例である。

本稿が異なるのは三点である。(i) 天井の源泉が拡散でも密度でもなく、「辺長1の正則測地セルが整数個収まるか」という**離散的存在条件** $d\sin^2(1/2R)\le 1$ であること(Böröczky は密度上界であって次元の天井ではなく、本稿の対象とは異なる)。(ii) これが次元の共役量 $\kappa=\sin^2(1/2R)$ と**容量関係** $d\cdot\kappa\le 1$ として、$\nu\lambda=1$・$\sum\nu^2=\mathcal{N}^2$(論文1)と同型の「保存予算の飽和」に乗ること。共役変数の既存の文献は Fourier/シンプレクティック対(位置-運動量、エネルギー-時間)に限られ、次元を曲率予算と共役にする定式化は見当たらない。(iii) 幾何天井と論文11 の検閲天井という**独立2機構が** $d=4$ **で臨界一致**すること((√4−1)²/2=½ の等号)。「4次元を導く」議論は多数あるが [Ehrenfest 1917; Tegmark 1997 ほか]、独立な二機構の臨界交差として4を導く構成は見当たらない。

本稿の主張はこの**連鎖の具体性**にあり、「次元の曖昧さ」一般の主張ではない。上記の先行研究との差は、いずれも本シリーズの数え上げ的出発点(整数セルの存在)に由来する。

**English**:
### Relation to prior work

The idea that "dimension need not be unique and can depend on scale" is itself not new. The spectral dimension measured by diffusion is known to flow from $\sim 2$ at short distances to $\sim 4$ at long distances in quantum gravity (CDT, asymptotic safety, LQG) [Carlip 2019; Ambjørn–Jurkiewicz–Loll 2005]; this is a continuous, diffusive effective dimension. Techniques treating dimension as a continuous parameter (dimensional regularization) and non-integer dimensions defined from a measure (fractal dimension) are also related in that dimension is not restricted to integers. Where a curvature radius gives a radius-dependent geometric limit in constant-curvature space, the closest precedent is Böröczky's packing-density upper bound [Böröczky 1978].

This paper differs in three respects. (i) The source of the ceiling is neither diffusion nor density but a **discrete existence condition** $d\sin^2(1/2R)\le 1$ — "do integer-many regular geodesic cells of edge 1 fit?" (Böröczky is a density bound, not a ceiling on dimension, a different object). (ii) This rides, as the conjugate quantity of dimension $\kappa=\sin^2(1/2R)$ and the **capacity relation** $d\cdot\kappa\le 1$, on the same "saturation of a conserved budget" as $\nu\lambda=1$ and $\sum\nu^2=\mathcal{N}^2$ (Paper 1); the existing literature on conjugate variables is confined to Fourier/symplectic pairs (position–momentum, energy–time), and a formulation making dimension conjugate to a curvature budget is not found. (iii) The geometric ceiling and Paper 11's censorship ceiling — **two independent mechanisms — coincide critically at** $d=4$ (the equality (√4−1)²/2=½). Arguments that "derive four dimensions" are numerous [Ehrenfest 1917; Tegmark 1997, etc.], but a construction deriving 4 as the critical crossing of two independent mechanisms is not found.

The claim of this paper lies in the **specificity of this chain**, not in any general claim about "dimensional ambiguity." Each difference from prior work stems from the counting-based starting point of this series (the existence of integer cells).

---

## 挿入3:参考文献節を新設(論文0 末尾、付録Bの後)

**冒頭注記(日英)**:

**日本語**:
> **参考文献について**:本シリーズ(論文1〜16・総説)は外部文献を意図的に引用しない方針である(総説 参考文献節)。本稿(論文0・基礎篇)はこの方針の例外であり、隣接領域(実効次元・球パッキング・次元性)との関係を明示する性格上、最小限の外部参照を持つ。これは方針変更ではなく、基礎篇が隣接領域との境界を画定する役割に由来する論文0 固有の扱いである。

**English**:
> **On references**: This series (Papers 1–16 and the survey) deliberately cites no external literature (survey, references section). This paper (Paper 0, Foundations volume) is an exception: by its nature of making explicit the relation to adjacent areas (effective dimension, sphere packing, dimensionality), it carries a minimal set of external references. This is not a change of policy but a treatment specific to Paper 0, arising from the Foundations volume's role of demarcating the boundary with adjacent areas.

**書誌(確定)**:
- [Coxeter 1973] H. S. M. Coxeter, *Regular Polytopes*, 3rd ed., Dover, 1973.
- [Carlip 2019] S. Carlip, "Dimension and Dimensional Reduction in Quantum Gravity," *Class. Quantum Grav.* 34 (2017) 193001; arXiv:1904.04379.
- [Ambjørn–Jurkiewicz–Loll 2005] J. Ambjørn, J. Jurkiewicz, R. Loll, "Spectral Dimension of the Universe," *Phys. Rev. Lett.* 95 (2005) 171301.
- [Böröczky 1978] K. Böröczky, "Packing of spheres in spaces of constant curvature," *Acta Math. Acad. Sci. Hungar.* 32 (1978) 243–261.
- [Ehrenfest 1917] P. Ehrenfest, "In what way does it become manifest in the fundamental laws of physics that space has three dimensions?" *Proc. Amsterdam Acad.* 20 (1917) 200.
- [Tegmark 1997] M. Tegmark, "On the dimensionality of spacetime," *Class. Quantum Grav.* 14 (1997) L69–L75.

> **注(書誌の検証)**:Carlip 総説の巻号・年(Class. Quantum Grav. の掲載年と arXiv:1904.04379 の対応)は、Claude Code が公開前に原典で最終確認すること。他5件も DOI/巻号を原典照合のうえ確定。claude.ai は文献の*役割と差別化*の正確性までを保証し、書誌情報の最終照合は Claude Code 管轄とする(二者分担)。

---

## 反映後の作業(Claude Code)

1. 挿入1〜3 を ja/en 双方に反映、付録A数式修正(済)を含めて v1.4 として再ビルド。
2. 書誌6件を原典照合(特に Carlip の巻号・年)。
3. related_identifiers 拡張:論文5・9・11(既定)に加え、本版は外部文献を参考文献に持つため、Zenodo メタデータの references フィールドに6件を登録(任意だが推奨)。
4. newversion 投稿(DOI 20680270 系列の更新版)。Concept DOI 20680269 は最新版へ自動転送。
5. 投稿後、Version DOI をヘッダに更新。

## claude.ai 査読確定事項

- §4.6 段落:Claude Code 案を採用、「次元正則化/フラクタル次元=同系」を「関連の質を分けて」精密化(技法 vs 非整数次元)。
- 差別化三点(離散的存在条件/曲率予算共役/二天井臨界一致)は調査でも先例が見当たらず、新規性主張として堅い。
- 「次元の曖昧さ=新概念」とは打ち出さない方針を確定(先例で反証されるため)。主張は連鎖の具体性。
- 物理的同一視(ビッグクランチ等)は本文に書かない。R(t) 動力学は §4.7 のとおり範囲外を維持。
