# Paper 13: Position Had Never Been Read — Sign Sectors, Duality Theorems, the Readout Hierarchy, and the Derived Clock

**Author**: Noriaki Kihara
**Date**: June 12, 2026
**Version**: v1.1 (content and results unchanged from v1.0; nomenclature unified to Paper 0.5’s displacement-record definition — ledger/register/bookkeeping mapped to conserved quantity / invariant / accounting)
**DOI (this version)**: 10.5281/zenodo.20690058 / **Concept DOI**: 10.5281/zenodo.20665633
**License**: CC BY 4.0
**Series**: Dual Geometry of Wavelength Space and Frequency Space (Paper 13, sequel to Papers 1–12)

---

## Abstract

"Can discrete states be mapped back and forth to position and time?" — we answer the oldest question posed to this series' system (three axioms: νλ=1, zero point ½, one asymmetric bit, plus configuration reading and the record principle) with **exactly zero additions** to the axiom inventory. Main results: (1) **Position is not a new degree of freedom — it was the unread part of the dictionary.** Discrete position = the sign sector of the wave (quarter-period offsets: cos→sin→−cos→−sin), which is merely a change in how the data already present in the dictionary of Paper 5 is read out — the addition to the state space is exactly zero (Theorem 1). (2) **Space is derived**: position space = the Pontryagin dual $T^4$ (character group) of the conserved-quantity lattice (charge-like) $\mathbb{Z}^4$, and the three symmetries — translation, rotation, inversion — are **theorems of duality**, not assumptions (Definition 2, Theorems 2a–2b). No substantial background spacetime is assumed — position space is a chart on relational data (line coefficients). (3) **The readout has three strata**: the scalar data of intensity fringes (spacings) collapse the 21 relation classes to only 11, line-position data to 13, while **branch channels (complex line coefficients) separate all 21 orbits completely** (Theorem 3; exhaustive counts, 11 ⊂ 13 ⊂ 21) — the quantitative demonstration that "position had never been read," locating exactly which information is lost in intensity-only readings. (4) **The quarter-displacement theorem**: a fragment's quarter-lattice displacement is read exactly from branch-channel phase classes (Z₄), with the general component rule $m\bmod4$ (Theorem 4) — continuous position is the hierarchical completion of these quarter digits (the explicit map is Paper 15, to be published concurrently). (5) **Time is resolved — but asymmetrically to position**: the targets $\nu_t^2$ of the clock frequency $\nu_t=\sqrt{s}$ run over **all odd integers** with spacing 2 — a **proven theorem of unlimited range via Gauss's triangular-number theorem** (Theorem 5: for odd $s$, $(s-1)/2$ is a sum of four triangular numbers, always solvable by Gauss's theorem that every natural number is a sum of three triangular numbers). The most regular spectrum possible stands as a proven claim. The forward map is exact on all four axes, while **t admits no independent inverse map** (Theorem 6): $\nu_t$ is derived as the norm of the spatial configuration, so t has no independent scale to round. This is not a defect but the theorem form of the all-axes-symmetric principle — "t merely appears to have been selected by observation." t is **easier** than xyz. This paper claims no re-derivation of standard theory. No physical identification is made.

---

## 1. Introduction

### 1.1 The fifty-year question

How do continuous position $x$ and time $t$ coexist with a discrete state description? For the present system the question takes the form:

> Can position and time be **read out** of the system's state data (occupied cells, branches, conserved quantities)? If so, does the readout require an addition to the state space (a new degree of freedom)?

This paper's answer: **they can be read, with zero additions.** Position is a re-reading of the sign sectors the dictionary (Paper 5) has carried from the start; time is a derived quantity of the conserved norm. What was needed was neither an axiom nor an extension but a **change of readout principle** — hence the title.

### 1.2 What is not claimed

(i) The explicit construction of the map to continuous coordinates (mixed-radix expansion, readout operation, inverse map) is carried out in Paper 15 (concurrent publication) — this paper is its **theorem-level foundation** (what can be read; what is derived). (ii) The choice of time direction (why some axis looks like t) and the status of amplitudes are Paper 14. (iii) **This paper does not select an aggregation rule (measure) for path amplitudes** (the measure line = Paper 16, concurrent publication). (iv) No physical identification is made.

## 2. Preliminaries (minimal recap)

The dictionary of Paper 5: cell $k\in\mathbb{Z}^4$ ↔ product wave $\Phi_k=\prod_j\varphi_{k_j}(x_j)$, with branch $k_j>0$ = cos and $k_j<0$ = sin. Configuration reading (Paper 7): the state is the set of occupied cells = the wave itself. The record (Papers 6/12): the hologram $I=\Psi^2$, whose physical content is the line coefficients (relational data). Conserved quantities: $s=\sum_j(|k_j|+\tfrac12)^2$ and $\varepsilon=(-1)^{\Sigma|k_j|}$.

**Gauge data held fixed (complete enumeration)**: the marking $u$ (time-read axis), the per-axis orientations $Z_2^{4}$ (positive branch directions), and the reference fragment (origin). The sign-sector / Z₄ readings of Theorems 1 and 4 are **gauge-relative quantities** with respect to these (limitation 1 of §6 connects here).

**Internal terminology declaration**: "position" and "time" in this paper are names of readout quantities internal to the system and carry no identification with physical quantities (naming discipline: conserved quantities carry the "-like" flag, e.g. charge-like).

## 3. The resolution of position

### 3.1 Theorem 1 (sign sector = discrete position; zero additions)

> **Theorem 1**: the quarter-period translation $T_{1/4}$ of the fundamental wave cycles the dictionary's (branch, sign) data: $\cos\to\sin\to-\cos\to-\sin$ ($d\to d+1\bmod4$). Hence **the (branch, sign) pair of the dictionary = quarter-period offset = discrete position**, and introducing position requires not a single new state.

**Proof** (two lines): $\cos(\theta-\pi/2)=\sin\theta$ and $\sin(\theta-\pi/2)=-\cos\theta$ — the quarter translation is the cycle of these trigonometric identities. ∎ Machine confirmation (②): exhaustive 4/4 (Fig. 1a). The sign sectors hitherto treated as "alternative presentations of the same state" were position data all along — only the readout principle had to change.

### 3.2 Space = the dual of the conserved-quantity lattice (split into definition, general theorem, and consistency theorem)

**Definition 2**: position space is **defined** as the character group (Pontryagin dual) $T^4=\mathrm{Hom}(\mathbb{Z}^4,U(1))$ of the conserved-quantity lattice (charge-like) $\mathbb{Z}^4$ (the dual chart).

> **Theorem 2a (general — a matter of proof, not of verification)**: under this definition, (i) translation symmetry = the group structure of the character group, (ii) rotation (B₄) symmetry = the dual action of lattice automorphisms, (iii) inversion symmetry = complex conjugation — all standard theorems of Pontryagin duality.

> **Theorem 2b (model-specific consistency)**: the physical content of the record (line coefficients) is a function on $\mathbb{Z}^4$, and the sign-sector reading of Theorem 1 embeds into this dual chart as the per-axis $Z_4\subset U(1)$ (quarter translation = quarter rotation of characters). Machine confirmation (②): the 4/4 cycle of Theorem 1 is precisely the consistency check of $Z_4\hookrightarrow U(1)$.

Consequence (background independence): position space is not a substantial container assumed in advance but the dual chart of relational data (lattice-indexed line coefficients). Continuous position is obtained as the hierarchical completion of quarter digits (Fig. 1b; explicit construction in Paper 15, Theorem 2, concurrent publication).

### 3.3 Theorem 3 (the readout hierarchy 11 ⊂ 13 ⊂ 21)

The relation classes (B₄×swap orbits) of all pairs of the 64 shell-9 cells (2016 pairs) number **21**. Separating power per readout stratum (exhaustive counts; Fig. 2):

| Readout stratum | Readable data | Separated classes |
|---|---|---|
| Scalar fringe | spacings $(|\Delta|^2,|\Sigma|^2)$, unordered | **11** (18 even when ordered) |
| Line-position fringe | unsigned line pair $\{[u{+}v],[u{-}v]\}$ | **13** |
| Branch channels | complex line coefficients (cos/sin resolution) | **21 = complete separation** |

> Intensity-only reading (the traditional fringe) irrecoverably loses **10 classes' worth of distinctions at the scalar stratum (21→11; the pair structure of these mergers is unverified) and 8 pairs at the line-position stratum (verified to be pairs)** out of the 21 orbits. **The only reading that recovers all information is the branch channel** — the quantitative form of "position had never been read." All 8 merged pairs of the line-position stratum are of the type "sum line confused with difference line = unreadable relative sign."

### 3.4 Theorem 4 (readability of quarter displacements — the component rule)

> **Theorem 4**: a fragment's quarter-lattice displacement shifts the phase class (Z₄) of the line of frequency component $m$ by **$m\bmod4$** (general rule). For the fundamental ($m=1$) the shift is +1; for odd components ($m\bmod4\in\{1,3\}$) the shift is invertible, so the displacement is readable.

"+1" is restricted to the fundamental (per-component $m\bmod4$ for multi-component fragments). This is consistent with Paper 15's digit conventions ($d^{(0)}$, the $\min(c,k{-}1)$ clip) as the $m=1$ section of the general rule. Machine confirmation (③, sampled): the two fringe-shift tests (tripod quarter rotation, ring inversion) PASS — sampled, not exhaustive. This is the operational readout of discrete position; the hierarchical extension (peeling decode) and the continuum are Paper 15, §5.

## 4. The resolution of time

### 4.1 Theorem 5 (clock targets = all odd integers — with proof)

> **Theorem 5**: the range of the squared clock frequency $\nu_t^2=s=\sum_j(|k_j|+\tfrac12)^2$ is **exactly the set of all odd integers** (spacing 2).

**Proof**: (range ⊆ odds) $s=\sum_j|k_j|(|k_j|+1)+1$, and each $|k_j|(|k_j|+1)$ is a product of consecutive integers, hence **even** — so $s$ is always odd. (odds ⊆ range) For any odd $s$, $(s-1)/2=\sum_j|k_j|(|k_j|+1)/2=\sum_j T_{|k_j|}$ is a decomposition into **four triangular numbers**, always solvable by Gauss's triangular-number theorem (every natural number is a sum of three triangular numbers; the 1796 "Eureka"), adding $T_0=0$ as the fourth. ∎ (Identification of the proof: first review, point #1, claude.ai.)

Machine confirmation (②): all 100/100 odd integers ≤ 199 reached (Fig. 3a). The clock's targets are not "sparse and irregular" but the most regular possible — **as a proven claim**. The tick at each $s$ is $\Delta t=1/\sqrt{s}$ (Fig. 3b) — the same construction as the internal time $t=\sum1/\sqrt{s}$ of Paper 8.

### 4.2 Theorem 6 (no independent inverse map for t)

> **Theorem 6 (factorization of t)**: the forward map (state → readings of four axes) exists exactly on all four axes. In the reverse direction: **any t-value readable from the record factorizes through a function of $s$** — there exists no readout granting t a resolution independent of $s$. The independent rounding (inverse map, Paper 15, D4) that each xyz axis possesses has **no counterpart** for t.

**Proof sketch (inventory audit)**: the accounting of state data is exhausted by occupied cells (four integer axes), branches, and signs (configuration reading, §2). The only time-involving readings are the norm clock $\nu_t=\sqrt{s}$ and the event count (Theorem 5; Paper 8); every t-dependent quantity in the record is of the form "tick count × rate $1/\sqrt{s}$" — both factors functions of $s$ and the count. No fifth accounting entry exists in the inventory that could supply a t-scale independent of $s$ (finite enumeration). ∎ (sketch)

The composite round-trip check of Paper 15 (500/500) is the **confirmation of the realization** of this theorem (the distinction between properties of a construction and theorems).

This is a signature, not a defect: **t is not a fundamental degree of freedom**. All axes (including R and Q) are symmetric and interchangeable; t merely appears to have been selected by observation — the present theorem is the exact expression, at the level of round-trip maps, of that principle (whose gauge form — which axis looks like a clock — is the marking theorem of Paper 14). Ironically, t is **easier** than xyz precisely because there is nothing to round.

## 5. Consequence: the present state of the fifty-year question

| Question | Answer | This paper |
|---|---|---|
| Can position be read out? | Yes (branch channels) — zero additions | Theorems 1, 3, 4 |
| Where do space and its symmetries come from? | Theorems of the dual of the conserved-quantity lattice (charge-like) | Definition 2, Theorems 2a–2b |
| Can time be read out? | Forward map exact; t is derived (no independent inverse) | Theorems 5, 6 |
| Explicit map to continuous coordinates | Constructible | Paper 15 |
| Choice of time direction; amplitudes | Theorematizable | Paper 14 |

## 6. Honest limits

1. "Position" here is the reading of a single fragment (relative to a reference fragment). Common time for multiple fragments (synchronization) belongs to the gauge-fixing mechanism of the record interface and is not constructed.
2. The counts of Theorem 3 are exhaustive at shell 9 (s=9). The separation counts of the readout hierarchy at general s are uncounted (mechanically extensible).
3. Convergence and error guarantees of the continuum completion are deferred to Paper 15, Theorem 2.
4. No physical identification is made.

## Figures

- **Figure 1** (`paper13_fig1_sign_sector.png`): (a) the four quarter-translates of the fundamental (cos/sin/−cos/−sin) — position = sign sector. (b) nested hierarchical quarter lattices (refinement of index 3 per level; depth 0–2 shown)
- **Figure 2** (`paper13_fig2_readout_layers.png`): separating power of the three readout strata, 11 ⊂ 13 ⊂ 21 (exhaustive count over all 2016 pairs)
- **Figure 3** (`paper13_fig3_derived_time.png`): (a) clock targets = all odds (100/100 below 199, spacing 2). (b) tick $\Delta t=1/\sqrt{s}$ — no independent scale for t

All figures are exact computations (no schematics).

## Appendix A: Verification summary (labels: ① true by construction / ② theorem + confirmation / ③ check that could have failed)

| Claim | Label | Verification | Method | Script |
|---|---|---|---|---|
| Theorem 1 (quarter translation = digit shift) | ② (proof given) | 4/4 | exhaustive confirmation | supplement62_projection_formalization.py |
| Theorem 2a (duality, general) | — (proof matter; not subject to verification) | — | — | — |
| Theorem 2b ($Z_4\hookrightarrow U(1)$ consistency) | ② | same 4/4 as Theorem 1 | exhaustive confirmation | same |
| Theorem 3 (11/13/21) | ③ | all 2016 pairs | **exhaustive** | supplement79_scale_free_and_stripe_table.py, supplement81_mirror_lemma_check.py (scalar fringe), paper13_14_figures.py |
| Complete separation of 21 orbits (branch channels) | ③ | 21/21 | exhaustive | supplement44_branch_channel_verification.py |
| Theorem 4 (readability of displacements) | ③ | 2 fringe-shift tests PASS | **sampled** | supplement44_branch_channel_verification.py |
| Theorem 5 (all odds) | ② (proof by Gauss's theorem) | 100/100 odds ≤ 199 | in-range confirmation | paper13_14_figures.py |
| Theorem 6 (t factorization) | ② (inventory-audit sketch) | Paper 15 round trip 500/500 as realization check | cross-check | supplement62_inverse_roundtrip_test.py |

---

**Acknowledgments / history**: verification was carried out under the two-party independent verification protocol of Claude Code (local machine verification) and claude.ai (independent re-computation). Precise record of verification scope: the three-stratum counts on the 64-cell domain (11/13/21) are Claude Code verifications; claude.ai confirmed method equivalence on the 48-cell domain (13/8) — the independent 64-domain recomputation is registered in claude.ai's queue. The proof of Theorem 5 is due to the identification of Gauss's triangular-number theorem (claude.ai, first review #1). The prototype of Theorem 6 is Supplement 48 (June 12, 2026). Source material: Supplements 43, 44, 47 (branch-channel part only), 48, 79, 81 (June 12, 2026). v0.2: 14 points of the first review (major revision). v1.0: points R1–R5 of the second review (accept with minor conditions).

No physical identification is made.
