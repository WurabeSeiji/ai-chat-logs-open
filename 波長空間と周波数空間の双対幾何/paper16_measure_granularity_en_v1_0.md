# Paper 16: An Exact Reduction of the Measure Problem — Finite Decomposition into Interference Granularity, Event Structure, and Conditioning

**Author**: Noriaki Kihara
**Date**: June 12, 2026
**Version**: v1.0 (final; v0.2 incorporated the 6 mandatory points of the first review [major revision; the entrance examination — "report of the selection problem, not a pre-emption of the selection" — was passed], v1.0 the points R1–R6 of the second review [accept with minor conditions])
**Series**: Dual Geometry of Wavelength Space and Frequency Space (Paper 16, sequel to Papers 1–15)

---

## Abstract

The system supplies path amplitudes (transport phases η, Paper 14) with no free parameters — the **raw material** of probability is fully in stock. What is not yet fixed is the **aggregation rule** that bundles the material into probabilities; this is precisely what standard theory installs as the Born rule, by axiom. This paper does **not** select the aggregation rule. Its main result is an **exact reduction of the selection problem itself**: (1) **The historical proximity 0.0007 was not a meaningful agreement** (Theorem 1): the s=9 twin observable is doubly compressed — exclusion forces a two-outcome space, and both measures pile ~0.98 onto the dominant outcome — collapsing the total-variation distance; in the three-outcome space at s=11 the divergence **jumps three orders of magnitude** to 0.24. A one-line compression lemma (TV $=1-\sum\min(P,Q)\le\varepsilon$ when both distributions give the same dominant outcome mass $\ge1-\varepsilon$; the tightened form is due to the reviewer) makes the mechanism exact. (2) **Structural divergence** (Theorem 2; stated for the verified domain $s\in\{9,11,13\}$, exhaustive): derived (amplitude) measures do not classicalize to counting measures — the divergence is persistently O(0.1–0.9), robust to the convention width of the capacity factor; **the two sides share their support (the zeros coincide) and diverge only in weights**. (3) **Interference granularity is not reducible to convention** (Theorem 3): the choice of which alternatives are summed coherently — channel-consistent (the historical W) versus configuration granularity (W′) — changes values algebraically. The exact mirror example: $(-3,0,0,0)$ has $W=0$ exactly under channel-consistent aggregation yet decays perfectly normally at configuration granularity ($W'=768=W'(A)$ exactly). The effect is not scale-free (it fails at m=2), and at m=4 the **unit of cancellation is the channel** — the mirror does not vanish as a whole ($S_+-S_-\neq0$, supporting the parity hypothesis), but channel (5,7,9) **cancels by itself** ($W_A=147456$, $W_B=0$, $W'$ both 3072): per-channel $W=0\wedge W'>0$ occurs even at even m, so the supply of granularity divergences is **richer** than the odd-m subsequence; moreover **W′ mirror-equality holds for all channels at m≥3** (m=2 appears to be a small-s exception). The occurrence condition has been reduced to **algebra that is the object of proof**: the four-value lock $z^-_K/z^+_K\in\{0,\infty,\pm1\}$ (all 36 configurations, zero exceptions; exclusive families weigh equally, 384=384). (4) **Conditioning demands the branch-channel stratum** (Theorem 4): the twin measure is a function of the parent pair's relation class; the check that could have failed is that pairs indistinguishable by intensity fringes (the sign of dot) have **different** P(X) — the measure cannot be supported on strata 1–2 of the readout hierarchy 11⊂13⊂21 (Paper 13, Theorem 3) and requires stratum 3. (5) **Event structure is exactly one bit**: one-shot (merged-key) versus sequential (split-key) translates, via the record principle, into interference granularity; the historical counting dichotomy is its phase-free shadow — the 2×2 table is numerically complete. (6) **Main theorem** (Theorem 6, a classification within the family): the inventory-defined candidates organize without duplication into the product of three independent axes — granularity × event-structure bit × conditioning — each axis independently changing the measure (with explicit examples); no claim is made that candidates outside the family cannot exist. The decision experiments are **computable inside the model** (decision table, §9; items m=4 and ι executed before publication — Supplement 84; m=5@s31 honestly registered as open, the load-bearing point for the Born question). Under this reduction the Born rule has the status of an **effective law to be merged with in a suitable limit**, with the three-part acceptance criterion **convergence, suppression (statistical censorship ~1/√N), and prediction**. No physical identification is made.

---

## 1. Introduction

### 1.1 The question, and why this paper does not select

Papers 5–15 established states, conservation, time, position, phases, and coordinates as theorems or machine verifications. What remains is the rule that turns the set of amplitudes η into **probability** — which alternatives to add coherently (granularity), what to condition on, and under which event structure to read. Standard theory settles this with one axiom: the Born rule. The present system has the opposite problem: aggregation candidates are **over-derivable** — sequential and batch counting, channel-consistent and configuration-granularity amplitudes, one-shot and sequential joint amplitudes are all definable from stock, and they genuinely differ at finite s. The honest deliverable at this stage is therefore not a verdict but an **exact reduction of the selection problem with internally computable decision experiments**. That is this paper.

### 1.2 What is not claimed

(i) Selection of the aggregation rule (no result here depends on it). (ii) A derivation of the Born rule (only its status is organized — §10). (iii) The reality of deviations (only the location of their structure). (iv) No physical identification.

## 2. Definitions (the measure family at general s)

For twins $(s,s)$ (capacities $c(m)=|SH[m]|$; channels = `all_finals(s)`):

- **D1 (sequential counting)**: first-channel probability ∝ configuration count on the remaining lattice; the second is the conditional update. Well-definedness rests on the lemma that conditional counts depend only on per-shell occupancies (products of binomial capacities).
- **D2 (batch counting)**: equal weight on joint final configurations.
- **D3 (derived twin)**: $P\propto \mathrm{mult}\cdot W_{F_1}W_{F_2}\cdot[$joint capacity compatibility$]$, with $W_F=|z_F|^2$ the channel-consistent amplitude under the canonical counting of Paper 14, Theorem 5.
- **W′ (configuration granularity)**: $W'_F=\sum_K|z_K|^2$ (square per final configuration $K$, then sum).
- **D4 (tracking observable)**: $\Delta(s)=\mathrm{TV}(\text{D3},\text{D1})$; $\Delta'(s)$ is the W′ version.

Anchors (s=9), all exactly reproduced: sequential 396/403, batch 60/61, derived 0.98195, $\Delta(9)=0.00068$, single decay 24/31.

## 3. Theorem 1: the identity of 0.0007 is two-outcome compression

**Compression lemma (one line)**: if two probability distributions both give the same dominant outcome mass $\ge1-\varepsilon$, then $\mathrm{TV}=1-\sum\min(P,Q)\le\varepsilon$ (the tightened form is due to the reviewer). In the s=9 twin, (i) exclusion (uniqueness of the origin cell) forces a two-outcome space, and (ii) both measures concentrate ~0.98 on the dominant outcome ($\varepsilon\approx0.018$) — this double compression had crushed the TV (Fig. 1a). When the three-outcome space opens at s=11, $\Delta(11)=0.24136$ — a **three-order jump** (Fig. 1b). Corroboration: at the single-decay level even s=9 shows 0.774 vs 0.965; sector B has $\Delta_B(9)=0.030$. **The proximity of 0.98195 and 0.98263 was an accident of s=9 and sector A**; the historical significance attached to it is corrected here (the values themselves stand).

## 4. Theorem 2: structural divergence (no classicalization)

**Verified domain made explicit**: this theorem is a statement at the exhaustive level for $s\in\{9,11,13\}$ (a three-point series). Over this series $\Delta(s)$ runs 0.0007 → 0.24 → up to 0.90 (strong parent-type and sector dependence; Fig. 2). The internal width of the counting side (TV(sequential, batch)) is always an order of magnitude smaller. The graded variant of the capacity factor does not remove the O(0.1) divergence — **robust to convention width**. The independent general-s observable D(s) likewise fails to shrink (0.135→0.140→0.185). Meanwhile, **support coincidence**: capacity-impossible channels and amplitude-zero channels coincide bidirectionally at s=9/11/13 (only the $(1,1,s{-}2)$ type) — the two sides **share their zeros and diverge only in weights**.

## 5. Theorem 3: the non-conventionality of interference granularity

**Internal definition (convention)**: in this series a "convention" is a choice of presentation that does not change the bookkeeping values (e.g., the diagonal and branch treatments — the granularity audit of Paper 14, Theorem 5). The choice between the channel-consistent sum (historical W: coherent across record-distinguishable final configurations) and configuration granularity (W′: square per configuration) is **not** a convention in this sense — it splits values algebraically (Fig. 3a):

| Parent | W (channel-consistent) | W′ (configuration granularity) |
|---|---|---|
| $(+3,0,0,0)$@13 | 9216 | 768 (30 configurations) |
| $(-3,0,0,0)$@13 (mirror) | **0 (exact)** | **768 = W′(A) exact** |

The mirror parent decays perfectly normally at configuration granularity; only the channel-consistent aggregate bends a healthy distribution to exact zero. Moreover (the m=4 execution of Supplement 84, propagated into this section — R1): (i) the cancellation is **not scale-free** (at m=2, $W(B)=32\neq0$ and W′ also differ) — the occurrence condition depends on m. **At m=4 (s=21) the unit of cancellation turns out to be the channel**: the mirror does not vanish as a whole ($S_+-S_-\neq0$, supporting the parity hypothesis), but **channel (5,7,9) cancels by itself** ($W_A=147456$, $W_B=0$, $W'$ both 3072) — per-channel $W=0\wedge W'>0$ occurs even at even m, so the supply of granularity divergences is **richer** than the odd-m subsequence. Also, **W′ mirror-equality holds for all channels at m≥3** (the m=2 discrepancy appears to be a small-s exception). (ii) The m-dependence is reduced by the **mirror-flip lemma** $\eta_B/\eta_A=-i(-1)^{\sigma_0}$ (394 paths, zero exceptions) to the vanishing of a **σ-parity indicator sum** $S_+=S_-$, and at configuration level is pinned to the **four-value lock** $z^-_K/z^+_K\in\{0,\infty,+1,-1\}$ (all 36 configurations at m=3: 12+12+6+6, zero exceptions; exclusive families weigh equally, 384=384; Fig. 3b) — algebra that is **the object of proof**. (iii) Sector (holonomy) dependence nearly vanishes at configuration granularity ($W'$ A≈B) — **holonomy information is localized in the relative phases between configurations**.

## 6. Theorem 4: conditioning demands the branch-channel stratum

**Terminology note**: the "A/B" of mirror pairs in this paper (parents $\pm m$) and the "A/B" holonomy classes of Paper 14, §5 are **distinct** — beware when reading both.

The joint amplitude of twins depends strongly on the relative configuration of the parent pair (including completely cancelling pairs). The three-stratum discipline (review #4): **class-constancy (s=9 sector A, 23 pairs, all 9 classes, zero exceptions) is ① true by construction** — the invariant classes are the exact orbit classification (B₄×swap), so constancy follows from the equivariance of the machinery and has no refutation power (a consistency check). **The check that could have failed (③) — the fang of this section — is this**: pairs with dot=±1, indistinguishable by intensity fringes, have **different** P(X) (0.972 vs 0.994 — they could have been equal). Hence the measure cannot be supported on strata 1–2 of the readout hierarchy 11⊂13⊂21 (Paper 13, Theorem 3) and **requires stratum 3 (branch channels)**. The principle "interference granularity = record distinguishability" holds only with the stratum made explicit. The floor (pairs with disjoint axis supports = the configuration-granularity value 0.96644) is also ③: it agreed with the granularity audit independently.

## 7. Event structure: exactly one bit

One-shot (no intermediate append → outcome = merged occupied set → splits indistinguishable → coherent across splits) versus sequential (the first append fixes the split) — the record principle translates this single bit into interference granularity. The historical counting dichotomy (batch 60/61, sequential 396/403) is its phase-free shadow; everything fits the 2×2 table (numerically complete):

| | one-shot | sequential |
|---|---|---|
| counting | 0.98361 | 0.98263 |
| amplitude | pair-dependent: mean 0.989, range [0.966, 1.000] + completely cancelling pairs (surviving both keyings) | pair-dependent: mean 0.781, range [0, 1] |

The amplitude side retains pair dependence under either event structure — conditioning (§6) is an inseparable component.

## 8. The proven cancellation theorems (with scope qualifications)

This line of inquiry has proven the exact cancellations on channel-consistent sums in the **style of involutions and indicator sums**: (i) $\Sigma\eta^2=0$ (the piecewise involution σ\*\* — completeness, no fixed points, $\eta\to\pm i\eta$; **at the exhaustive level of s=9**). (ii) A1 ($\Sigma i^{\delta_1}=0$) — the shell-level involution $x\leftrightarrow x\mp2$ (borrow ↔ lend). (iii) B2 ($Z_4$ equidistribution of sin-branch counts) — torsor residues (4-solution orbits are complete residue systems) × stage factorization; **based on the fiber structure of the 531 path set; the 333 final is not equidistributed**. (iv) The mirror-flip lemma and the indicator-sum reduction (§5). These are the **internal structure** of one aggregation candidate (channel-consistent sums), complementary to the adjudication of its status (this reduction).

## 9. Main theorem (classification of the candidate family) and the decision table

> **Theorem 6 (classification within the family — the main theorem)**: for the inventory-defined family of aggregation candidates — {counting D1/D2; amplitude channel-consistent W; configuration granularity W′} × {one-shot, sequential} × {with/without relation-class conditioning}: (i) any two candidates of the family **can** differ by total variation O(0.1) at finite $s$ (exhaustive at $s\le13$) — they are genuinely distinct measures. (ii) All candidates share support coincidence (§4) and canonical counting (Paper 14, Theorem 5). (iii) **The independently defined candidates organize without duplication into this product, and each of the three axes is independent** — every axis has an explicit example where changing it alone changes the measure (granularity: the mirror $W=0\wedge W'>0$; event structure: 0.98361 vs 0.98263, and joint-amplitude means 0.989 vs 0.781; conditioning: the dot=±1 inequality). This is a **classification theorem within the family**; no universal claim is made that candidates outside the family cannot exist (limitation per review #1; the non-trivial content of (iii) is the duplication-free organization and the independence of the axes — R2).

### Decision table (experiments computable inside the model)

| # | Experiment | What it decides | Status |
|---|---|---|---|
| 1 | Parity hypothesis: m odd ⟺ $S_+=S_-$ | the occurrence condition of mirror cancellation | **m=4 executed (Suppl. 84): supported** — $S_+-S_-\neq0$; however channel (5,7,9) cancels even at even m (the unit refines to per-channel indicator sums). **W′ mirror-equality is generic for m≥3** (all channels exact). m=5@s31 remains (the Born load-bearing point) |
| 2 | Mutual convergence of candidates at large s ($\Delta/\Delta'$ series, per granularity) | whether measure selection matters only in the deep discrete regime or everywhere | open |
| 3 | The ι construction (axis-0 branch flip of the hidden vx) | proof component of the four-value lock | **executed (Suppl. 84): complete at m=3** — configuration-preserving, σ-flipping, 384/384; multiplier exactly ±1 (192/192). **At m=2 the closure of the path set breaks (2/10)** — the structural address of "why odd m" = ι-closure. Remaining: the distribution law of the multiplier within configurations |
| 4 | Gauge independence (invariance of the measure under marking and orientations) | the falsifiable content of principle v2 | open |
| 5 | Comparison of isomorphic pairs in the unsaturated regime at s=11 | whether the dot=±2 equality is saturation or an extra symmetry | open |

## 10. The status of the Born rule and the acceptance criterion

Under this reduction, the Born rule is not "an exact target to reproduce" but an **effective law to be merged with within limits and censorship precision**. The acceptance criterion is the triple: **convergence** (the selected measure merges with Born in a suitable limit — the limit map itself is a construction task), **suppression** (in domains corresponding to verified regimes, deviations below **statistical censorship** ~1/√N — the statistical version of Paper 15's censorship quantum), and **prediction** (verifiable structure in the deep discrete regime — granularity, parity, quantization). The floor is secured: if internal principles do not decide, one aggregation rule can be installed as a postulate, placing the axiom count level with standard theory (never below). The upside is that the record principle (the branch-channel stratum) fixes the aggregation uniquely — then an analogue of the Born rule becomes a theorem.

## 11. Honest limits

1. **No selection has been made** (§1.2) — the reduction is preparation for selection; which of the three components is fixed first depends on executing the decision table.
2. The numerical series are exhaustive facts at s≤13 (partly s=11). Large-s behavior is decision-table item 2.
3. The parent-pair sweep of joint amplitudes is s=9, sector A (23 pairs). Sweeps at s≥11 and mixed-sector pairs are unexecuted.
4. The four-value lock and the weight pairing are exhaustive facts at m=3 (the proof is in progress along decision-table item 3).
5. No physical identification is made.

## Figures

- **Figure 1** (`paper16_fig1_compression.png`): two-outcome compression — joint distributions at s=9 (2 outcomes, Δ=0.00068) and s=11 (3 outcomes, Δ=0.24136)
- **Figure 2** (`paper16_fig2_divergence.png`): the divergence landscape — the Δ(s) series over all parent types, and configuration-granularity Δ′ (collapse of sector dependence)
- **Figure 3** (`paper16_fig3_granularity.png`): non-conventionality of granularity — the mirror-cancellation m-series (absent at m=2; full at m=3; per-channel (5,7,9) at m=4) and W′ mirror-equality (m≥3); the four-value lock 12+12+6+6 and the weight pairing 384=384

All figures are presentations of machine-verified values (no schematics).

## Appendix A: Verification summary (labels: ① true by construction / ② theorem + confirmation / ③ check that could have failed)

| Claim | Label | Verification | Method | Script |
|---|---|---|---|---|
| Anchor reproduction | ③ | all points | exact (rational arithmetic) | supplement73_handoff_execution.py |
| Δ(9/11/13), all parent types | ③ | table (Fig. 2) | exhaustive | same |
| Support coincidence | ③ | s=9/11/13 bidirectional | exhaustive | same |
| Graded convention width | ③ | s=9/11/13 | exhaustive | same |
| Mirror parent W=0∧W′ / m=2,4 | ③ | m=1–4 | exhaustive | supplement75/79/84 scripts |
| Class constancy | **① (consequence of equivariance — consistency check)** | 23 pairs × 9 classes | exhaustive | supplement75_granularity_execution.py |
| Fringe-class test (stratum-3 requirement) | **③ (the fang)** | dot=±1 inequality | exhaustive | suppl. 75/77 |
| 2×2 table (incl. A_seq) | ③ | complete | exhaustive / sweep | supplement75_granularity_execution.py |
| Cancellation theorems (σ\*\*, A1, B2) | ② (involution constructions = proofs) | s=9 level | exhaustive + involutions | supplement69/70/71 scripts |
| Mirror-flip lemma, indicator-sum reduction | ② (three-line arithmetic) | 394 paths, 0 exceptions | exhaustive | supplement81_mirror_lemma_check.py |
| Four-value lock, weight pairing | ③ (ι is upgrading it toward ②) | 36 configs, 0 exceptions; 384=384 | exhaustive | suppl. 83 (verification record), suppl. 84 (ι) |

---

**Acknowledgments / history**: verification was carried out under the two-party independent verification protocol of Claude Code (local machine verification) and claude.ai (independent re-computation; proposal and mutual correction of several propositions). Source material: Supplements 68–84 (June 12, 2026; suppl. 83 = the four-value-lock verification record, suppl. 84 = the m=4 and ι executions). In the course of this line, predictions were refuted (scale-freeness), explanations corrected (σ₀ constancy), and definitions reconciled (the three fringe strata) across seven exchanges between the two parties — the reliability of the reduction rests on this history of mutual correction. v0.2: 6 mandatory points of the first review. v1.0: points R1–R6 of the second review — propagation of Supplement 84 into §5 and Fig. 3, the non-trivialized Theorem 6(iii), the compression lemma, distribution ranges, terminology note.

No physical identification is made.
