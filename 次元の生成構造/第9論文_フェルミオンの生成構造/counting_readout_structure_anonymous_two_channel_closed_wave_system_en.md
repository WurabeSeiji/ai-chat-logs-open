# Structure of Counting Readouts in an Anonymous Two-Channel Closed Wave System
## Exact Rationality of Equal-Weight Points, the Meta Counting Rule, the Convention Audit, Complete Verification of the Niven Intersection Points, and Exact Lock States Constructed by Counting Alone

**Author:** Noriaki Kihara<br>
**ORCID:** 0009-0004-6753-4020<br>
**Date:** August 3, 2026<br>
**Version DOI:** `10.5281/zenodo.21763998`<br>
**Concept DOI:** `10.5281/zenodo.21763997`<br>
**Position:** "Wave Information Readout" series, two-channel exchange scattering system, supplementary paper v1 (ninth-paper series; **Part II of a trilogy**. Part I = Two-Grammar Decomposition [4], Part III = Lock Dynamics [5])

---

## Abstract

**(Decomposition theorem)** The readout of an anonymous two-channel closed wave system decomposes exactly into "counting + amplitude deformation." At the point where the harmonics are given equal weights (the equal-weight point), the readout degenerates to a counting of the existence set of harmonics and becomes an **exact rational number with a small denominator** (examples: $61/128$ for odd 1..63, $29/64$ for odd 1..31, $15/31$ when 3 is removed, exactly $0$ for the even control), and its value can be predicted solely from the measured bin table of each harmonic and the accounting of the mask (**meta counting rule**: error $\le8\times10^{-16}$ across all 20/20 conditions of 7 held-out supports × 5 convention variants, with the predictions fixed as in-code constants before measurement). Amplitude distributions, on the other hand, contaminate the readout continuously (spread $0.407$ for an identical existence set). **The only carrier of structural information is the existence set of harmonics.**

**(Convention audit)** This structure survives across 5 readout convention variants (changes of the mask definition). What does not survive — the concrete form of the counting rule, the concrete values of the equal-weight points, and the role assignment of odd/even — is bookkeeping. In particular, under the parity-inverted mask, the support that carries the lock moves to **30 even harmonics**: "odd harmonics = fermionic" is a convention label, and the structure (counting → lock) is independent of the label.

**(Complete verification of the Niven intersection points)** The intersection of the rational numbers $R$ reachable by counting with rational-angle locking $\theta/\pi\in\mathbb{Q}$ is, by Niven's theorem, strictly limited to the 5 points $R\in\{0,\,1/4,\,1/2,\,3/4,\,1\}$. At every reachable intersection point we verified the direct link "counting construction (equal weights, **zero amplitude search**) → exact lock," with all predictions fixed before measurement, and scored **6/6 hits**: the mixed-parity support $\{2,5\}$ and its kin land exactly on $R=1/4$ with period 12 and $(m,n)=(1,3)$; the $\{1,3\}$-free support gives $R=1/2$, period 8, $(1,4)$; odd 5..63 under the variant mask gives $R=3/4$, period 6, $(1,6)$; and $R=0$ is invariant (all errors $\le3.9\times10^{-16}$, recurrence residuals $\sim10^{-15}$). **The realized orders $n=3,4,6$ coincide exactly with the complete menu of nontrivial orders of the two-dimensional crystallographic restriction** — Niven's theorem and the crystallographic restriction are two faces of the same rational-cosine constraint.

**(Where the integers live, and the limit)** The integer pair $(m,n)$ can be read out from the order and winding number of the evolution (measured: exact reconstruction in all 3 state families). By contrast, a static state decomposition carries no winding number: winding-number spectroscopy on a cyclic register showed that the state spectrum coincides completely with mapping artifacts and that the occupancy of the target winding number is at random level (a refutation experiment). Hence the generation rule for integer pairs **outside** the 5 Niven points (including the physical $\,(23,124)$) can only be sought in dynamics, not in static structure — that verdict is delivered by Part III [5].

**Keywords:** counting readout, rational degeneration, meta counting rule, convention audit, Niven's theorem, crystallographic restriction, exact recurrence, reproducible computation

---

## 0. Conclusion

$$
\boxed{
\begin{aligned}
&\text{Readout = counting + amplitude deformation. Only the existence set carries structure (exact rationals at equal-weight points, }\le8\times10^{-16}\text{ on 20/20).}\\
&\text{Role labels (odd = fermionic, etc.) are convention; the structure is label-independent (under parity inversion the even side locks).}\\
&\text{The direct counting}\to\text{lock link scores full hits at the 5 Niven intersection points, with realized orders }\{3,4,6\}\text{ = the crystallographic menu.}\\
&\text{Integer pairs can be read only from dynamics — static decomposition carries no winding number (refuted).}
\end{aligned}
}
$$

## Position of This Paper (Three-Part Structure)

Part I [4] treats the binary identity decomposition of the flow and the two grammars; this paper (Part II) treats the structural layer of the readout — what is counting and what is appearance; Part III [5] treats the dynamics (the reality of locking and the limits of selection). The three parts share a single reproduction package.

## 1. Research Question

How much of a readout is the structure of the system, and where does the product of the way of reading begin? From where can the integers be read? [4] established the compositional grammar of readouts, but the status of the readout **value** itself — the boundary between the structural layer (independent of convention) and the bookkeeping layer (a product of convention) — remained unclassified. This paper completes that classification experimentally and verifies the full range of lock states that counting can build (the Niven intersection points).

## 2. System Definition and Design Boundaries

The system, state construction, standard readout, and rotation are identical to Section 2 of [4]. In all experiments the scattering core and the readout formula are unmodified, target values are used only for state construction, predictions are fixed as in-code constants before measurement, and each runner begins with a reproduction of known results (an anchor).

## 3. Experiment I: Two-Memory Battery (Discovery of the Decomposition)

**Condition A (identical existence set, different amplitude distributions, same pair norm)**: fixing the support to odd 1..63 and varying the amplitude distribution over 5 types (equal weights, linear decay, $1/k$, exponential decay, random), the readout $R$ moves across $0.069$–$0.477$ (**spread $0.407$**). A power-ratio-type readout is strongly contaminated by amplitude (the visible layer).

**Condition B (same norm, different existence sets, equal weights)**: at the equal-weight point the readout degenerates to an exact rational number.

| Support | $R$ (exact value) |
|---|---|
| Odd 1..63 (32 harmonics) | $61/128$ |
| Odd 1..31 (16 harmonics) | $29/64$ |
| Odd step 4 (16 harmonics) | $15/32$ |
| Odd with 3 removed (31 harmonics) | $15/31$ |
| Even control | $0$ (exact) |

The denominator tracks the support size, and the exact zero of the even control is a structural readout of parity. Control test: the equal-weight single-harmonic sum agrees exactly with the template construction ($1.7\times10^{-16}$). Invariance under evolution (100 collisions) has drift $\le3.3\times10^{-16}$.

**Consequence: readout = counting readout + amplitude deformation. The equal-weight point is the structural representative of each support class.**

## 4. Experiment II: The Meta Counting Rule and Held-Out Prediction Verification

The counting rule decoded from the 5 points of Condition B (its concrete form under the canonical mask),

$$
R_{\mathrm{eq}}(S)=\frac{2|S|-c(S)}{4|S|},\qquad
c:\{k{=}1\to2,\ k{=}3\to1,\ \text{even }k\to2,\ \text{odd }k{\ge}5\to0\}
$$

was verified as a **prediction fixed before measurement** on 7 held-out supports ($\{1,3\}$, $\{5\}$, odd 5..63, $\{1,5,7\}$, $\{3,7\}$, odd 1..15, even-mixed $\{1,2,5\}$), and every one was a hit (error $\le5\times10^{-16}$, including the generalized $c$ rule for the even-mixed case).

However, the concrete form $c$ is bookkeeping of the canonical mask (Section 5). The convention-independent core is the **meta counting rule** — "the $R$ of an equal-weight point can be predicted solely from the measured bin table of each harmonic alone and the accounting of the mask" — and this survives under all conventions.

## 5. Experiment III: Convention Audit

We swept the fermion mask over 5 variants (canonical even$\ge$4 / even$\ge$6 / even$\ge$2 / any$\ge$4 / parity-inverted odd$\ge$3) and measured all 20 combinations with 4 state classes (a control test confirmed that the variant readout agrees exactly with the original implementation in the canonical setting, difference $0.0$).

**Physics-grade (survives all conventions)**: the meta counting rule (20/20, $\le8\times10^{-16}$), the rationality of the equal-weight points ($61/128$, $59/128$, $95/128$, $61/124$, $3/4$, $1/4$, etc.), the amplitude dependence (excluding degenerate classes where the readout is identically zero), and the counting → lock bridgehead. **In particular, under the parity-inverted mask, 30 even harmonics lock** (period 8, $(1,4)$) — the convention relativity of role assignment.

**Bookkeeping-grade (products of convention)**: the concrete form of the $c$ rule, the concrete values of the equal-weight points, the role assignment of odd/even, and the existence condition of the bridgehead (under even$\ge$2, A is no longer purely bosonic and the counting lock disappears).

All claims are stated in meta form (physics-grade). All audit measurements are collected in a figure (`audit_figure_v1`).

## 6. Experiment IV: Complete Verification of the Niven Intersection Points

### 6.1 Intersection Theorem

The equal-weight points of counting give rational $R$. On the other hand, exact locking requires $\theta/\pi\in\mathbb{Q}$ (a rational angle). Compatibility of $R=\cos^2\theta$ is, by Niven's theorem (if $\theta/\pi\in\mathbb{Q}$ and $\cos\theta\in\mathbb{Q}$, then $\cos\theta\in\{0,\pm\tfrac12,\pm1\}$ [8]), limited via the rationality of $\cos2\theta$ to the 5 points

$$
R\in\{0,\ 1/4,\ 1/2,\ 3/4,\ 1\}
$$

The corresponding orders are $n\in\{3,4,6\}$ (the nontrivial points) — **the complete menu of nontrivial orders of the two-dimensional crystallographic restriction** (two faces of the same rational-cosine constraint).

### 6.2 Full Prediction Verification (6/6 Hits)

All predictions were fixed as in-code constants before measurement, and every reachable intersection point was verified:

| Case | Support | Readout | $R$ (prediction / error) | Period (prediction / measured) | $(m,n)$ (prediction / measured) |
|---|---|---|---|---|---|
| Anchor | Odd 5..63 | Canonical | $1/2$ / $3.9\times10^{-16}$ | 8 / **8** | $(1,4)$ / **$(1,4)$** |
| P1a–c | Mixed parity $\{2,5\}$ and 2 kin (3 types) | Canonical | $1/4$ / $\le3.3\times10^{-16}$ | 12 / **12** | $(1,3)$ / **$(1,3)$** |
| P2 | Odd 5..63 | even$\ge$2 variant | $3/4$ / $3.3\times10^{-16}$ | 6 / **6** | $(1,6)$ / **$(1,6)$** |
| P3 | $\{1\}$ | Canonical | $0$ / $0.0$ | Invariant ✓ | — |

All of these are pure counting constructions with equal weights and **zero amplitude search** ($R=1$ is unreachable by counting because the A fundamental wave is bosonic under all masks. Period-4 recurrence in the amplitude-tuned family has been confirmed in the existing R-sweep experiments — we record this honestly). The mixed-parity supports (the P1 series) are a state class opened for the first time by this verification.

## 7. Experiment V: A Reader for the Integer Pair $(m,n)$

For an exact lock state, we measure the order (the collision count $P$ of the first exact recurrence) and the winding number (the accumulated rotation angle over one period $/2\pi$), and obtain $(m,n)$ by rational reconstruction of $x=1/2-\bar\theta/\pi$. In all lock state families the result agreed exactly with the prediction: $(1,3)$, $(1,4)$, $(1,6)$ (winding-number error $\le6\times10^{-16}$). Integer pairs can be read out from the evolution.

## 8. Control and Refutation: Static Decomposition Carries No Winding Number

On a cyclic register of non-common harmonics 31 × 4 phases = 124 slots, we measured the winding-number spectrum of the antisymmetric component of the state under 4 canonical orderings (winding-number spectroscopy). The spectrum of the equal-weight state coincided completely with the "mapping-only control" under all orderings, every observed winding structure ($k=31,1,93$) was an artifact of the slot mapping, and the occupancy of the target winding number $k=23$ was at random level under all conditions. **A static initial-state decomposition carries no winding number** (recorded as a refutation experiment).

Therefore, the generation rule for integer pairs outside the 5 Niven points — including the physical $(23,124)$ [2] — can only be sought in dynamics. The verdict on the dynamics side (the reality of locking and the limits of selection) is delivered by Part III [5].

## 9. Claims

**Claim 1 (Decomposition theorem).** The readout decomposes exactly into counting (the structural layer: a function of the existence set, exactly rational at equal-weight points, predictable by the meta counting rule) and amplitude deformation (the visible layer: spread 0.407). The only carrier of structural information is the existence set of harmonics.

**Claim 2 (Convention relativity).** The structural layer survives across 5 readout convention variants, and the role labels (odd = fermionic, etc.) are convention. Structure alone remains without labels — the readout-convention version of the anonymity principle of [4].

**Claim 3 (Niven intersection points).** The direct link counting → exact lock holds at the 5 Niven intersection points, and only there. The realized orders $\{3,4,6\}$ are the complete menu of the crystallographic restriction. All predictions scored 6/6 hits under pre-measurement fixing.

**Claim 4 (Where the integers live).** The integer pair $(m,n)$ can be read out from the evolution (order and winding number) and cannot be read out from static decomposition (refuted). The generation rule for general $(m,n)$ is restricted to the dynamical type.

## 10. Open Problems

1. **A closed-form general theory of the counting rule**: generalization of the meta counting rule (bin table × mask) to arbitrary state constructions and arbitrary carrier structures
2. **Outside the 5 Niven points**: dynamical generation of $(23,124)$-type integer pairs — already transformed, by the elimination argument of Part III [5], into the question of "the commutativity of the observation clock"
3. **The $R=1$ intersection point**: counting reachability depends on the construction of the A anchor. A general theory including the anchor degree of freedom

## 11. Reproducibility

We follow the same discipline as [4]. All runners, all measurement CSVs, and figures are committed:

| Experiment | Folder | Commit |
|---|---|---|
| Two-memory battery | `two_memory_readout_battery_pre_v1` | 5444fbdb |
| Counting rule verification + lock link + $(m,n)$ readout | `counting_rule_lock_bridge_pre_v1` | 8dff131e |
| Convention audit | `convention_audit_pre_v1` | 5a8f4b0a |
| Audit figure | same folder, `run_audit_figure_v1` | 47e2da02 |
| Niven intersection-point lock | `niven_points_lock_pre_v1` | 8f99d5c6 |
| Winding-number spectroscopy (refutation) | `winding_spectroscopy_pre_v1` | 641048ad |

---

# References

## Self-citations

1. Noriaki Kihara, "Anonymous Equal-Amplitude Composite-Wave Model: Basic Axiom System v6," Version DOI: `10.5281/zenodo.21465984`, Concept DOI: `10.5281/zenodo.21315735`, 2026.
2. Noriaki Kihara, "Discovery of Finite-Order Resonance in Iterated Exchange Scattering," Version DOI: `10.5281/zenodo.21421367`, Concept DOI: `10.5281/zenodo.21421366`, 2026.
3. Noriaki Kihara, "Future Phase-Position Acceleration Map and the Inverse-Square Law via Harmonic Closure in an AB Two-Body Closed Phase System v4," Version DOI: `10.5281/zenodo.21468270`, 2026. (Precedent for the three-step procedure)
4. Noriaki Kihara, "Two-Grammar Decomposition of Interaction in an Anonymous Two-Channel Closed Wave System" (Part I of the trilogy), Version DOI: `10.5281/zenodo.21763996`, Concept DOI: `10.5281/zenodo.21763995`, 2026.
5. Noriaki Kihara, "Lock Dynamics in an Anonymous Two-Channel Closed Wave System" (Part III of the trilogy), Version DOI: `10.5281/zenodo.21764000`, Concept DOI: `10.5281/zenodo.21763999`, 2026.
6. Noriaki Kihara, "Exchange-Interference Scattering Matrix … Summary of Preliminary Locality-Exchange Experiments v1," Version DOI: `10.5281/zenodo.21333768`, 2026. (Code provenance of the engine)
7. Noriaki Kihara, "Numerical Experiments on Exchange-Scattering Coefficient R Concentration and Fine-Structure-Constant Correspondence Candidates v1," Version DOI: `10.5281/zenodo.21396761`, 2026. (Code provenance of System A)

## External References

8. I. Niven, *Irrational Numbers*, Carus Mathematical Monographs No. 11, MAA (1956). (The theorem on the rationality of cosines of rational angles. The common root of the intersection set and the crystallographic restriction.)

---

**Addendum (Record of a Correction)** In the course of this research, an error was once recorded that took the intersection of the counting rational family and the rational-angle lock family to be "$R=1/2$ only" (both the occurrence and the correction are appended and preserved in the experiment folder README). The correction to the 5 points by Niven's theorem is Section 6 of this paper, and the verification of the remaining 2 points arising from the correction (P1 and P2) was completed with all predictions hitting. We record, together with its history, how the error led to the discovery of a new state class (the mixed-parity supports).
