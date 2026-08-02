# Lock Dynamics in an Anonymous Two-Channel Closed Wave System
## Demonstration of the Reality of the Quantization Mechanism (Arnold Tongue) and an Elimination Program for Order-Selection Mechanisms — Measured Tongue Widths, the Fixed-Point Map, and Clock Inheritance of Observational Selection

**Author:** Noriaki Kihara<br>
**ORCID:** 0009-0004-6753-4020<br>
**Date:** August 3, 2026<br>
**Version DOI:** `10.5281/zenodo.21764000`<br>
**Concept DOI:** `10.5281/zenodo.21763999`<br>
**Position:** Additional paper v1 of the two-channel exchange-scattering system in the "Wave Information Readout" series (ninth-paper series; **Part III of a trilogy**. Part I = Two-Grammar Decomposition [A]; Part II = Structure of Counting Readouts [B])

---

## Abstract

**(Zeroth-order no-go)** The standard readout of the anonymous two-channel closed wave system (phase-blind per-bin power ratios) is an exact constant of the motion (drift $\le2.2\times10^{-16}$; theorem in [A]), and neither dynamics of the coupling angle $\theta$ nor locking exists.

**(Reality of the quantization mechanism)** Under a phase-sensitive readout (the power ratio of the symmetric channel $X=A+B$, introduced as a conservative working hypothesis), $\theta$ evolves and **a genuine Arnold tongue is real**: over a finite interval of initial amplitudes $[2.518,\ \ge3.44]$ (width $\ge0.92$), the mean of the late-time winding number locks exactly to $\theta/\pi=1/3$ (machine precision over 800 collisions; $\theta$ itself oscillates on a periodic orbit — rational mean winding). Under the root convention, $(m,n)=(1,6)$, order 6. Together with the locks of orders 3, 4, and 6 obtained by the counting construction ([B], the complete menu of the crystallographic restriction), **the mechanism of quantization plateaus is real in an endogenously coupled system in which the coupling is not externally imposed**.

**(Elimination of selection mechanisms)** Locking, however, is not a mechanism of **selection**. (i) In-model measurement of tongue widths: against the width $\ge0.92$ for denominator 3, tongues with denominator $\ge4$ in the search band have width $<0.002$ — **a width ratio $>461$**. The weight of denominator 124 is negligible. (ii) Immediately below the lower edge of the tongue lies a supercritical region in which the winding number varies non-monotonically; the closest approach to the target $\rho=39/124$ is $3.3\times10^{-4}$, with no plateau detected. (iii) The attractors of the second dynamical candidate, the fixed-point-type readout (cross spectrum), are distributed over a continuous band $[0.016,0.450]$ with **no rational concentration** (0/41) — eliminated. (iv) Observational selection (repeated observation confirming recurrence to the initial state) yields fidelity exactly 1 at exact finite-order roots, directly showing that **selection is inherited by the divisor structure of the observation clock**: an observer with clock $J=8$ cannot see the order-124 recurrence ($F_8=0.0026$), while an observer with $J=248$ can ($F_{248}=1.000000000000001$). However, the fidelity landscape is smooth in the neighborhood of the roots, and finitely many observations can narrow the state down only to the recurrence **basin** (concentration ratio $1.001$ at $n=20$, recorded as a refutation of prediction P2).

**(Consequence)** Quantization (why discrete) is explained by locking; selection (why this value) is not explained by locking. With the elimination program complete, the order-124 selection problem has been transformed into a single falsifiable question: **"why is the observation clock commensurable with 248?"** This final form is taken up by the parent paper (the main paper of the ninth-paper series).

**Keywords:** Arnold tongue, mode locking, devil's staircase, supercritical circle map, return fidelity, observational selection, quantum recurrence, reproducible computation

---

## 0. Conclusion

$$
\boxed{
\begin{aligned}
&\text{Under the phase-blind readout, no }\theta\text{ dynamics exists (no-go, }2.2\times10^{-16}\text{).}\\
&\text{Under the phase-sensitive readout, an Arnold tongue is real: }\rho=1/3\text{ exact, width}\ge0.92\text{ (mechanism of quantization).}\\
&\text{By the width ratio}>461\text{, fixed-point elimination, and clock inheritance of observational selection,}\\
&\text{neither locking nor fixed points select order 124.}\\
&\text{Final form of the selection problem: "why is the observation clock commensurable with 248?"}
\end{aligned}
}
$$

## Position of This Paper (Three-Part Structure)

Part I [A] establishes the binomial identity decomposition of the flow and the two grammars (addition/product, sign, neutralization, hierarchy); Part II [B] establishes the structure of counting readouts (meta accounting rules, convention audit, Niven intersections, integer-pair readouts). This paper (Part III) treats dynamics: the reality of the mechanism of quantization plateaus, and an elimination program for candidate order-selection mechanisms. The three parts share a single reproduction package.

## 1. Research Problem

By [B], the locks reachable by the counting construction are limited to the Niven intersections (orders 3, 4, 6). Meanwhile, the root $R_{124,23}$ corresponding to the physical $\alpha$ has order 124 [2], and this selection lies outside counting. The central problem of the ninth-paper series is "what selects order 124?", and the candidates are organized into four: (a) static state structure (already refuted in [B]), (b) dynamical locking, (c) fixed-point-type attractors, (d) selection by observation (recurrence confirmation). This paper adjudicates (b), (c), and (d) by experiment.

## 2. System and Design Boundary

The system, the state construction, the standard readout, and the rotation are identical to Section 2 of [A] (scattering core and readout unchanged; target values are used only for state construction). The introduction of the phase-sensitive readout is made explicit as a conservative working hypothesis (Section 4), and the "observation" of observational selection is defined as the evaluation of the return fidelity to the initial state pair,

$$
F_J=\frac{|\langle a_0|a_J\rangle+\langle b_0|b_J\rangle|^2}{(\langle a_0|a_0\rangle+\langle b_0|b_0\rangle)^2}
$$

(no back-action on the dynamics; a pure readout).

## 3. Zeroth-Order No-Go

The standard readout is a function only of the per-bin rotation invariants $|A_k|^2+|B_k|^2$ and is exactly conserved under real orthogonal rotations (theorem in [A], Section 4). Numerical confirmation: over 400 collisions × 4 amplitudes, the maximum drift of $\theta$ is $2.2\times10^{-16}$. Hence dynamics, locking, and selection of $\theta$ do not exist at zeroth order. All of the dynamics below rests on an explicit choice of a readout that breaks this no-go.

## 4. Readout Classification (Dynamical Aspect)

The seven $\theta$ readouts were classified by their behavior under the dynamics they themselves drive (the dynamical aspect of the table in [A], Section 4). Only the diagonal sum is conserved. Within the dynamical class, the cross-spectrum readout (R4) is of **fixed-point-convergent type** (it moves substantially away from the initial value and then freezes; initially misidentified as conserved and corrected by the classification experiment — see the Addendum), while the $X$-power readout (V1) is of **oscillatory, locking type**. In what follows, V1 is used for the lock probe and R4 for the fixed-point probe.

## 5. Reality of the Arnold Tongue (Quantization Mechanism)

### 5.1 Discovery and Refinement

Sweeping the initial B amplitude under V1, there exists a finite interval in which the tail average of the late-time winding number sticks exactly to a rational number:

- A coarse scan (61 points × 400 collisions) detected a plateau, and a fine scan (13 points × 800 collisions, tail 300) confirmed it: **for amplitudes 2.6–3.3, $\overline{\theta/\pi}=1/3$ to machine precision (distance 0.0, maximum $5.6\times10^{-17}$)**
- The lower edge was fixed at 2.518 by a 39/124 micro-zoom (27 points × 1500 collisions), and the upper edge at $\ge3.44$ by an additional scan — **width $\ge0.92$**
- $\theta$ itself oscillates (tail standard deviation $\sim7\times10^{-2}$), and only the mean winding number locks rationally on a periodic orbit — the standard behavior of mode locking in circle maps

Under the root convention $\theta=\pi/2-\pi m/n$, $(m,n)=(1,6)$. Together with the counting locks of [B] (orders 3, 4, 6), the mechanism of quantization plateaus in an endogenously coupled system is established.

### 5.2 In-Model Measurement of Tongue Widths (First Negation of Selection)

Not a single tongue with denominator $\ge4$ was detected in the search band $[2.470, 2.522]$ (resolution 0.002, 1500 collisions). Therefore, as an in-model measurement,

$$
\frac{\mathrm{width}(q=3)}{\mathrm{width}(q\ge4)}>\frac{0.92}{0.002}=461
$$

and the width collapses rapidly with the denominator. **The tongue weight of denominator 124 is negligible on our own data, without even invoking external theory (the general theory of circle maps).**

### 5.3 Supercritical Diagnosis and the 39/124 Upper Bound (Second Negation)

Immediately below the lower edge of the 1/3 tongue, the winding number varies non-monotonically by $\sim10^{-2}$ between adjacent amplitude points — the symptom of a supercritical region in which tongues overlap. The closest approach to the target $\rho=39/124=0.3145161\ldots$ is $3.3\times10^{-4}$ (amplitude 2.504, 1500 collisions), with no plateau detected. A naive sweep search for high-order tongues is in principle unsuitable with this readout, and the coupling-interpolation campaign originally designed became **unnecessary for the conclusion** owing to the measured ratio of 5.2 (the history is recorded in the Addendum).

## 6. Fixed-Point Map (Elimination of the Third Candidate)

The attractor $\theta^*$ under the self-driven R4 readout was measured at 41 initial-amplitude points (400 collisions). The attractors are distributed over a continuous band $\theta^*/\pi\in[0.016, 0.450]$, with **no exact concentration on rational angles whatsoever** (0/41). Fixed-point-type attractors select no special value — the third candidate is eliminated.

## 7. Observational Selection (Fourth Candidate: The Mechanism Is Real, Selection Is Inherited by the Clock)

### 7.1 Design and Predictions

For a state family (61 amplitude-family points plus 2 exact-root states: $R=1/2$ and $R_{124,23}$, constructed by inverse search), we measure the return fidelity $F_J$ at clocks $J\in\{8,248\}$ and the survival weight $F_J^n$ after $n=20$ repeated observations. Predictions (fixed before measurement): P1 — $F=1$ only at exact roots ($R=1/2$ at both clocks by $8\mid248$; the 124 root only at $J=248$); P2 — concentration ratio onto recurrent states $>10^3$ at $n=20$.

### 7.2 Results

**P1 holds in full**: $F_8(R{=}1/2)=0.999999999999999$, $F_{248}(\text{124 root})=1.000000000000001$, and **clock selectivity** — $F_8(\text{124 root})=0.0026$. That is, **which recurrences are "visible" is determined by the divisor structure of the observation clock**: an observer with $J=8$ is blind in principle to the order-124 recurrence, while an observer with $J=248$ selects both orders 8 and 124. Selection is inherited not by the state but by the observation clock.

**P2 is refuted** (recorded as designed): the concentration ratio at $n=20$ is $1.001$. The fidelity landscape is smooth in the neighborhood of the roots, and neighboring states also survive with $F\approx1-\varepsilon$, so that **finitely many observations can select the recurrence basin, but the exact point only asymptotically**. This is a quantitative limit of observational selection, consistent with the statement in the published paper [2] that "the residual vanishes in ideal arithmetic, so the depth is formally infinite."

### 7.3 Consequence: Coordinate Transformation of the Selection Problem

The adjudication of the four candidates is complete: static structure (no — [B]), dynamical locking (no — width ratio 461), fixed points (no — continuous band), observational selection (the mechanism is real, but selection is inherited by the commensurability $P\mid J$ of the observation clock). Hence "why order 124?" has been transformed into the single question

$$
\boxed{\text{"Why is the observation clock commensurable with 248?"}}
$$

This is not a retreat but an advance: the question has moved from a search of the state space to a problem of the structure of the observer (register, resolution). This final form is taken up by the parent paper (the main paper of the ninth-paper series: the generation structure of fermions and wave-packet collapse).

## 8. Claims

**Claim 1.** Under the phase-blind readout, no dynamics of the coupling angle exists (no-go; theorem plus numerics $2.2\times10^{-16}$).

**Claim 2.** Under the phase-sensitive readout, an Arnold tongue is real: $\rho=1/3$ exact, width $\ge0.92$, $(m,n)=(1,6)$. Together with the counting locks (orders 3, 4, 6 [B]), **the mechanism of quantization (why discrete) is real in an endogenously coupled system**.

**Claim 3.** Locking is not a mechanism of selection (why this value): by in-model measurement, tongue widths collapse by a factor $>461$ between denominator 3 and denominator $\ge$4, and the weight of order 124 is negligible.

**Claim 4.** Fixed-point-type channels do not select (continuous band of attractors; rational concentration 0/41).

**Claim 5.** Observational selection is real as a mechanism (fidelity 1 at exact roots) and **inherits selection into the divisor structure of the observation clock** (direct demonstration of clock selectivity). However, finitely many observations can narrow the state down only to the recurrence basin (P2 refuted). Hence the final form of the selection problem is "why is the observation clock commensurable with 248?"

## 9. Relation to Known Systems

Mode locking, the devil's staircase, and the denominator dependence of tongue widths are established mathematics of circle maps [13], and quantization plateaus produced by phase locking exist physically as superconducting Shapiro steps [12]. The contributions of this paper are: (i) demonstrating locking in an endogenously coupled system in which the coupling is read from the state rather than being an apparatus constant; (ii) completing the elimination program for selection candidates by measuring tongue widths, fixed points, and observational selection in the same system; (iii) directly showing the inheritance of selection into the divisor structure of the observation clock.

## 10. Open Problems

1. **Commensurability of the observation clock**: the final form of "why 248?". Model the observer itself as a closure (register) within the system, and search for a mechanism by which its intrinsic period becomes commensurable with the dominant recurrence period — passed to the parent paper
2. **Derivation of the V1 readout**: the phase-sensitive readout is a working hypothesis. The conserved-readout classification ([A], Section 4) narrowed the candidate space, but there is no first principle for the readout that the dynamics selects
3. **Staircase in the subcritical region**: measurement of the complete staircase via the coupling interpolation $P_\lambda$ became unnecessary for the conclusion (5.2), but remains as a quantitative theory of tongue-width scaling

All of the above are open problems internal to the model. The overall picture of the problems that remain as connections to observational physics is organized collectively, as a separate series, in the next section.

## 11. Physical-Connection Problems Beyond the Trilogy — From the Three Directions to 3+1-Dimensional Spacetime, and from the Two Grammars to Gravity and the Coulomb Force

Since this is the final paper of the trilogy, we conclude by making explicit, in a single problem table, the range established inside the model and the range that remains as connection to observational physics.

First, we distinguish the logical hierarchy. The **internal selection problem** left by the preceding sections ("why is the observation clock commensurable with 248?") asks why a particular finite-order root is selected inside the model. The **physical connection problems** of this section ask how the internal structure thus selected is realized as the observed 3+1-dimensional spacetime, gravity, and the Coulomb force. The two are non-competing, separate series of problems.

There is no need to rediscover what has already been derived: the spontaneous emergence of the three directions (the metastabilization process reported by the earlier stage of this series and by the parent paper), the inverse-square law via harmonic closure [3], the unique differentiation of the gravity-type and charge-type two grammars ([A] Section 6.5, exclusion of the inverted assignment), and the structure of counting, finite order, locking, and the observation clock ([B] and this paper). What remains is the construction of the **realization map** that carries these internal structures onto observational physics — concretely, the following problems.

**(1) Local spatialization of the three directions.** Let $T_1,T_2,T_3$ be the generating actions of the three directions that arise at metastability, and show that they form the basis of a local three-dimensional space rather than remaining three internal modes. The zeroth-order candidate is the Gram metric $g_{ij}(\Psi)=\mathrm{Re}\langle T_i\Psi|T_j\Psi\rangle$; verify whether $\mathrm{rank}\,g=3$ and positive definiteness hold for the family of metastable states. No physical meaning may be assigned to the numbering of the three directions (anonymity); the requirement is that the space be definable as equivalent under basis rotations $SO(3)$. **The existence of three directions and the establishment of a three-dimensional space are not the same thing** — this is the core of what must be examined throughout this section.

**(2) Globalization of the local three directions and integrability.** Do the three directions obtained at each closure and at each local state form a common space? Examine the structure of $[T_i,T_j]$, the parallel transport of local bases and its path dependence (holonomy), and Frobenius-type integrability conditions, and determine whether the local three directions can be connected into a single global space.

**(3) The time direction and the construction of 3+1-dimensional spacetime.** This system has multiple time candidates: collision count, phase advance, finite-order recurrence, and the observation clock $J$. Determine which of these (or which combination) carries physical time, and verify whether, together with the three positive-definite directions, a Lorentz-type metric $ds^2=-c_{\mathrm{eff}}^2 d\tau^2+g_{ij}dx^i dx^j$, causal structure, and a finite propagation speed arise endogenously.

**(4) Identity of distances.** Derive the identity of the harmonic-closure distance, the geometric distance constructed from the three directions, and the distance read by the observer,

$$
r_{\mathrm{harmonic}}=r_{\mathrm{geometric}}=r_{\mathrm{observed}}
$$

(the completed form of the $\Delta\theta\leftrightarrow r$ dictionary of [A] Open Problem 2). A scalar $r$ alone does not suffice: one must construct the direction $\hat{\mathbf r}$ joining two closures and a gradient operation, and show that a **directed** inverse-square field corresponding to $\nabla(1/r)=-\hat{\mathbf r}/r^2$ can be defined on the emergent three directions.

**(5) Map from the gravity-type grammar to physical gravity.** Map the phase-blind, additive source readout of the magnitude term onto physical mass-energy sources, and formalize how sources deform the emergent metric and connection (magnitude-type source → $\delta g_{\mu\nu}$ or $\delta\Gamma^\mu{}_{\nu\rho}$ → motion of closures). The first verification targets are: Newton-type acceleration in the weak-field, low-velocity limit; universal free fall independent of the internal structure of the test closure; agreement of the inertial-side and gravitational-side readouts; and the endogenous derivation of the hierarchy exponent linking the coupling to the scale ratio $R/R_0$ (the scale dictionary of [A]).

**(6) Map from the charge-type grammar to the physical Coulomb force.** Read the overlap term $\mathrm{Re}\langle a|b\rangle\sim|q_A||q_B|\cos(\phi_B-\phi_A)$ as a charge-type source vertex, connect it to the mediation by harmonic closure [3] and to the response on the receiving side, and derive $\mathbf F_{AB}=K\,q_Aq_B\,\hat{\mathbf r}_{AB}/r_{AB}^2$ on the emergent space (the structure source vertex × mediator × response vertex). Items to be settled individually: whether the sign of the relative phase maps to approach or recession (Bjerknes type or electromagnetic type — the automatic-adjudication prediction of [A]); whether the neutralization of non-shared channels maps to physical electrical neutrality; the conservation law of charge-type sources and a Gauss-type flux law; the path by which $1-R_\ast=\sin^2(23\pi/124)$ maps, unit conventions included, to the observed elementary charge; and whether the introduction of time-dependent phases produces magnetic and radiative components beyond the static Coulomb field.

**(7) A spacetime dictionary common to the two interactions.** Separate distances and separate spaces must not be assigned to the gravity type and the charge type. On the same emergent metric, the same $r$, and the same $\tau$, the magnitude-type grammar → gravitational response and the overlap-type grammar → charge response must hold **simultaneously**, and this correspondence must be invariant or covariant under basis rotations of the three directions, A/B channel exchange, renaming of harmonics, changes of readout convention, and the observer's choice of coordinates.

**(8) Falsification conditions for the physical connection map.** The connection is not to be certified on grounds of similarity. If any of the following holds, the identification with physical gravity and the Coulomb force is rejected: no stable rank-3 metric can be constructed from the three directions / the local three directions are not globally integrable / the harmonic-closure distance and the emergent geometric distance do not agree / magnitude-type sources produce no universal geometric response / the phase sign of the overlap type does not map to consistent attraction and repulsion / the same charge magnitude does not appear on both the source side and the response side / the two grammars do not act simultaneously on the same 3+1-dimensional spacetime.

The above are not problems that return the internal results of this trilogy to an unsettled state; they are the next-stage connection problems of mapping the derived pre-geometric structures onto observables. The current position of this series is organized as follows:

$$
\boxed{
\begin{aligned}
&\text{Spontaneous emergence of the three directions} && \text{derived (earlier stage of this series)}\\
&\text{Inverse-square mediation} && \text{derived [3]}\\
&\text{Gravity-type and charge-type two grammars (unique assignment)} && \text{derived and verified [A]}\\
&\text{Structure of counting, rational quantization, and the observation clock} && \text{derived, limits verified ([B], this paper)}\\
&\text{Three directions}\to 3+1\text{-dimensional spacetime} && \text{unconnected (Problems 1–4)}\\
&\text{Two grammars}\to\text{physical gravity and Coulomb force} && \text{unconnected (Problems 4–7)}
\end{aligned}
}
$$

The conclusion of the trilogy closes in three steps: **The two grammars differentiate spontaneously, and the assignment is unique ([A]). Their readout carries the structure of counting, finite order, locking, and the observation clock ([B], this paper). The realization map to physical spacetime and physical forces is taken up, as the eight problems above, by the next stage — the parent paper and its successors.**

## 12. Reproducibility

The runners for all experiments, all measurement CSVs, and the figures are committed, following the same discipline as [A] (SHA-verified core references, pre-measurement fixing of predictions, anchor reproduction, and preservation of refutations and corrections).

| Experiment | Folder | Commit |
|---|---|---|
| Mode-locking probe (Z₆ discovery, no-go) | `mode_locking_probe_pre_v1` | 56453653 |
| 1/3-tongue fine scan | same as above, `run_one_third_tongue_fine_scan_v1` | 15ee20cb |
| Readout classification + staircase + 39/124 zoom | `tongue_spectroscopy_pre_v1` | 913756fe, 15ee20cb |
| Tongue-width measurement (upper edge, width ratio 461) | `tongue_width_measurement_v1` | 47e2da02 |
| Fixed-point map | `partial_sharing_fixed_point_v1` (E7 part) | 47e2da02 |
| Observational selection (clock selectivity) | `observation_selection_v1` | 47e2da02 |

---

# References

## Self-citations

1. Noriaki Kihara, "Anonymous Equal-Amplitude Composite-Wave Model: Basic Axiom System v6," Version DOI: `10.5281/zenodo.21465984`, Concept DOI: `10.5281/zenodo.21315735`, 2026.
2. Noriaki Kihara, "Discovery of Finite-Order Resonance in Iterated Exchange Scattering — Identifying the Origin of Peaks near Fine-Structure-Constant Values 137 and 128 with a Reproducible Wave-Packet Model," Version DOI: `10.5281/zenodo.21421367`, Concept DOI: `10.5281/zenodo.21421366`, 2026.
3. Noriaki Kihara, "Future Phase-Position Acceleration Map and the Inverse-Square Law via Harmonic Closure in an AB Two-Body Closed Phase System v4," Version DOI: `10.5281/zenodo.21468270`, Concept DOI: `10.5281/zenodo.21441081`, 2026. (Precedent for the three-step procedure: zeroth-order no-go → conservative working hypothesis → independent verification)
4. Noriaki Kihara, "Two-Grammar Decomposition of Interaction in an Anonymous Two-Channel Closed Wave System" (Part I of the trilogy), Version DOI: `10.5281/zenodo.21763996`, Concept DOI: `10.5281/zenodo.21763995`, 2026.
5. Noriaki Kihara, "Structure of Counting Readouts in an Anonymous Two-Channel Closed Wave System" (Part II of the trilogy), Version DOI: `10.5281/zenodo.21763998`, Concept DOI: `10.5281/zenodo.21763997`, 2026.
6. Noriaki Kihara, "Preliminary Summary of Acceleration-Basis and Localization Exchange in Exchange-Interference Scattering-Matrix Fermion-Like Collisions v1," Version DOI: `10.5281/zenodo.21333768`, 2026. (Code provenance of the engine)
7. Noriaki Kihara, "Selection of the State-Exchange Weight G_R = 1 − R and Candidate Correspondence with the Fine-Structure Constant: A Numerical Experiment v1," Version DOI: `10.5281/zenodo.21396761`, 2026. (Code provenance of System A)

## External References

12. S. Shapiro, *Phys. Rev. Lett.* **11**, 80 (1963). DOI: `10.1103/PhysRevLett.11.80`.
13. M. H. Jensen, P. Bak, and T. Bohr, *Phys. Rev. Lett.* **50**, 1637 (1983). DOI: `10.1103/PhysRevLett.50.1637`.

---

**Addendum (record of corrections and history)** (i) The cross-spectrum readout R4 was misidentified as "a second conserved readout" → corrected after the classification experiment showed it to be of fixed-point-convergent type. (ii) The concentration prediction P2 of observational selection ($>10^3$ at $n=20$) was refuted — recorded in the main text as a principled limit of finitely many observations due to the smooth fidelity landscape. (iii) The originally designed subcritical-staircase campaign via coupling interpolation ($P_\lambda$) became unnecessary for the conclusion because the in-model measurement of tongue widths (width ratio $>461$) gave the same conclusion more directly. The decision history is recorded together with the judgment.
