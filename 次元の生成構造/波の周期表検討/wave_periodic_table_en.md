# The Periodic Table of Waves: A Hypothesis for Particle Classification by Winding Addresses and Observation Clocks

**Author**: Noriaki Kihara
**Date**: August 6, 2026
**Type**: Hypothesis-proposal paper (supported by numerical experiments; not a proof paper)
**Version DOI**: 10.5281/zenodo.21822359
**Concept DOI**: 10.5281/zenodo.21822358

---

## Abstract

By what are particles classified? From numerical experiments on wave dynamics that assume no particle concept whatsoever (two published systems, used as-is), this paper proposes the **"periodic table of waves"** — a hypothesis for particle classification by winding addresses and observation clocks. The starting point was a sign-vector thought experiment [1] classifying particles by which of the six axes (x, y, z, t, R, Q) carries their winding. This paper rebuilds that classification from wave forms alone, by measurement, arriving at five pillars. **Pillar 1 (periodic law)**: the foundation of classification is a measured universal clock — the collective clock ω = π/72 per step holds across the entire range N = 4 to 144 (long-window ±0.1%) — and species are classified by rational addresses on a single U^144 clock lattice and their relation to the clock. **Pillar 2 (charge rule)**: charge is Q_em = m/3 (m the raw winding number), read by the dominant observation clock (denominator 3 — based on the previously measured tongue-width ratio of 461) as "3 windings = 1 elementary charge." It follows that the u-type m = +2 reads +2/3, the d-type m = −1 reads −1/3, the electron m = −3 reads −1, and **confinement and fractional charge become two faces of the single fact of "mod-3 readability"** (the bookkeeping checks for the proton uud = +3 → +1 and the neutron udd = 0 pass automatically). In the numerical experiments, quark-type addresses remain exactly unreadable in isolation (readable fraction 0.0000 persists), while in the sea they become readable only in integer charge units, to the extent that hadronization (formation of m ≡ 0 composites) proceeds. **Pillar 3 (statistics rule)**: the Fermi/Bose distinction lives not in spatial rotation but in the **double cover of the clock**. The state phase is quantized to {0, π} at every observational recurrence (measured to 3 digits); fermionic species (odd harmonics) visit both sheets (covering degree 2), bosonic species do not (degree 1). A pure m = 0 fermionic-band species also shows covering degree 2, establishing the row of the neutral fermion (neutrino). **Pillar 4 (mass rule)**: mass is not an attribute of the address but a relational quantity with the sea (vacuum) — isolated addresses are all unitarily equivalent under winding shifts (measured), and in-sea properties depend only on the divisor class gcd(m, 16) (odd addresses agree to 5 digits). **Pillar 5 (the table)**: the table is a double signboard — an isolated stable-element table and an in-sea lifetime table — on which a 62-species Standard-Model assignment is proposed. Under this hypothesis, the problem left open by the previous paper [5] — "charge does not stabilize at ±1" — turns from a contradiction into a prediction, and the 1/3-charge problem of quarks connects to confinement through a single mechanism. This paper claims no proof. What it claims is that seven independently measured structures connect, under a single classification hypothesis, simultaneously to the particle structure predicted by the Standard Model, to the open problem of the previous paper, and to the sign-vector thought experiment. All experiments are published as deterministic scripts, and the entire process — including 13 refutations — is reproducible. The closing section fixes the requirements that the universal interaction function of the next paper must satisfy in order to verify this classification hypothesis dynamically.

---

## 1. Introduction — why a "periodic table," and why a hypothesis proposal

### 1.1 Starting point: the sign-vector thought experiment

A preceding thought experiment [1] proposed classifying particles by which of the six axes (x, y, z, t, R, Q) carries their winding — photon γ = (0,0,0,0,R,0), W± = (0,0,0,±t,0,0), gluon g = (0,0,0,0,0,Q), graviton G = (0,0,0,±t,±R,0). The classification is attractive but remained at the level of a thought experiment. Can particle species be classified as **measurements** of wave dynamics? That is the question of this paper.

### 1.2 Central proposition

All pillars of this paper are bound by one proposition:

> **Particle species are classified not as fixed attributes of the waves themselves, but by how an observation clock reads winding states.** Charge, confinement, statistics, mass, and lifetime are not independent labels; they decompose into **state-side address × observation clock × relation to the sea**.

Charge quantization exists neither in the state space alone nor in the measuring apparatus alone — it exists in **the relation between the cyclic structure of the state and the divisor structure of the clock** (§4). Statistics lives in the double cover of the clock, not in spatial rotation (§5), and mass is differentiated by the sea (§6). This is the most important theoretical proposal of this paper.

### 1.3 Methodological declaration: no particles assumed

This paper introduces no new dynamics. It uses two published systems — the N-body relational-wave dynamics [2] (used for measurements on the three-directional condensate) and the closed-form exact solution of the two-body universal inelastic map [3] (used for all charge, statistics, and confinement experiments) — both built on the anonymity principle: "wave forms are the only input, and there is no special-casing inside." The concepts of particle, charge, and spin do not enter the dynamics. We measure how they arise **as readouts**.

### 1.4 Declaration of character, and historical lessons

This is a **hypothesis-proposal paper**. Attempts to classify particles have a history of both success and failure. Gell-Mann's Eightfold Way [9] succeeded, predicting the Ω⁻ from SU(3) classification. Kelvin's vortex atoms [8] — atoms as knots of vortices, a topological classification — were beautiful, but the dynamics was wrong. Preon models [14] attempted 62-species-scale full assignments and faded. **A beautiful classification does not guarantee correct dynamics.** That this paper remains a hypothesis proposal rather than a proof, strictly separates readout-level from dynamics-level claims (§2.3), and states explicit falsification conditions (§11), is due to this lesson. The value of the hypothesis should be measured by the breadth of its connections and the clarity of its falsification conditions.

## 2. Methods and conventions

### 2.1 Dynamics (read-only)

- **N-body system** [2]: closure dynamics of N relational waves (edges of the complete graph), forming a three-directional metastable condensate. Used for the clock-universality scan (§3).
- **Two-body system** [3]: the closed-form exact solution of the pointwise interaction δa = −2R·Im(b̄a)·b of a two-channel wave (a, b). States live on a (χ, η) grid, and the η winding number m is the object of charge readout. Used for all experiments of §4–6.

Both are deterministic; every script and result can be re-executed standalone according to the attached inventory (実験一覧_v1.md).

### 2.2 Machine-verified conventions

To avoid convention collisions, two points were machine-checked (attachment #25): (i) bundle-construction bin numbers and χ frequencies correspond with a ±1 shift and inverted parity (even-bin bundles → odd frequencies and vice versa, power fraction 1.0000). The text standardizes on the **χ-frequency representation**. Fermionic = even χ frequency (≥4) = odd harmonics [7]. (ii) Bundles carry exactly one unit of η winding m = +1 by construction (η occupancy (1, 1.000)). A species of arbitrary winding m is constructed by the winding shift e^{i(m−1)η}, not by projection.

### 2.3 The most important convention: readout level vs. dynamics level

The claims of this paper are restricted to the **readout level** — "particle classification arises as a relation between observation clock and address." The **dynamics level** — "62 species are generated, decay, and interact as the spectrum of a single dynamics" — is outside the scope of this paper and is delegated to the next paper (the universal interaction function). The closing section (§13) fixes its requirement table on the basis of measurements.

## 3. Pillar 1: the foundation of the periodic law — the universal clock (measured)

The collective clock of the N-body condensate was scanned from N = 4 to 144 (28 points).

**Measurement**: ω_clock = π/72 per step holds across the entire range — ratio 0.988 ± 0.028 in short windows (T = 4000), and 0.9992–1.0013 (**±0.1%**) in long windows (T = 42000). One clock revolution = 144 steps, independent of the effective energy (N).

![Figure P1](figs/fig_p1_clock_universality_v1.png)

**Figure P1**: Clock universality: ω/(π/72) = 1 for all N.

The long-window census further showed that in the low-energy sea (N ≤ 16) exactly **one** species is long-time stationary: the massless ground species locked to the clock, ρ = 1/1 (|ρ−1| < 9×10⁻⁴; mass degree 10⁻¹¹). The massive sidebands seen in short windows are all finite-lifetime transients, entrained into the clock (Figure P2). The refutation process (short-window apparent addresses and mass laws vanishing in long windows) is fully recorded in §12.

![Figure P2](figs/fig_p2_stable_vs_resonance_v1.png)

**Figure P2**: Short-window sidebands (resonances) collapse to 1/1 in long windows.

**Hypothesis (Pillar 1)**: the classifier of particle species is the rational address of U^n = I, and the foundation is the single U^144 clock lattice. Species are classified by "address (winding) × relation to the clock."

## 4. Pillar 2: the charge rule — Q_em = m/3 and confinement = mod-3 readability

### 4.1 Multivaluedness of charge (state-side measurement)

The charged species (winding q = +1) is longer-lived than neutral transients by orders of magnitude (τ ≈ 1.3×10⁴ collisions; 85% retention), but not permanent. Its decay is not disappearance but transport to the partner (+3) by the walk of the sum rule m* = 2m_B − m_s [4] (Figure P3). Moreover, isolated pure winding species **replicate themselves exactly for any m** (retention 1.000; τ numerically infinite), and even the sea-free mixture {+1, +3} does not leak — **the walk is driven only by the m = 0 sea** (leak destinations {0, −1, +2}/{±4} agree completely with the sum-rule products). The previous paper's [5] §9.5 finding — "charge does not stabilize at ±1" — is a direct consequence of this sea-driven walk.

![Figure P3](figs/fig_p3_charged_lifetime_walk_v1.png)

**Figure P3**: Quasi-stability of the charged species and the sum-rule walk (+1 → +3).

### 4.2 The conservation law of winding charge is cyclic

The true conservation law of winding charge is not integer but **cyclic, mod ne (the η register order, 16 in this system)**. The pointwise interaction is a cyclic convolution of η spectra (in the continuum limit δQ = R∮∂_η(s²) = 0), and conservation of the integer charge Q_wind is merely a corollary for in-band states (exact conservation CV ≈ 0 measured in the charged census construction, Figure P6a). The apparent violation is the fold-back of the doubling walk reaching Nyquist (37% edge-power accumulation; correlation −0.93, Figure P6b). On the finite register the doubling ladder neutralizes in 4 steps: 1 → 2 → 4 → 8 → 0.

![Figure P6](figs/fig_p6_cyclic_conservation_v1.png)

**Figure P6**: Cyclic conservation and Nyquist fold-back.

### 4.3 Readout rectification (readout-side measurement)

Against all products of the sea-driven walk (the doubling-and-negation ladder 2^b·(±1)), **only the observation clock J = 3 reads the entire charged content as |q| = 1** (concentration 1.000 — in the +1-seeded universe and in the +2-seeded universe alike). J = 4 erases charge; J = 5, 6 scatter it into multiple values (Figure P4). The algebraic basis is 2^b mod 3 = ±1. That the observation clock has denominator 3 rests on prior measurements [6][7] — the locking tongue of denominator 3 dominates with a width ratio > 461. The value of the elementary charge itself connects to the published measurement [6] as the rational address sin²(23π/124) = √(4πα) (agreement 16 parts in 10⁸).

![Figure P4](figs/fig_p4_rectification_v1.png)

**Figure P4**: Readout rectification: only the denominator-3 observation clock reads all charged content as |q| = 1.

The bookkeeping also closes: in conserving constructions, ΔQ3 = −3ΔW holds to a precision of 7×10⁻¹⁰ (Q3 = net charge readable mod 3; W = the ledger of charge hidden in mod-3-neutral composites). **Readable charge never disappears — its decrease equals, exactly, the transport into composites that the denominator-3 clock reads as neutral** (Figure P5).

![Figure P5](figs/fig_p5_ledger_v1.png)

**Figure P5**: The bookkeeping identity ΔQ3 = −3ΔW.

### 4.4 Hypothesis (Pillar 2): Q_em = m/3 and confinement

We propose the charge rule that binds the above: **charge is Q_em = m/3, and the dominant observation clock (denominator 3) reads "3 raw windings = 1 elementary charge."** Consequences:

- u-type m = +2 → +2/3; d-type m = −1 → −1/3; electron m = −3 → −1; ν m = 0 → 0.
- **Confinement = mod-3 readability**: the observer's space is the Z₃ quotient of the η circle, and only species single-valued on it (m ≡ 0 mod 3) can be read out freely. Quarks (m ≢ 0) are unreadable alone — **the 1/3 charge and confinement are two faces of the same fact** and require no separate explanations.
- Checks: proton uud = 2+2−1 = +3 → Q = +1. Neutron udd = 2−1−1 = 0 → Q = 0. π⁺(ud̄) = 2+1 = +3 → +1. **Integer charge for all hadrons follows automatically.**

**Confinement test (measured)**: quark-type (m = +2) + sea starts at readable fraction f_read = 0.0000 and becomes readable only to the extent that the walk forms m ≡ 0 composites (hadron analogues) (→ 0.2504 over 4000 collisions). Electron-type (m = +3) + sea starts at f_read = 1.0000 (free from the outset). **The isolated quark-type system remains at f_read = 0.0000 even after 4000 collisions** — "a single quark cannot be observed" is realized as a consequence of readout (Figure P9).

![Figure P9](figs/fig_p9_confinement_v1.png)

**Figure P9**: The confinement = mod-3 readability test.

Correspondence with the standard theory: this is the same structure as the superselection by the SU(3) center Z₃ — "only triality-0 states are free" [10][11]. But whereas in the standard theory Z₃ enters as an axiom on the dynamics side (the center of the gauge group), in the present hypothesis mod 3 comes from the **readout side**, as the measured dominant denominator of the observation clock — one of the novelties of this paper (§9).

**Two precisions**: (i) two logics in this section must be distinguished — "the denominator-3 clock rectifies the charged walk to ±1" is a measurement; "identifying the physical charge as m/3" is a hypothesis. To promote the identification to a law requires a dynamical test of **response** to charge — whether the interaction strength Γ_int of species with different m scales as m/3 (or (m/3)²) (§13). (ii) The confinement of this paper is **readout-level confinement** (single-species unreadability); the separation-energy growth, flux tubes, and asymptotic freedom contained in standard confinement are not yet derived. Precisely, it is a **readout precursor structure of confinement**; the test of whether the cost of separating two unreadable species grows with distance (§13) is the promotion condition to dynamical confinement.

## 5. Pillar 3: the statistics rule — the Z₂ clock cover

### 5.1 Statistics does not live in spatial rotation (measured)

The global rotation of the channel doublet is an exact symmetry of the dynamics (all observables machine-zero flat in the control experiment; the rotation invariance of Im(b̄a) verified analytically), and the spatial rotation charge is ℓ = ±1 for all observed modes (after calibration, deviations within ±0.01 at N = 6, 8) — neither half-integers nor ℓ = 2 appear among occupied modes. The two-valuedness of statistics lives neither on the channel side nor on the spatial side.

### 5.2 Statistics lives in the double cover of the clock (measured)

Tracking the state phase of a species along the sequence of observational recurrences (peaks of the autocorrelation of the observable field s = Im(b̄a)), **the phase of the charged species is exactly quantized to {0, π} at every recurrence point** (Φ/π = ±0.999 or +0.002; 3-digit precision; no continuous values) — the state commutes between the +1 and −1 sheets (Figure P8). The neutral non-projected bundle shows a continuous phase; the presence or absence of quantization separates the species.

![Figure P8](figs/fig_p8_z2_quantization_v1.png)

**Figure P8**: Z₂ phase quantization: the charged species has Φ ∈ {0, π}; the neutral bundle is continuous.

The cross-table experiment (4 cells of χ parity × winding {1, 2}) showed that **the covering degree depends only on χ parity, not on winding (charge)**: fermionic species (even χ frequency = odd harmonics) have covering degree 2 (Z₂, both sheets); bosonic species have degree 1. Furthermore, a pure m = 0 fermionic-band species (winding-shift construction; η occupancy (0, 1.000)) also showed covering degree 2 with Qz2 = 1.00 — **the row of the neutral fermion (neutrino) is established** (Figure P10).

![Figure P10](figs/fig_p10_nu_spinstat_v1.png)

**Figure P10**: Establishment of the ν row and the spin-statistics correspondence.

### 5.3 Hypothesis (Pillar 3) and prior contrast

We propose **Fermi/Bose distinction = clock covering degree** (1 = bosonic / 2 = fermionic). This is the counterpart of the structure by which Finkelstein–Rubinstein [12] derived soliton spin-statistics from the double cover of configuration space — but here **the double cover lives not in space (rotation/exchange) but in the clock (state period = 2 × observation period)**. The same structure connects to the previously measured fidelity of order 248 / observation 124 [6][7] (recurrence invisible at F₁₂₄; F₂₄₈ = 1).

## 6. Pillar 4: the mass rule — mass is a relation to the sea

### 6.1 Equivalence theorems (three, measured)

- **Winding-shift equivalence**: e^{imη} is an exact symmetry of the dynamics (b̄a invariant); isolated pure address species are all unitarily equivalent — measured: mass², polarization, and s_z agree to machine precision across all addresses. **Mass and spin are not attributes of the isolated address.**
- **Divisor-class theorem**: in the sea, properties differentiate — but **only down to the divisor class gcd(m, 16)** (mass² of the odd addresses {1,3,5,7} agrees to 5 digits; {2,6} agree; {4} is distinct; Figure P7). The reason is exact: the unit automorphisms m → um of Z₁₆ fix the sea (m = 0), so sea coupling cannot distinguish unit addresses in principle. Consequence: **the distinction between ±1 and ±5, ±7 does not exist in the dynamics; only the mod-3 readout breaks the unit class** — an independent corroboration of Pillar 2 (readout rectification).
- **χ-band translation symmetry**: placing the same-winding species in χ bands (10,12,14)/(30,32,34)/(50,52,54) leaves mass² identical to 5 digits — band position does not differentiate species either (§10, weakness 1).

![Figure P7](figs/fig_p7_divisor_class_v1.png)

**Figure P7**: The divisor-class theorem: in-sea properties depend only on gcd(m, 16).

### 6.2 Hypothesis (Pillar 4)

**Mass is a relational quantity with the sea (vacuum).** Only the species locked to the clock (ρ = 1/1) is massless and stable (the photon-like ground species; mass degree 10⁻¹¹ measured); content away from the clock carries mass as resonances, and the sea differentiates the individuality of species (lifetime, polarization, mass) at the level of divisor classes. The view of the vacuum-as-sea as co-determinant of particle properties resonates with the lineage of emergent gauge fields and emergent fermions [13] (string-net condensation; superfluid ³He), but this paper differs in having measured the sea dependence of properties in the form of equivalence theorems.

## 7. Pillar 5: the periodic table of waves — double signboard and the 62-species assignment

### 7.1 The double signboard

- **Table A (isolated stable-element table)**: pure address species on the U^144 clock lattice, all stable by self-replication (measured). Algebraic classifiers: readable charge (mod 3), doubling orbit, neutralization steps, η parity. A second parity structure: odd addresses resist doubling-neutralization maximally (4 steps), even ones fall fast (±4 in 2 steps, +8 in 1).
- **Table B (in-sea lifetime table)**: sea coupling selects species — charged long-lived (τ ~ 10⁴), neutral transients short-lived (τ ~ 10²⁻³), composites neutralized. **Observed "particles" are this quasi-stable hierarchy.**

**The native periodic table (model-specific, measurement-based master table)**: we first write the table using only this paper's classifiers, without Standard-Model particle names. The Standard-Model assignment (§7.2) is a hypothesis on top of this table.

| Raw winding class | mod-3 readable | Neutralization steps | Clock cover (F-band/B-band) | In-sea mass² (divisor class) | Isolated stability | In-sea lifetime |
|---|---|---|---|---|---|---|
| m = 0 | 0 (free, neutral) | 0 | 2 / 1 | — (sea, ground species) | stable | ground species stable |
| odd (units) {±1,±3,±5,±7} | ±1 or 0 | 4 (max) | 2 / 1 | 0.787 (5-digit agreement) | exactly stable | τ ~ 10⁴ (charged) |
| {±2,±6} | ∓1 or 0 | 3 | 2 / 1 | 0.764 | exactly stable | retention 0.226 |
| {±4} | ±1 | 2 | 2 / 1 | 0.722 | exactly stable | retention 0.354 (longest) |
| {+8} | −1 | 1 (shortest) | — | unmeasured | — | unmeasured |

The covering degree is determined by χ parity, not by the raw winding class (cross-table measurement), hence an independent column.

### 7.2 The 62-species Standard-Model assignment

Confidence tags: **E** = measurement anchor (directly supported by measurements in this or published papers) / **S** = structural correspondence (mechanism matches; quantitatively unverified) / **H** = hypothesis (future verification target).

| Particle | States | Winding m (Q_em = m/3) | Cover | Conf. | Remarks |
|---|---|---|---|---|---|
| u, c, t | 18 | +2 (+2/3) | 2 | **E** | m ≢ 0 → confinement (Fig. P9 measured). Color = mod-3-unreadable raw-winding residue (ledger W) |
| d, s, b | 18 | −1 (−1/3) | 2 | **E** | same |
| e, μ, τ | 6 | −3 (−1) | 2 | **E** | free (m ≡ 0). Elementary-charge address sin²(23π/124) [6] |
| ν ×3 | 6 | 0 (0) | 2 | **E** | neutral fermion (established by measurement, §5.2) |
| γ | 1 | 0 | 1 | **E** | ground species ρ = 1/1 (massless, stable, measured) |
| g | 8 | color pairs (mod-3 neutral, raw winding nonzero) | 1 | S | color non-singlet → unreadable alone = confinement. The octet's generation not yet derived |
| W±, Z | 3 | ±3, 0 | 1 | S | t-winding (own-clock detuning) → massive (consistent with entrainment measurements; dynamics unverified) |
| H | 1 | 0 | 1 | S | collective mode of the sea (condensate). ℓ = 0 |
| G | 1 | 0 | 1 | H | only as a composite of two ℓ = 1 quanta (§8.4 prediction). ℓ = ±2 |
| (generation axis) | — | — | — | H | the axis differentiating 3 generations is unidentified (§10) |

Total **62 species** (36 quarks + 12 leptons + 8 gluons + γ + 2 W + Z + H + G). Antiparticles are m → −m (charge conjugation = η reversal; pair production measured [3]). The generation axis g = 1, 2, 3 is unresolved (§10). Per-row confidence is spelled out in the attached periodic table v0.4.

### 7.3 The relation between the two finite structures (open structural problem)

Two finite structures are measured in this paper: the universal clock U^144 (§3) and the winding register Z₁₆ (§4.2). The factor structure is 144 = 2⁴·3², 16 = 2⁴, 144/16 = 9 = 3², and this paper's classifiers correspond to the factors — **2 = statistics cover (§5), 3 = charge readability (§4), 16 = raw winding cycle (§4.2)**. Whether this coincidence is accidental, or a resonance structure in which both follow from the update rule (a resonance classification of clock order × internal winding order), is the largest open structural problem of this paper, registered in §13 as a task to be tested **as an operator recurrence structure**, not as numerology. Formalization: bring every readout operation into the same operator representation and measure the order ord(U|_{H_i}) per subspace — the true periodic law may be the family of "orders that a single operator has on different readout subspaces." The true periodic law of the "periodic table" probably lives here.

## 8. Explanatory power — connections to existing problems

### 8.1 Resolution of the previous paper's open problem, "charge does not stabilize at ±1"

The previous paper [5] §9.5 recorded as an open problem that "charge (winding) spreads over multiple quantized values and does not stabilize at ±1." Under the present hypothesis this is no contradiction: on the **state side**, the sea-driven walk (measured, §4.1) keeps dispersing windings. On the **readout side**, the denominator-3 observation clock folds all its products back to |q| = 1 (measured, §4.3). "Multivaluedness of charge" and "universality of the elementary charge" coexist without contradiction — an open problem turns into a prediction of the hypothesis.

### 8.2 The 1/3-charge problem of quarks

Under Q_em = m/3, fractional charge (why 1/3?) and confinement (why never alone?) are a single mechanism — mod-3 readability (§4.4; measured, Fig. P9).

### 8.3 Antimatter

Charge conjugation m → −m is a symmetry of the spectrum, and pair production is a necessity of the walk (published measurement [3]).

### 8.4 The weakness of gravity (a bold prediction — classified as a prediction, not a pillar)

The measurements are two: the rotation-generator spectrum of the condensate has no linear ℓ = 2 slot (calibrated ratio at most 1.000, exact); and the spectrum of the squared readout carries a coherent 2ω quadrupole line at all N (ratio 2 ± 6%). The strong hypothesis that follows: **if a gravitational spin-2 mode exists, it appears not as a first-order degree of freedom but only as a second-order composite of two ℓ = 1 quanta.** Beyond that — the quadratic suppression of the coupling (consistent with the prior derivation α_G = (m/M)² [15]), universal attraction, and the isomorphism with the double copy [16] — are registered as **predictions** for future dynamical verification. If the prediction is right, the reason gravity is weak and the reason it is hard to verify are one and the same structure.

## 9. Relation to prior work, and the delimitation of novelty

Each component of the hypothesis has a strong prior lineage: winding = charge (Skyrme [17]; Kaluza–Klein [18]); mod-3 confinement (triality [10][11]); the double cover of statistics (Finkelstein–Rubinstein [12]); vacuum = sea (Wen; Volovik [13]); gravity = square (KLT/BCJ [16]); full-assignment attempts (the Eightfold Way [9]; preons [14]). **This paper claims no novelty for the components.** The novelty is limited to three points:

1. That these components come out **simultaneously from measurements of a single anonymous dynamics** that assumes no particle concept (confluence).
2. The presentation of a counterpart in which charge unit, confinement, and elementary-charge uniqueness come not from a dynamics-side axiom (the center of a gauge group) but from the **divisor structure of the observation clock (the readout side)**.
3. The measurement that the double cover of statistics lives not in space (rotation of configuration space) but in the **clock** (state period = 2 × observation period).

## 10. Weaknesses and unverified parts (honest registration)

1. **The origin of generations is unresolved**: the naive hypothesis generation = χ-band position was refuted (χ-band translation symmetry, §6.1). Band ratios, crystallographic orders {3,4,6}, and composite hierarchy remain candidates.
2. **The W/Z rows**: direct production is impossible in the low-energy laboratory (N ≤ 16, small sea) — O(1) detuning decays at the register edge (structural). The 1/δ² scaling test of the effective contact vertex (the Fermi-theory road) is specified but not yet executed.
3. **The graviton row**: the quadrupole 2ω line exists at all N, but the pre-registered criterion (|ratio−2| < 0.02) passed only at 1 of 3 — frame calibration is incomplete.
4. **Dependence on the register order ne = 16**: the details of Table A (neutralization steps etc.) are functions of ne. Re-examination with variable ne is needed (also a testable prediction: "the period of the periodic table is a function of the register order").
5. **The unification of Q_em = m/3 with the transmission charge √(4πα) [6][19] is incomplete** — the relation between the "unit" (3 windings) and the "value" (0.302822) of the elementary charge is a task for the next paper.
6. Part of the sea constructions (artificial projected seas) was found to lie outside the conservation class (Nyquist fold-back, §4.2). The corresponding quantitative results are treated with reservation; qualitative conclusions were cross-checked in conserving-class series.

## 11. Falsification conditions

1. If an observation clock with a denominator other than 3 is constructed that rectifies the entire set of charged-walk products to a unique elementary-charge unit, Pillar 2 is rejected.
2. If, in the same condensed phase, with the same update rule and the same settling criterion, a series is obtained whose ω deviates systematically from π/72 in the long-window limit, Pillar 1 is rejected.
3. If a species is constructed whose covering degree varies independently of χ parity, Pillar 3 is rejected.
4. If a measurement shows species properties (mass, lifetime, polarization) differentiating without a sea, Pillar 4 is rejected.
5. If a construction is found in which an m ≢ 0 (mod 3) species becomes readable alone, the confinement rule is rejected.

## 12. Record of refutations and corrections (13 items, summarized)

All hypotheses and instruments refuted in the course of this work are recorded: (1) short-window sideband addresses are transients, not stable species (vanish in long windows). (2) The mass = detuning² law was an apparent relation on transients. (3) The naive form of ±1 fixed-point uniqueness (every pure m is a fixed point; uniqueness does not follow). (4) No mass-width correlation (r = −0.07). (5) The sea-construction convention (χ-analytic projection) hypothesis failed — the true cause of conservation violation is Nyquist fold-back (§4.2). (6) ±1 uniqueness by escape rates is undecidable in principle, by the divisor-class theorem. (7) The kinematic version of SU(2) recurrence discrimination is spinorial for all states (trivial). (8) Single-point covering judgments pick up accidental proximity to π (calibrated to recurrence-sequence judgment). (9) No linear ℓ = 2 slot for the graviton (reinterpreted as composite). (10) Generation = χ-band position refuted. (11) The band-shift method of W/Z detuning injection failed (effectively linear dispersion). (12) The projection method of ν construction failed (bundle hair m = +1 was the cause; solved by the shift method). (13) The earlier projected seas were amplified numerical residue (legitimate as m = 0 fields; qualitative results unaffected; recorded). Details are in the attached analysis notes.

## 13. Requirements handed to the unified dynamics (specification for the next paper)

To verify this classification hypothesis at the dynamics level, the two currently separate dynamics must be unified into a single **universal interaction function**, and the 62 species must be shown to be generated, to decay, and to interact as its spectrum. The measurements of this paper fix the requirements that the function must satisfy:

1. **Anonymity**: wave forms are the only input; no special-casing by species (no IF statements).
2. **All channels simultaneously**: electromagnetic (winding exchange), weak (clock detuning), strong (mod-3-unreadable raw-winding residue), and gravitational (second-moment coupling) must come out not as separate terms but as different readouts of the same operation, driving the three planes (xy = spin / zR = mass / tQ = charge) [20] simultaneously.
3. **Automatic conservation**: cyclic winding conservation (§4.2), closure conservation, and norm conservation must hold from the structure of the operation, not by design.
4. **Reproduction of the universal clock**: the U^144 clock lattice with ω = π/72 (§3) must appear as spectrum.
5. **Emergent selection rules**: the sum rule m* = 2m_B − m_s, the doubling-and-negation ladder, Z₂ phase quantization, and the divisor-class structure must come out without hand placement.
6. **Exact computability**: preserve anonymity ⇒ computation linear in state size [4] (closed form or integer-exact).
7. **Specialization to the two systems**: the two-body closed-form solution [3] and the N-body relational-wave dynamics [2] must be re-derived as restrictions/limits of the function.
8. **Falsifiable form**: if any of the 62 assigned species (§7.2) fails to appear in the spectrum, the periodic-table hypothesis is rejected row by row.

### 13.2 Verification program (top four)

In parallel with constructing the universal function, we register, in order of priority, experiments that can directly strengthen or refute this paper's hypothesis:

1. **Dynamical derivation of denominator 3, and charge response**: the dominance of denominator 3 currently rests on locking measurements [20]. If an analytic derivation is obtained of why the denominator-3 Arnold tongue is maximal in this dynamics, the "3" of Q_em = m/3 advances from an empirical selection to a dynamical necessity. In addition, measure whether the exchange quantities and phase shifts of two species (m₁, m₂) scale as the product m₁m₂/9 (a Coulomb-type charge-product test) — if it passes, m/3 is promoted from bookkeeping to the physical charge appearing in interactions.
2. **Clock covering degree and exchange sign in a single experiment (the decisive experiment to be performed before constructing the universal function)**: connect covering degree 2 ⟺ exchange sign −1 of two identical species in one experiment. If the exchange operation E and the one-clock-revolution operation T are shown to belong to the same nontrivial class of the same double cover (E ≃ T, E² = T² = I), then half-integer recurrence, exchange sign, exclusion, and clock covering degree contract to a single Z₂, and Pillar 3 advances from "resembling statistics" to "a reconstruction of the structure of statistics."
3. **Deformation law of the table under variable register order**: with ne = 12, 18, 24, 32, …, measure whether divisor classes, neutralization steps, readability, and covering structure follow gcd(m, ne)/Div(ne). If they follow, the divisor-class theorem generalizes; if not, 16 hides a structure of its own.
4. **Blind predictions**: without using Standard-Model particle names, predict from the model alone the unoccupied classification cells, lifetime rankings, charges, covering degrees, and decay destinations, and only afterwards compare with the known particle table — a test format with the Ω⁻-type decisive power of the Eightfold Way.

## 14. Reproducibility

All experiments are published as deterministic scripts (28) in the same folder as this paper. The correspondence table (実験一覧_v1.md) gives the complete mapping of scripts ↔ result JSONs (25) ↔ figures (10) ↔ claims. The three items originally executed interactively were formalized as scripts, and reproduction agreement with the saved results was machine-verified. Judgment criteria were fixed and recorded in each script header before execution.

## References

**Self-citations (sources of dynamics and measurements)**
[1] N. Kihara, Thought Experiment "On the R Axis," Zenodo (2026). doi:10.5281/zenodo.19902677
[2] N. Kihara, Causal Separation of the Metastable Phase by Two-Stage Seed Removal (Paper 8), Zenodo (2026). doi:10.5281/zenodo.21614402
[3] N. Kihara, Generation of Fermionic Structure — the Universal Inelastic Map, Zenodo (2026). doi:10.5281/zenodo.21808091
[4] N. Kihara, How Inflation Ends Is How Matter Generation Begins (N-Body Connection), Zenodo (2026). doi:10.5281/zenodo.21809814
[5] N. Kihara, Genesis of the Three Spatial Axes and Proper Time, Zenodo (2026). doi:10.5281/zenodo.21816651
[6] N. Kihara, Counting Readout in an Anonymous Two-Channel Closed Wave System (Trilogy B), Zenodo (2026). doi:10.5281/zenodo.21763997
[7] N. Kihara, The Generation Structure of Fermions (Paper 9), Zenodo (2026). doi:10.5281/zenodo.21766706
[15] N. Kihara, Re-Reading Centered on Reality — the Derivation-Replacement Edition, Zenodo (2026). doi:10.5281/zenodo.21765367
[19] N. Kihara, Two-Grammar Decomposition in an Anonymous Two-Channel Closed Wave System (Trilogy A), Zenodo (2026). doi:10.5281/zenodo.21763995
[20] N. Kihara, Locking Dynamics in an Anonymous Two-Channel Closed Wave System (Trilogy C), Zenodo (2026). doi:10.5281/zenodo.21763999

**External references**
[8] W. Thomson (Lord Kelvin), "On Vortex Atoms," Proc. Roy. Soc. Edinburgh 6, 94 (1867).
[9] M. Gell-Mann, "The Eightfold Way," Caltech Report CTSL-20 (1961); M. Gell-Mann and Y. Ne'eman, The Eightfold Way (Benjamin, 1964).
[10] K. G. Wilson, "Confinement of quarks," Phys. Rev. D 10, 2445 (1974).
[11] G. 't Hooft, "On the phase transition towards permanent quark confinement," Nucl. Phys. B 138, 1 (1978).
[12] D. Finkelstein and J. Rubinstein, "Connection between spin, statistics, and kinks," J. Math. Phys. 9, 1762 (1968).
[13] M. A. Levin and X.-G. Wen, "String-net condensation," Phys. Rev. B 71, 045110 (2005); G. E. Volovik, The Universe in a Helium Droplet (Oxford, 2003).
[14] H. Harari, "A schematic model of quarks and leptons," Phys. Lett. B 86, 83 (1979); M. A. Shupe, Phys. Lett. B 86, 87 (1979).
[16] Z. Bern, J. J. M. Carrasco, and H. Johansson, "Perturbative quantum gravity as a double copy of gauge theory," Phys. Rev. Lett. 105, 061602 (2010); H. Kawai, D. C. Lewellen, and S. H. H. Tye, Nucl. Phys. B 269, 1 (1986).
[17] T. H. R. Skyrme, "A unified field theory of mesons and baryons," Nucl. Phys. 31, 556 (1962).
[18] T. Kaluza, Sitzungsber. Preuss. Akad. Wiss. (1921) 966; O. Klein, Z. Phys. 37, 895 (1926).
Note: S. Coleman, "Quantum sine-Gordon equation as the massive Thirring model," Phys. Rev. D 11, 2088 (1975), representative of the bosonization lineage, is referenced via [7].
