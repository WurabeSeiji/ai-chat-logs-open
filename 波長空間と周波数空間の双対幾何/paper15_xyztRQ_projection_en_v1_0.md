# Paper 15: Projection onto xyztRQ — Construction and Machine Verification of the Map and Inverse Map from the Discrete Dual System to Spacetime Coordinates

**Author**: Noriaki Kihara
**Date**: June 12, 2026
**Version**: v1.0 (final; v0.2 incorporated the 14 referee points of the claude.ai review, v1.0 the 3 residual points of the second-round review (accept))
**Series**: Dual Geometry of Wavelength Space and Frequency Space (sequel to Papers 1–14)
**Note**: originally reviewed and accepted as "Paper 13"; renumbered to **Paper 15** upon insertion of the new Papers 13 (position and time resolved) and 14 (time direction and the amplitude stage) — the only content change is the redirection of citations. The numbering is independent of other series.

---

## Abstract

From the state data of the discrete dual system constructed in Papers 1–12 (three axioms — reciprocal duality νλ=1, zero point ½, asymmetric single bit — plus the operating principles of configuration reading and the record principle), we construct an **explicit projection and inverse projection** onto the spacetime coordinates $(x,y,z,t)$ and the ledger coordinates $(R,Q)$, and machine-verify every stage. Main results: (1) **Spatial coordinates** are given by a mixed-radix $[4;k,k,\dots]$ quarter-digit expansion ($k$ = the odd genealogy ratio); the digits are read from the phase classes (Z₄) of the transcription lines of the record by peeling decode. The only genealogy for which the readout closes within the digital alphabet of the record is the **ternary one ($k{=}3$)** (information count $2\,\mathrm{bit}\ge\log_2 k$). (2) **Time is a derived quantity**: $\nu_t=\sqrt{s}$ (the norm clock), $t=\sum 1/\nu_t$ (event accumulation). It has no independent inverse map — not a defect, but the map-level realization of the theorem (Paper 13, Theorem 6) that $t$ is not a fundamental degree of freedom. (3) **R and Q require no map**: $R=\sqrt{s}$ and $Q=\varepsilon$ are B₄ invariants shared by all markings, and the image of the projection lies on the constraint surface $\Sigma=\{R^2=s,\ Q=\varepsilon\}$ (codimension 2) — the formal statement of "subjective 4D within six background axes." (4) **Inverse map**: space is constructively invertible (rounding loss below the censorship quantum; composite round trip 500/500), time admits only event counting, R/Q are constraints. (5) **Compatibility with conservation laws**: the value of $s$ read from the record alone is additive across decays (60/60 — a check that could have failed), and the multiplicative conservation of the charge-like $Q$ is a **theorem of parity arithmetic on the support** (Theorem 8), enforced mechanically by transport: 39.2% of raw cell combinations violate ε, yet zero of the amplitude-carrying paths do. (6) **4D lift**: product waves possess no pure-axis lines (sum/difference lines only), and the simultaneous half-wave shift of an even number of active axes is an identity of the state — the readout decodes exactly to that quotient (4-axis exhaustive: 32 classes, one-to-one), while parity-mixed composites break the degeneracy completely (all 256 configurations uniquely decoded). (7) **Kinematics**: uniform motion is a mixed-radix odometer of the position digits; the amplitude spectrum, wavelength content, R, and Q do not evolve in time (deviation $10^{-15}$) — only the phase evolves, linearly, and its rotation rate is the velocity (exactly $f\cdot v$). **Velocity exists nowhere in a single record; it exists only in the sequence of records.** This paper does not claim a re-derivation of standard theory. No physical identification is made.

---

## 1. Introduction

### 1.1 The question

Papers 1–12 constructed a discrete system from the reciprocal duality $\nu\lambda=1$, the zero point ½, and a single asymmetric bit, and showed that position (Paper 13, Theorem 1 and Definition 2 / Theorems 2a–2b), time (Paper 13, Theorems 5–6), and the stage for amplitudes (Paper 14, Theorems 1–3) are handled by the existing inventory alone. What remained was the **map itself**:

> Can the state data of the system be mapped to, and inverted from, the standard physical coordinates $(x,y,z,t)$ and the ledger coordinates $(R,Q)$ as an **actual operation**?

There is a distance between "it should be possible" and "the formalization is exhibited and machine-checked." This paper closes that distance: the projection is given as a five-stage operation, each stage machine-verified (methods listed in Appendix A), and beyond "values round-trip" we demonstrate on the output side that the map **commutes with the conservation laws** — that R is a derived quantity governed by an energy-like conserved ledger and that t is derived from R.

### 1.2 Position within the series

This paper is the **coordinate-map installment** of the correspondence-dictionary arc (Paper 12, §8); the measure line is outside its scope and remains in active development — the arc is not yet complete. The transport principle — "theorems of standard theory transport without re-derivation when their premises have verified counterparts in this system" — becomes operational only when the coordinate map connecting the two sides has been constructed. The materials are the state structure of Papers 5–7, the internal time of Paper 8, the odd-nesting theorem of Paper 9, and the readout theory of Paper 12; no new axiom is added.

### 1.3 What is not claimed

We do not claim (i) a re-derivation of standard theory, (ii) a construction of common time for multi-fragment systems (synchronization), or (iii) a derivation of dynamics (why uniform motion occurs). (ii) belongs to the gauge-fixing mechanism of the record interface (a separate line); (iii) to the effective interpolation of the wave sector (Paper 12, §4.2). No physical identification is made.

## 2. Preliminaries (minimal recap)

State data: occupied cells $k\in\mathbb{Z}^4$ (configuration reading: the state *is* the set of occupied cells = the wave itself), branch (cos/sin) and sign per level, the genealogy chain (odd-ratio nesting $k_\ell$, Paper 9), and the ledger $s=\sum_i(|k_i|+\tfrac12)^2$, $\varepsilon=(-1)^{\Sigma|k_i|}$. The record: the hologram $I=\Psi^2$ — its physical content is the **line coefficients (relational data)**; the position-space picture is their Pontryagin-dual chart (space = the character group of the conserved-charge lattice; Paper 13, Theorems 1–2). Gauge data of the chart: the marking $u$ (choice of the time-read axis; 16 choices forming a single B₄ orbit (Paper 14, Theorem 2)) and the **per-axis orientations** $Z_2^{\,4}$ (choice of positive branch direction — the position-reading counterpart of the orientation fixing of Paper 14, §3). These are held fixed below.

**Emphasis**: a particle does not "remember" its state. Position = the phase offsets of its wave; R = the norm of its wavelength content; Q = the parity of its content — all are aspects that *constitute* the state. The only memory in the system is the record; every readout in this paper takes the record alone as input.

## 3. Definition of the projection

Fix a marking $u$. The three transverse axes provide $(x,y,z)$.

**D1 (spatial coordinates)**: for each axis $j\neq u$, the mixed-radix quarter-digit expansion

$$
x_j=\frac{1}{4}\Bigl[d_j^{(0)}+\sum_{\ell\ge1}c_j^{(\ell)}\prod_{m\le\ell}k_m^{-1}\Bigr]\bmod 1,
\qquad d^{(0)}\in Z_4,\ c^{(\ell)}\in Z_{k_\ell},\ k_\ell\ \text{odd}
$$

where $d^{(0)}$ is the level-0 quarter digit (cos=0, sin=1, −cos=2, −sin=3) and $c^{(\ell)}$ the position of the level-$\ell$ child comb within the parent quarter gap.

**D2 (time coordinate)**: $t=\sum_{\text{record events}}1/\nu_t$ with $\nu_t=\sqrt{s}$ — identical in construction to the internal time of Paper 8, §3.2 (numerically cross-checked).

**D3 (R, Q)**: $R=\sqrt{s}=\nu_t$, $Q=\varepsilon$. No map required (identity reads); complete B₄ invariants shared by all markings.

**D4 (inverse map)**: below we specialize the genealogy to the canonical $k_\ell=3$ (Theorem 4). For $x_j$: greedy rounding to depth $L$,

$$
d^{(0)}=\lfloor4x\rfloor\bmod4,\qquad
c^{(\ell)}=\Bigl\lfloor 3^{\ell}\bigl(4x-d^{(0)}-\sum_{1\le m<\ell}c^{(m)}3^{-m}\bigr)\Bigr\rfloor\bmod 3
$$

with the error guarantee $|x-\hat x_L|\le\frac14(\prod k_m)^{-1}\frac{k_L}{k_L-1}$ (below the censorship quantum = no loss of physical fact). In implementation, floating-point boundary points may yield a digit equal to $k$; we clip with $\min(c,k{-}1)$ (a measure-zero boundary convention that does not affect the bound). For $t$: no independent inverse (only the event count $n=\lfloor t\sqrt s\rfloor$). For $R,Q$: constraints (no freedom).

## 4. Basic theorems

**Theorem 1 (quarter-digit consistency)**: the translation $T_{1/4}$ shifts the digit $d\to d+1\pmod4$ (exhaustive, 4/4).

**Theorem 2 (completeness and rounding bound)**: the chain of position lattices $P_\ell=(4k_1\cdots k_\ell)^{-1}\mathbb{Z}/\mathbb{Z}$ is nested for any odd genealogy ($[P_\ell:P_{\ell-1}]=k_\ell$); D1 is a complete mixed-radix expansion on $T$ (20,000 random points within the bound). An important distinction: **the hierarchical scale of position is the wavelength genealogy ratio $k$, not the container ratio $1/(2R')$** (the latter makes the digit chain lacunary — rejected by the machine). Oddness originates not from density but from wave consistency (the odd-nesting theorem).

**Theorem 3 (constraint surface)**: the image of the projection is the codimension-2 surface

$$
\Sigma=\bigl\{(x,y,z,t,R,Q): R^2=s(k),\ Q=\varepsilon(k)\bigr\}.
$$

R and Q carry no inverse-map freedom — the formalization of the "invisible directions" (Paper 12, Suppl. 39/40). Virtual displacements (±1 borrowing) are off-surface moves permitted only between records (§7).

## 5. The readout operation (five stages of the forward map)

For axis $j$ and frequency $f=3^\ell$ (henceforth $k=3$), the transcription line:

1. **Line extraction**: $Z^{(\ell)}=2\langle I,\cos(2\pi 3^\ell x_j)\rangle+2i\,\langle I,\sin(2\pi 3^\ell x_j)\rangle$
2. **Peeling**: $\tilde Z^{(\ell)}=Z^{(\ell)}-\sum_{m<\ell}(\text{known transcription contributions of decoded shallower levels})$. **Even-drop theorem**: comb×comb interference falls entirely on even frequencies (odd×odd = even; leakage $2\times10^{-15}$), so odd readout lines consist of transcription terms only and the peeled amounts are exactly known.
3. **Quarter class**: $r^{(\ell)}=\frac{2}{\pi}\arg\tilde Z^{(\ell)}\bmod 4\in Z_4$ (all phases land exactly on the quarter lattice; quantization error $3.4\times10^{-15}$)
4. **Digit decode (recursion)**: $c^{(0)}=r^{(0)}$, $c^{(\ell)}=(r^{(\ell)}-3r^{(\ell-1)})\bmod4$
5. **Assembly**: substitute into D1.

Verification: two levels 12/12; three levels 36/36; full combs (all odd harmonics, cutoff 81) three levels 36/36 — all uniquely decoded.

**Theorem 4 (digital preference for the ternary genealogy)**: a quarter channel carries $Z_4=2$ bits per level. Unique decoding of a digit $c\in Z_k$ requires $k\le4$; among odd $k\ge3$ only $k=3$ qualifies ($k{=}1$ is a degenerate genealogy and excluded; $k{=}5$ exhibits the degeneracy $c{=}4\equiv c{=}0$, machine-demonstrated 4/20). **The only genealogy whose positions can be read out entirely within the record's digital alphabet is the ternary cascade** — coinciding, by an independent principle, with the genealogy selected by minimality in Paper 8 (the Z₂ prohibition of binary).

## 6. The 4D lift

A product wave $\Phi=\prod_j\varphi(x_j-a_j)$ has no pure-axis lines (its spectrum consists of sum/difference lines only); digits are obtained by solving the quarter classes of sum/difference lines mod 4.

**Theorem 5 (diagonal quotient)**: the simultaneous half-wave shift of an even number of active axes ($d_j\to d_j+2$) is an **identity of the wave** (2-axis 16/16; 4-axis: state classes = quotient by the even-sign-flip group $Z_2^3$, 256/8 = 32). The readout decodes exactly to this quotient (32 signatures, one-to-one) — a redundancy of coordinates, not an ambiguity of states (configuration reading: same wave = same state).

**Theorem 6 (full decoding via parity mixing)**: composites containing cells with an odd number of active axes break the degeneracy. The condition is a **covering of every nontrivial even flip by an odd-parity cell intersecting it oddly**; constructions satisfying it decode uniquely in full enumeration — 2-axis 16/16 and **4-axis 256/256**.

## 7. Compatibility with conservation laws — R and t are derived quantities

Of the six coordinates xyztRQ, only $xyz$ plus the event count are free. The checks in this section fall into three distinct strata (made explicit at referee point #2):

**(i) Bookkeeping consistency (true by construction)**: the decay path table enforces $s_{\rm kept}+s_c+s_d=s_{\rm parent}$ at records and the ±1 borrowing cancellation at virtual stages by definition. Its holding over all 118,944/118,944 contributing paths at s=9/11/13 (eight representative parents per s, covering all three orbit types) is a **consistency check that the path enumeration and the amplitude machinery implement the same bookkeeping** — it has no refutation power.

**(ii) A theorem (which could have failed, and is proven)**:

> **Theorem 8 (charge-like parity selection rule; the general-lemma form is canonical as Paper 14, Appendix A-3 — this theorem is its selection-rule expression)**: if a transport coefficient is nonzero, $\varepsilon$ is multiplicatively conserved. Proof: $C_c\neq0$ requires the target $v_p=n_a+n_b$ (with $n$ on the coefficient support, $|n_{a,j}|=|v_{a,j}|$), and integer-sum parity gives $|v_{p,j}|\equiv|v_{a,j}|+|v_{b,j}|\pmod 2$ per axis; summing over axes, $\varepsilon(v_p)=\varepsilon(v_a)\varepsilon(v_b)$. Applying this to both stages, $\varepsilon_{\rm parent}=\prod\varepsilon_{\rm daughters}$. ∎

Machine confirmation: zero violations over all 118,944 contributing paths. Control: among raw cell combinations (s=13, 11.4 million) **39.2% violate ε** — violating candidates abound on the support, and the readout dictionary annihilates their amplitudes exactly. Charge-like conservation appears neither as a postulate nor as an empirical selection effect, but as a **theorem of parity arithmetic on the support**.

**(iii) A check that could have failed (via the record)**: the occupied cell set is recovered from the hologram alone (exact configuration reading, 60/60) → additivity of $s_{\rm read}$ (60/60) → the chain "configuration → $R=\sqrt{s_{\rm read}}$ → $\nu_t=R$ → $t=\Sigma1/\nu_t$" holds at every stage **on the output side**. This is independent of construction (the decoding could have failed; the additivity could have broken).

> The projection is not a picture of coordinates but a **faithful display of the conservation ledger**. The conservation laws and the projection are not merely compatible — they are two faces of the same dictionary.

## 8. The inverse map and the composite round trip

Continuous $x$ → (D4 rounding) → digits → **actual construction of the state** (full combs, three levels) → (§5 forward map) → digits, $\hat x$: over 500 random points, **complete digit agreement 500/500**, maximum error 0.0277 ≤ bound 0.0417.

$t$: no independent inverse exists — $\nu_t$ is derived from the spatial configuration (the realization of the theorem of Paper 13 (Theorem 6); D2 is its implementation). This is consistent with the principle that "t merely appears to have been selected by observation" (the all-axes-symmetric principle). $R,Q$: six-tuples violating the constraints lie outside the image (= the invisible directions).

## 9. Condition matrix and robustness

| Scenario | Configurations | Result |
|---|---|---|
| Three heterogeneous fragments (fundamentals 1/3/5; peeling cascade) | 64 | **64/64** |
| Three fragments at continuous positions (no digit alignment; direct phase read) | 200 | **200/200** (error $2.5\times10^{-16}$) |
| Harmonic amplitude jitter ±20% + additive noise σ=0.5 | 192 | **192/192** |
| Three same-species fragments (worst case: shared lines) | **all 220** | without half-wave pair **160/160** unique / with **60/60** state-identical (Theorem 7) |
| Three 4D composites (configuration reading) | 60 | **60/60** |

**Theorem 7 (half-wave pair = state identity)**: when a same-species pair sits at half-wavelength offset, the sum of the square combs vanishes identically (maximum $1.5\times10^{-14}$) — i.e., **a configuration containing such a pair is the same wave, hence the same state, as the configuration without it** (like Theorem 5, an element of the kernel of the configuration→wave map). The readout returns the state correctly; the "degeneracy" is mere label redundancy. This separates the present theorem from the open problem of §12-2 ($k\ge5$ deep digits = **different states** that cannot be read = censorship): the former is state identity, the latter a limit of readability.

## 10. Kinematics: the map of uniform motion

This section verifies uniform motion by **two complementary constructions** (they are not a single object — referee point #10):

**O1 (digit odometer: readout of hierarchical states)**: each tick ($\Delta t=1/\sqrt s$) advances the deepest quarter quantum by one, and the **hierarchical state, re-assembled digit-aligned at every tick**, is read. The hierarchical decode returns digit = counter value (base [4;3,3]) at all 36 ticks; $\hat x_n=n/36$ exactly. The world line is an exact straight line through lattice points, $v=\Delta a\cdot\sqrt s$.

**O2 (rigid drift: spectrum of a single comb)**: a single comb is rigidly translated continuously and the time evolution of its line spectrum is measured:

- **The amplitude spectrum, wavelength content, R, and Q do not evolve** (maximum deviation $1.9\times10^{-15}$ over 36 ticks)
- Only the phase evolves: exactly linear rotation at rate $f\cdot v$ (residual $\sim10^{-15}$; ratio $f{=}3/f{=}1=3.000000$)

The two are distinct verifications: a rigidly drifting hierarchical composite leaves the quarter lattice at intermediate ticks, and in digit-aligned staircase motion the phases of shallow lines advance stepwise (exact linearity holds only for the deepest line). The "rigid drift" of Figure 4 refers to O2.

> **Velocity appears nowhere in a single snapshot** — both |spectrum| and support are identical to the static state; velocity exists only in the sequence of records. On the ledger: |spectrum| = the conserved quantity of free motion (momentum-like), phase = position, phase rotation rate = velocity. Readable velocities lie on a rational lattice; continuous velocity is its completion (the same architecture as for position).

## 11. Figures

- **Figure 1** (`paper15_fig1_readout.png`): the readout in both spaces — transcription lines in frequency space (odd lines only), superposed combs in position space, fragments on the constraint surface in xyztRQ space, condition-matrix results
- **Figure 2** (`paper15_fig2_xyztRQ.png`): xyztRQ space — the $(x,y,R)$ 3D arrangement (positions machine-decoded), world lines (tick spacing $1/R$, branching at decay), the $(R,Q)$ ledger plane, the constraint surface $R^2=s$ with virtual ±1
- **Figure 3** (`paper15_fig3_conservation.png`): conservation — spectra before/after the decay $9\to(5,1,3)$, the $R^2$ ledger (virtual +1 borrowing), the ε control (39.2% vs 0), verification summary
- **Figure 4** (`paper15_fig4_uniform_motion.png`): uniform motion — the x–t world line, the spacetime record (rigid drift, O2), spectral invariance, linear phase rotation

All figures are exact computations (no schematic drawings).

## 12. Limits and outlook

1. **Common time for multiple fragments**: D2 is the proper clock of a single fragment. Synchronization = the (1,3) aggregation (gauge fixing by the record interface) is not yet constructed — the system's most important remaining item.
2. **Analog digits for $k\ge5$**: deep digits of non-ternary genealogies require sub-quarter amplitude-ratio reads, outside the digital alphabet. Whether "an unreadable digit exists as position" touches the principle of configuration reading (an open problem; the canonical genealogy is ternary, so the main line is unaffected).
3. **Dynamics**: the motion of this paper is kinematics (the map and readout of a schedule); deriving "no force → uniform motion" belongs to the wave-sector interpolation. Acceleration, relative motion, and collisions are not constructed.
4. General-$s$ algebraic proofs of single-sign support, torsor residues, etc. share the series' proof queue (s=9/11/13 are exhaustive facts).
5. No physical identification is made. Standard theory is exactly correct; this paper constructs the map portion of "a different mapping of the same results."

## Appendix A: Verification summary (all reproducible)

| Claim | Verification | Method | Script |
|---|---|---|---|
| Quarter-digit shift consistency | 4/4 | exhaustive | supplement62_projection_formalization.py |
| Mixed-radix completeness and bound | 20,000 pts, k=3,5 | random | same |
| Constraint-surface invariance | B₄ samples | sampled | same |
| t = internal time of Paper 8 | numerical match | cross-check | same |
| Readout (fundamental combs) | 12/12, 36/36 | exhaustive | supplement63_digit_readout_protocol.py |
| k=5 degeneracy (preference theorem) | 4/20 | exhaustive | same |
| Full combs, peeling | 36/36, err 3.4e-15 | exhaustive | supplement63_full_comb_collision_test.py |
| Even drop | leakage 2e-15 | exhaustive (cutoff 81) | same |
| 2-axis diagonal quotient, parity mixing | 16/16, 8 classes, 144, 16/16 | exhaustive | supplement64_4d_product_wave_test.py |
| **4-axis full check** | 32 classes one-to-one, **256/256** | exhaustive | paper15_appendix_4axis_check.py |
| Composite round trip (inverse) | 500/500 | random | supplement62_inverse_roundtrip_test.py |
| Bookkeeping consistency (R² additivity, ±1) | 118,944/118,944 | exhaustive (true by construction, §7-i) | supplement65_conservation_tests.py |
| Theorem 8 (ε selection rule), machine confirmation | 118,944/118,944; control 39.2%→0 | exhaustive + theorem | supplement65_controls.py |
| Configuration reading, s_read additivity | 60/60 | random (final states) | supplement65_c3_exact_decoder.py |
| Condition matrix B/C/D | 64/64, 200/200, 192/192 | exh./random/random | supplement66_condition_matrix.py |
| Same-species classification | 160/160 + 60/60 | **exhaustive (220)** | same + this paper |
| Uniform motion O1/O2 | 36/36, 1.9e-15, exact f·v | exhaustive (36 ticks) | supplement67_uniform_motion.py |

---

**Acknowledgments / history**: the verification of this paper was carried out under the two-party independent verification protocol of Claude Code (local machine verification) and claude.ai (independent re-computation and review). Source material: Supplements 62–67 (June 12, 2026); foundational theorems now in Papers 13–14. v0.2: 14 referee points incorporated (review with independent re-implementation reproducing all computational claims); v1.0: 3 residual points of the second-round review (accept).

No physical identification is made.
