# Survey: From One Relation, One Constant, and One Bit — A Dictionary to Standard Theory and the Balance Sheet of an Axiom Exchange

**Author**: Noriaki Kihara
**Date**: June 12, 2026
**Version**: v1.0 (final, publication edition; strict counterpart of the Japanese v1.0. Source DOIs fixed upon the Zenodo publication of Papers 13–16; all entries of the dictionary in §5 are "published"; machine cross-checking of the citation map completed [Claude Code]. Repository snapshot Concept DOI acquired and inserted [10.5281/zenodo.20666114]. Drafted by claude.ai; DOI procedures and tex/PDF by Claude Code)
**DOI (this version)**: 10.5281/zenodo.20666133 / **Concept DOI**: 10.5281/zenodo.20666132
**License**: CC BY 4.0
**Series**: Dual Geometry of Wavelength Space and Frequency Space (survey; numbering undecided: candidate Paper 17 / standalone survey)

---

## Abstract

The starting point of this series is one relation (νλ=1), one constant (the zero point ½), and one literal bit (the asymmetric bit) — that is the entire information content of the axioms, supplemented by two operating principles that fix how states are read (configuration reading and the record principle). **What was not brought in**: the gauge principle, the complex number field and its conjugation, spacetime dimension, a time parameter, the symmetrization postulate, the measurement postulate, the Born rule, the action. Up to now the series has deliberately **refrained from speaking of any physical connection**, in order to avoid misreading and dismissal, building instead a system verifiable step by step. Yet the goal was constant from the beginning: to test whether a system with expressive power close to standard quantum theory and standard gauge theory can be derived from these extremely limited assumptions alone. This survey is the balance-sheet report on that goal. The method is not re-derivation but an **item-by-item verified correspondence dictionary** together with a **transport principle**. Main results: (1) the kinematic dictionary is filled, every entry machine-verified and **every source now published**. (2) **Five items that standard theory installs as axioms — symmetrization, time, measurement, dimension, charge conservation — this system possesses as theorems.** (3) Two entries remain blank — the measure and the action. Beyond this, the survey **dissects the points of incomplete connection themselves (§9) and identifies their structural causes**. A completely filled dictionary would amount to a re-notation and would carry no empirical content — **the points where the connection is incomplete are the only places where the two theories make different claims and observation can adjudicate**; the remainder is at once a deficiency and the seat of this system's decidability. **Conclusion**: what was initially aimed at was a mere toy model; the reader is asked to inspect the dictionary (§5) — we believe the system has connected, at a rather close level, to the basic requirements of standard quantum theory and of gauge theory. **Caution**: this series is in no way an attempt at a new theoretical physics or a unified theory. **It may merely be that structures close to standard quantum theory can be derived from a different set of simple assumptions** — that is the reading the author prescribes. No physical identification is made. As a postscript (§13) we include a retrospective of the starting point and a reply to the anticipated criticism that "this is just complex numbers in disguise." References and access routes to all sources are given in the bibliography — the claims of this survey are auditable not only in content but **in process** (end of §2).

---

## 1. How to read this document (first)

This survey and the series as a whole are **not**: a proposal of new theoretical physics, an attempt at a unified theory, or a rival hypothesis to standard theory. Standard theory is rigorously correct, and this series contests none of its standing.

What the series set out to test is a single question — **from how few, and how alien, assumptions can the structure of standard quantum theory re-emerge?** This is a question of the same shape as the reconstruction programs of quantum theory; what distinguishes this series is the extreme smallness of the assumptions' information content (§3: one relation, one constant, one bit) and the fact that every step proceeds in a discrete system that is machine-checkable. Hence if the claims here are right, standard theory loses nothing; if they are wrong, only this document loses. In that sense it is a safe document to read. "A model in which structures close to standard quantum theory might be derivable from a different set of simple assumptions" — this is the reading the author prescribes, and the series-wide disclaimer "no physical identification is made" is the technical expression of that reading.

## 2. Background — why the connection was not spoken of until now

Papers 1–16 and the 84 supplements deliberately avoided claiming any physical connection, for two reasons. **Avoiding misreading**: speak of the connection first and it will be conflated with an unverified identification. **Avoiding dismissal**: unverifiable claims are turned away at the door. The strategy was therefore split into two phases — Phase One: without speaking of connections, accumulate only a stepwise, machine-checkable system under two-party independent verification (complete). Phase Two is this survey: on the evidence of the accumulated dictionary, to account for the state of the connection squarely for the first time. The original goal (deriving rich expressive power from limited assumptions) was unchanged throughout Phase One — it was simply not spoken.

**Where the record lives (auditability of the process)**: the claims of this survey can be audited not only in content but in process. The entirety of Phase One — the 84 supplements, all verification scripts, the back-and-forth of two-party independent verification, and **the history of mutual correction including refuted predictions, retracted explanations, and reconciled definitions**, together with the review records of Papers 13–16 (two rounds each, eight reviews in total) — is preserved, with commit history, in the public repository (see References). The tables of this survey are thus not a sudden idea but the balance of an **accumulation of dated thought experiments and machine verifications**, and that accumulation is itself published as an object of third-party scrutiny. To keep the history, failures and corrections included, undeleted — that is the verification culture of this series, and the ground of this survey's credibility.

## 3. The accounting of assumptions — the lightness of the start is the first result

### 3.1 Axioms (three; about 2.5 when counted by information)

| # | Axiom | Information content |
|---|---|---|
| 1 | Reciprocal duality νλ=1 | one relation |
| 2 | Zero point ½ | one constant |
| 3 | The asymmetric bit | **literally one bit** |

The third axiom is a binary choice — the smallest input an axiom can be. The tally "effectively 2.5 assumptions" is the precise expression of this information count.

### 3.2 Operating principles (two)

**Configuration reading** (the state is the set of occupied cells = the wave itself) and the **record principle** (records are appended from absence to presence; an existing record is presently readable). These are rules of reading, not axioms of content; but in a critical accounting they should be counted, and we count them. What matters is that they double as the counterparts of the **state postulate and the measurement postulate** on the standard side — not an added cost but part of the exchange (§6).

### 3.3 The list of what was not brought in

The gauge principle. The complex field ℂ and its conjugation. The dimension of spacetime. An external time parameter. The symmetrization postulate. The measurement postulate (projection, state update). The Born rule. The action and the Lagrangian. — None of these is among the axioms. Of these, the gauge skeleton, the complex structure, the Z₄ phase, dimension, time, symmetrization, and measurement **emerge within the system** (§5, §7). The Born rule and the action are the remainder (§8, §9).

## 4. Method: dictionary and transport principle

The **dictionary** is a table of pairs {object of the system ↔ object of standard theory}; an entry counts as "filled" when (i) the system-side construction is established as a theorem or by machine verification, and (ii) the structural correspondence (algebra, order, transformation behavior) is confirmed. The **transport principle**: when all premises of a theorem T of standard theory correspond to filled entries of the dictionary, T is transported without re-derivation. The reason for not choosing re-derivation is fundamental: the proofs of standard theory are correct. What is needed is not a redoing of proofs but a verification of the counterparts of their premises. With the completion of the coordinate map (Paper 15), transport turned from a "principle" into an "operation." The dictionary asserts no physical identification — it asserts structural correspondence, and a different map producing the same results supports transport without identification.

## 5. The dictionary as it stands (the main table — the conclusion rests on this table)

Legend — status: **pub** = published on Zenodo (**every source is published as of this survey** — DOIs in the References). Verification: thm (proof given) / mach (machine verification).

| Object of the system | Standard-theory counterpart | Status | Source |
|---|---|---|---|
| Set of occupied cells (configuration reading) | state | pub · thm+mach | Paper 7 |
| Wave dictionary, R quantization | state↔wave correspondence | pub · mach | Paper 5 |
| Conservation of Σν² | energy-like conserved quantity | pub · thm+mach | Paper 6 |
| Exclusion (no double occupancy) | counterpart of the exclusion principle | pub · thm | Paper 7 |
| Internal time t=Σ1/√s | time parameter | pub · **thm (derived quantity)** | Papers 8, 13 |
| Position = sign sector (zero additions to states) | position observable | pub · thm+mach | Paper 13 |
| Space = Pontryagin dual of the conserved-quantity lattice | position representation | pub · thm | Paper 13 |
| Complex structure C⊂H (from the choice of time direction) | the i of complex Hilbert space | pub · thm | Paper 14 |
| Derived Z₄ transport phase | phase / skeleton of U(1) | pub · thm+mach | Paper 14 |
| Canonical-convention theorem (identical-species symmetrization) | **symmetrization postulate** | pub · thm | Paper 14 |
| Necessity of dimension 4 ({1,2,4,8} × squares) | **spacetime dimension (an input)** | pub · thm | Paper 11 |
| Measurement = appending + record theorem | **measurement postulate** | pub · thm+mach | Paper 12 |
| Charge-like parity selection rule | **charge conservation** | pub · thm+mach | Paper 15, Thm 8 |
| Coordinate map xyztRQ (forward, inverse, conservation-commuting) | coordinate system | pub · mach (all points independently reproduced) | Paper 15 |
| Gauge skeleton (emergent from counting) | gauge structure (**action not constructed**) | pub · mach | Paper 10 |
| Kinematics of uniform motion | kinematics of free motion | pub · mach | Paper 15 |

## 6. The axiom exchange balance — five items that standard theory installs as axioms, this system holds as theorems

| Status on the standard side | Status in this system | Source |
|---|---|---|
| Symmetrization postulate (identical particles) | **Theorem** (canonical-convention theorem) | Paper 14 |
| Time = external parameter | **Theorem** (t is derived: no independent ledger, no independent inverse map) | Papers 13, 15 |
| Measurement postulate (projection, update) | **Theorem** (appending + record theorem) | Paper 12 |
| Spacetime dimension = input | **Theorem** (necessity) | Paper 11 |
| Charge conservation (with symmetry as input) | **Theorem** (parity arithmetic of the support — no symmetry input) | Paper 15, Thm 8 |

Payment: one relation + one constant + one bit + two operating principles. Receipts: the five items above, plus the complex structure, the Z₄ phase, the gauge skeleton, and the coordinate map. The Born rule is a separate account (§8.1): both the fallback of importing it (even imported, the axiom count merely draws level with standard theory) and the upside of theoremization remain open.

## 7. The record of resolving, one by one, "what was not brought in comes out"

1. **Position cannot be read** → position = sign sector; resolved with zero additions to the state space (Paper 13).
2. **There are no complex numbers** → the choice of time direction selects the complex structure C⊂H; complex conjugation, too, appears internally as the inversion involution, not as an axiom (Paper 14).
3. **Where does the phase come from** → the Z₄ transport phase is derived (Paper 14); the necessity of orientation gauge-fixing emerged at the same time.
4. **The treatment of identical particles** → the canonical-convention theorem (the arbitrariness of a convention turns into a theorem; Paper 14).
5. **An obstacle to the transport principle itself** (audit finding) → resolved by making the conditions explicit and promoting them (collected in Paper 14).
6. **Touchdown onto coordinates** → Paper 15 constructs the map as a five-stage operation and machine-verifies it (round trip 500/500, conservation-commuting).

## 8. The two blank entries of the dictionary

### 8.1 The measure

The system supplies path amplitudes with no free parameters — the raw material of probability is complete. What is undecided is the aggregation rule, and this has been decomposed into a finite decision structure: **interference granularity × one event-structure bit × branch-channel conditioning** (Paper 16: the rigorous reduction of the selection problem, with its adjudication table — published; the adjudication experiments are in progress). The floor (fallback): import the Born rule as one axiom, and the axiom count draws level with standard theory — the lower bound of the connection is secured. The upside: if the record principle uniquifies the aggregation, a Born analogue becomes a theorem. The acceptance criterion is the triple — **convergence** (the selected measure merges into Born in the appropriate limit), **suppression** (in verified regimes, deviations below the statistical censorship ~1/√N), **prediction** (verifiable structure of deviations in the deep discrete regime).

### 8.2 Action and dynamics

Not constructed (as Paper 10 itself declares). The kinematics of uniform motion exists (Paper 15), but the derivation of "no force → uniform motion" and the dynamics of interaction do not.

## 9. Anatomy of the incomplete connections — why the remainder remains, and why its remaining is the key

We do not hide the incomplete points of the connection; we **dissect** them. First, a classification.

### 9.1 Two kinds

- **Under construction**: incomplete, but completable in principle (the isomorphism theorem; the construction of the continuum-limit map). A matter of time.
- **Structural mismatch**: points where complete agreement is **not to be expected in principle**. These are the subject of this section.

### 9.2 The anatomy table

| Locus | Content of the incompleteness | Structural cause | Inverted value (status as adjudication key) |
|---|---|---|---|
| **Measure (Born)** | The aggregation rule is unselected. The candidate measures (channel-coherent vs. configuration granularity) **differ rigorously** at finite s (e.g., the mirror parent with W=0 ∧ W′>0 — algebraic structure, not numerical error) | Standard theory possesses no selection mechanism either (Born is an axiom). This system contains a structure that **refuses to copy Born**: several candidates genuinely branch in the finite discrete regime, and at most one line can merge | The largest seat of this system's empirical content. If Born is an exact law, the upside closes (connection kept via the fallback). If Born is a **limiting effective law**, deviations of this system's type (granularity, parity, mirror-cancellation structure) live in the deep discrete regime — and only the side possessing a selection mechanism can specify the deviations' very **shape** |
| **Action / dynamics** | Not constructed | **Stopped at the premise check of the transport principle**: the action is formulated as an integral over external time t, but the dictionary's t entry reads "derived quantity" (no independent axis — a theorem). The premise's counterpart has a **different status**, so direct transport is impossible in principle. Dynamics can connect only through a reformulation into the language of record sequences | Standard theory itself is known to meet the **same difficulty** at the gravity interface (the problem of time: no external t). The reformulation this system is forced into — dynamics as the bookkeeping of records — may be isomorphic to the reformulation standard theory would be forced into at that interface. The point of failed connection coincides with standard theory's own point of failed extension |
| **Discrete ↔ continuous (transversal)** | Running through every row: this system is finite and discrete, standard theory continuous and infinite-dimensional. The connection extends to correspondence and inclusion; the limit map is unconstructed | Discreteness is an **axiom** (lattice, one bit), not a choice of approximation. Whether the system merges into the continuous side is a construction problem, not a premise | It poses the very adjudication problem "which is fundamental and which the approximation." If traces of discreteness (statistical censorship ~1/√N, granularity structure) are found at depth, the continuous description is the effective one; if not, this system is the effective model. **Adjudication is possible in both directions** |
| **Continuous symmetry** | The system's symmetries are B₄ (discrete) and Z₄ (derived phase). Connection to full U(1) and continuous rotations extends to the inclusion (Z₄↪U(1)) | A direct consequence of the discreteness axiom (a corollary of the transversal mismatch) | Z₄ quantization of phase and deep-regime anisotropy are structures that would not exist if continuous symmetry were exact — a seat of discernible structural prediction |
| **Dynamics of couplings** | There is a gauge skeleton, but no dynamics of coupling constants | Same root as the unconstructed action. Moreover, on the standard side coupling constants (e.g., α) are **measured inputs — the standard side has no theorem cell to fill at all** | The unique case where the absence of connection takes the form "**only this system's side can write the cell**." On the counting side a candidate identity exists (a separate series, published on Zenodo: an identity of the type α⁻¹=137+(π²/2)α, Concept DOI 10.5281/zenodo.19869266 — noted as an extra-series reference; a self-consistent equation at 8.7 ppb precision). For a quantity standard theory cannot compute in principle, the counting system holds a candidate formula |

### 9.3 The logic of inverted value (the thesis of this section)

Suppose the dictionary were filled completely. Then this system would be a **re-notation** of standard theory, and no experiment could distinguish the two — the academic value (a demonstrated minimality of axioms) would remain, but the empirical content would vanish. **Therefore, the only places where this system can "say" anything against standard theory are the points of incomplete connection.** The remainder (measure, action) and the transversal mismatch (discrete/continuous) are at once a table of deficiencies and the complete list of points where the two theories claim different things — that is, **a map of the places where observation may someday hand down a verdict**.

The direction of adjudication can be stated both ways. In a future where traces of discreteness, structured deviations from Born, or the need for a record-based dynamics are observed, this system's "incomplete connection" will turn out not to have been incomplete: it will stand as the prior description of the fact that **standard theory was the limiting effective law**. Conversely, in a future where exact Born, fundamental continuity, and external time carry everything, this system is demoted to an effective model or a re-notation. That a verdict stands either way — that is the reason for dissecting and publishing this incompleteness rather than hiding it.

### 9.4 Discipline (what this section does not claim)

This section does not claim the **reality** of deviations. It claims only the identification of the **structural locus** where deviations would live (granularity, parity, quantization, recordness). The identification turns into a "key" only when the "prediction" item of the triple criterion (§8.1) is quantified and translated into observable form — the work of the measure selection (Paper 16's adjudication table) and of the future dynamics line.

## 10. Conclusion

What was aimed at initially was a mere toy model. For the present state, inspect the dictionary of §5 — every major kinematic entry is filled, **every source is published**, on five items the system surpasses standard theory by holding as theorems what the latter installs as axioms, and the remainder has been pinned down to two entries: measure and action. As the author's conclusion: **we believe the system has connected, at a rather close level, to the basic requirements of standard quantum theory and of gauge theory.** Here "basic requirements" means the kinematic requirements (state space, observables, time, phase, symmetrization, conservation, measurement); on the quantum-theory side this covers that full range, on the gauge-theory side the skeleton (structural requirements) — the action remains unconstructed. And the incomplete points of the connection were dissected in §9: they are at once the table of deficiencies and the address book of this system's empirical content; future decidability lives there and only there. This conclusion is asserted only under the reading of §1: not as a new theory, but as an empirical study of the fact that **the structure of standard theory can re-emerge from assumptions this few and this alien**.

## 11. What is not claimed

(i) Re-derivation of standard theory (not chosen as a method). (ii) Physical identification (the central thesis: a different map with the same results). (iii) A new or unified theory (§1). (iv) Selection of the measure (in progress — this survey is measure-agnostic). (v) Dynamics. (vi) The reality of deviations (§9.4 — locus identification only). (vii) This survey contains no new numerical verification — every claim rests on the verifications of its sources (machine cross-checking of the citation map completed).

## 12. Limits and outlook

1. **Measure** (jurisdiction: Paper 16) — the decision structure is finite and the adjudication experiments live inside the model (in progress per the adjudication table: m=5@s31 and others). The quantification of §9's "prediction" item hinges here.
2. **Action** (jurisdiction: the future dynamics line) — reformulation of dynamics into the language of record sequences.
3. **The isomorphism theorem** (jurisdiction: the Paper 14 line or a standalone paper) — promoting the itemized dictionary into a single theorem: "full faithfulness of the comparison functor on the kinematic fragment." The ingredients (Pontryagin duality, C⊂H, Z₄↪U(1), B₄, the canonical convention) are complete; the proof is not.
4. **The continuum-limit map** (jurisdiction: a construction problem after the measure is settled) — the "under construction" part of §9's transversal mismatch.
5. **(Resolved) Publication order of sources** — Papers 13–16 were published on Zenodo on June 12, 2026, and every entry of §5 became "published." With this, the present document is the source-finalized edition (the course of the resolution is kept as a record).

## 13. Postscript: a retrospective of the origin — and on the criticism "isn't this just complex numbers?"

### 13.1 How it started (the author's retrospective)

In the earliest phase of this series, the starting assumptions were simple: the strong constraint that λ and ν are conjugate with λ=1/ν, and an energy-like conserved quantity, the conservation of a sum of squared frequencies Σνₙ² = N². That no structure appears without a displacement δ was clear from the start, so δ was introduced almost from the beginning. In the course of the work we noticed the possibility that ν and λ are orthogonal — that is, that they stand in the relation of sine and cosine waves — and further that δ=±1/2.

At that point we also noticed the possibility that **what we were handling was, in fact, nothing more than the polar-coordinate representation of complex-conjugate pairs of complex numbers, canonically normalized to amplitude 1**. In this series, however, we did not re-derive anything along that route. Rather, we believe that much of what was noticed was noticed precisely because we kept treating the objects, to the very end, as an amplitude-free (amplitude fixed at 1) ν and λ.

### 13.2 A reply to "isn't this just complex numbers written up plausibly?"

The criticism is to be expected, so we answer it in advance, in four points.

1. **The core of the observation is conceded and disclosed.** That the structure formed by amplitude-1 cos/sin pairs with phase corresponds locally to the unit circle U(1)⊂ℂ and its conjugate pair in polar representation — this is a fact, and the authors noticed it themselves mid-course. Far from hiding it, we state it here.

2. **But the direction of derivation is reversed.** Audit the axiom table (§3) — neither the complex field ℂ nor its conjugation operation **exists among the inputs**. The complex structure is an output: C⊂H appears as a theorem from the choice of time direction (§5), and conjugation appears internally as the inversion involution (§7-2). If "the result looks like complex numbers," that is the confirmation that **the complex structure re-emerged from one relation, one constant, and one bit** — not a failure of the series' claim but its very success case. The criticism bites only if ℂ has been smuggled into the inputs; the axiom table is short and auditable.

3. **The isomorphism is partial — too much is missing, and too much is left over, for "just complex numbers."** Missing: the amplitude degree of freedom (the system has none — fixed at 1); continuous phase (the phase appearing as state data is not continuous U(1) but the derived Z₄). The complex "plane" appears nowhere in the state space; complex line coefficients appear only in the readout of records, as derived objects. Left over: the lattice arithmetic of νλ=1; shell counting (Σνₙ²=N² is the tally of lattice points in frequency space, from which 1, 9, 137 emerge); the status of δ=±1/2 as a constant; the structures of parity, branches, and the readout hierarchy. The moment one says "complex numbers in polar form," all of these disappear into the shadow of the representation.

4. **The choice of representation was heuristically nontrivial — and the rigorous version of the criticism is the next theorem of this series.** Had we canonicalized early to "complex numbers in polar coordinates," the lattice arithmetic and counting structures above would most likely have been passed over as a "trivial rewriting." Even between equivalent representations, what can be seen depends on the representation. Keeping to the amplitude-free ν/λ representation to the end was the path that led to the discoveries of §5–§7. And if one tries to state "isn't it just a paraphrase of complex numbers?" **rigorously**, it becomes the task: "construct a structure-preserving correspondence between this system's kinematic fragment and a fragment of complex Hilbert space, and prove its full faithfulness" — that is, **precisely the isomorphism theorem of §12-3**. Rigorize the criticism and you obtain the theorem this series intends to prove next. In that sense this criticism is not an enemy but a roadmap.

## References (within the series)

**Policy**: this series deliberately cites no external literature. Its claims are correspondences at the structural level and include no claims of priority or doctrinal lineage (the mention of reconstruction programs in §1 is generic and relies on no specific work). All citations are to the series' published items and to the public repository.

**Published papers (Zenodo, Concept DOIs, all 16)**:

| # | Title | Concept DOI |
|---|---|---|
| 1 | Dual Geometry of Wavelength Space and Frequency Space | 10.5281/zenodo.20588036 |
| 2 | Radius Sweep of 4D Lattice Cell Counting | 10.5281/zenodo.20588038 |
| 3 | Closed 4D Structures and Lattice Counting (Addendum to Paper 1) | 10.5281/zenodo.20589261 |
| 4 | Reciprocal-Dual Cell Decomposition and the Hierarchical Vacuum | 10.5281/zenodo.20638962 |
| 5 | The Wave Dictionary and R Quantization | 10.5281/zenodo.20640454 |
| 6 | Conservation, Records, and Time | 10.5281/zenodo.20640456 |
| 7 | Configuration Statistics and Relation-Limited Readout | 10.5281/zenodo.20640458 |
| 8 | Two Ledgers, Internal Expansion, and the Area Law | 10.5281/zenodo.20640460 |
| 9 | Logic Waves, Censorship, and Stability | 10.5281/zenodo.20640462 |
| 10 | Emergent Gauge Structure | 10.5281/zenodo.20640464 |
| 11 | The Necessity of Dimension | 10.5281/zenodo.20640466 |
| 12 | The Minimal Twin-Space Model | 10.5281/zenodo.20640468 |
| 13 | Position Had Never Been Read | 10.5281/zenodo.20665633 |
| 14 | Time Direction and the Stage of Amplitudes | 10.5281/zenodo.20665661 |
| 15 | Projection onto xyztRQ | 10.5281/zenodo.20665688 |
| 16 | The Rigorous Reduction of the Measure Problem | 10.5281/zenodo.20665699 |

**Supplements and how to access them**: the supplements (84), all verification scripts, figure-generation code, the records of two-party verification exchanges, the review records of Papers 13–16 (eight reviews), and the drafts of this survey are in the series folder of the public repository. **Official URL (commit-pinned)**:

> https://github.com/WurabeSeiji/ai-chat-logs-open/tree/547af3c5c005fffc646ecf6ab513a200a3956ed2/%E6%B3%A2%E9%95%B7%E7%A9%BA%E9%96%93%E3%81%A8%E5%91%A8%E6%B3%A2%E6%95%B0%E7%A9%BA%E9%96%93%E3%81%AE%E5%8F%8C%E5%AF%BE%E5%B9%BE%E4%BD%95

The above is a reference pinned to commit 547af3c5 — **the immutable reference to the object this survey reports on: the record of Phase One** (a hash-pinned reference cannot, by construction, point to a version containing itself; this survey is therefore correctly absent from it). The immutable reference to this survey itself is this survey's own Zenodo DOI. For permanent citation of the repository as a whole, use the snapshot's **Concept DOI** below — stable across versions, it resolves also to subsequent releases that contain this survey. The **pinned reference to the state at completion of the publication work** (containing this survey v1.0, all sixteen papers, the indexes, and the inserted snapshot DOI; commit 0e96abb2):

> https://github.com/WurabeSeiji/ai-chat-logs-open/tree/0e96abb2ea995e318b0ccea9d940d230abbaabc3/%E6%B3%A2%E9%95%B7%E7%A9%BA%E9%96%93%E3%81%A8%E5%91%A8%E6%B3%A2%E6%95%B0%E7%A9%BA%E9%96%93%E3%81%AE%E5%8F%8C%E5%AF%BE%E5%B9%BE%E4%BD%95

Subsequent states are on the main branch of the same repository; the indexes are 論文一覧.md (all 16 DOIs, three reading routes) and 補遺一覧.md (the 84 supplements and the papers that absorb them). Supplements 1–41 are collected in Papers 5–12, 42–48 in Paper 13, 49–61 in Paper 14, 62–67 in Paper 15, 68–84 in Paper 16 (the measure adjudication experiments are in progress). **Repository snapshot Concept DOI**: 10.5281/zenodo.20666114 (this version: 10.5281/zenodo.20666115 — June 12, 2026, commit 6fae9b56; the complete Papers 1-16, Supplements 1-84, verification scripts and review records, 339 MB). Cite alongside the URL above for permanent reference.

**Extra-series reference (one item)**: the α identity of §9.2 (a counting identity of the type α⁻¹=137+(π²/2)α, published on Zenodo — Concept DOI 10.5281/zenodo.19869266).

## Appendix A: claims ↔ sources ↔ publication status ↔ verification methods

The table of §5 doubles as this appendix (all sources published; DOIs in the References). For the details of verification methods (exhaustive / sampled / random), follow the appendices of the respective sources. Machine cross-checking of the citation map was completed by Claude Code (June 12, 2026). The repository snapshot Concept DOI has been acquired and inserted (10.5281/zenodo.20666114).

**Acknowledgments / procedure**: this survey was produced under the two-party verification protocol of Claude Code (local machine verification) and claude.ai (independent re-computation, review, and drafting of this document). No physical identification is made.
