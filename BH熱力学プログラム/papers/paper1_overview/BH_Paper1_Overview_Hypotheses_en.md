# Black Hole Thermodynamics as Projection from a Higher Spacetime Structure: A Hypothesis-Driven Research Program (Paper 1: Overview and Hypotheses)

**Author**: Noriaki Kihara (WF System Co., Ltd.)
**ORCID**: [0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020)
**Date**: April 2026
**DOI**: [10.5281/zenodo.19839393](https://doi.org/10.5281/zenodo.19839393)

---

## Abstract

We propose a hypothesis-driven research programme that explores black hole thermodynamics from the standpoint of two conditional propositions: that black holes in an effective four-dimensional spacetime may be the projected images of higher-dimensional geometric structures (Hypothesis 1), and that, if spacetime is discrete at the Planck scale, then Hawking evaporation may correspond to the drift of discrete units across a higher-dimensional region's boundary (Hypothesis 2). Both hypotheses follow from three premises widely shared in the field: that black hole entropy is geometric, that evaporation is local to the horizon, and that spacetime may be discrete at the Planck scale. This first paper situates the programme within the existing literature, articulates the two hypotheses precisely, previews the structure of the verification across five companion papers, and addresses anticipated critiques. The paper does not derive black hole thermodynamics; it defines a research programme.

---

## §1. Introduction

The thermodynamic description of black holes, established by Bekenstein, Hawking, and Bardeen–Carter–Hawking in the early 1970s, remains one of the most striking and least fully understood developments in classical and semi-classical gravity. A black hole, viewed externally, behaves as a thermal object: it carries entropy proportional to the area of its horizon, radiates with a characteristic temperature inversely proportional to its mass, and evaporates over a finite lifetime. These features are derived within the framework of quantum field theory on a fixed curved background, and their underlying microscopic origin — what microstates count toward the Bekenstein–Hawking entropy, and what mechanism is responsible for evaporation at the Planck scale — is a problem that connects gravity, quantum theory, and information in ways that are still being negotiated by the field.

Several research programs address this problem. The holographic principle of 't~Hooft (1993) and Susskind (1995), realized concretely in the AdS/CFT correspondence of Maldacena (1997), proposes that the bulk degrees of freedom in a gravitational theory are encoded on its boundary. Loop quantum gravity, causal dynamical triangulations, and causal set theory each posit, in different ways, that spacetime is fundamentally discrete at the Planck scale, and seek to derive black hole thermodynamics from this discreteness. Higher-dimensional generalizations — the Tangherlini black holes (1963) and the Myers–Perry rotating solutions (1986) — extend the four-dimensional formulae to arbitrary spacetime dimension, and have become a standard testbed for holographic and string-theoretic ideas. Each of these programs offers partial insight; none has yet produced a complete and observationally testable account.

This paper proposes a hypothesis-driven research program that explores a complementary direction. We begin from three premises that we take to be widely shared within the field — that black hole entropy is geometric (Bekenstein), that evaporation reflects horizon-scale dynamics (Hawking), and that spacetime may be discrete at the Planck scale (the common premise of the discrete quantum gravity programs). From these we extract two hypotheses. The first is that a black hole, when viewed within an effective four-dimensional spacetime, may be the projected image of a higher-dimensional geometric structure rather than an independent four-dimensional entity. The second is that, if spacetime is built from discrete minimal units, the drift of these units across the boundary of a region in the higher-dimensional structure may correspond to what is observed as Hawking evaporation in the four-dimensional projection.

We emphasize at the outset that this paper does not derive black hole thermodynamics. It defines a research program. The mathematical and physical content — the metric arising under central projection, the combinatorial structure of unit-cube packings inside a four-dimensional ball, the asymptotic correspondence with five-dimensional Tangherlini entropy, and the self-consistent dynamical determination of the curvature radius — is developed in five companion papers, each self-contained and each making a single, narrow claim. The present paper provides the connective tissue: it situates the program within the existing literature, articulates the two hypotheses precisely enough to be tested, and previews the structure of the verification across the full sequence.

Three features distinguish the program. First, it does not assume the holographic principle or the AdS/CFT correspondence. The four-plus-one-dimensional structure here arises from a single elementary construction — central (gnomonic) projection from a four-sphere of radius $R$ embedded in $\mathbb{R}^5$ to a four-dimensional subjective coordinate chart — and the Schwarzschild-like behaviour of the induced metric is a consequence of that construction, not an input. Second, the discrete-spacetime side of the program does not posit any specific microscopic theory of quantum gravity. It uses only the requirement that unit cubes be packed inside a four-dimensional ball without protrusion, and the resulting integer-lattice combinatorics. Third, the natural comparison object for the resulting thermodynamics is the four-plus-one-dimensional Tangherlini black hole, not the familiar 3+1-dimensional Schwarzschild black hole. This is a deliberate strategic choice: it positions the program within the higher-dimensional gravity literature (Randall–Sundrum, ADD models, Myers–Perry) where the formal machinery of higher-dimensional black holes is already well-developed.

The contribution of this first paper is therefore modest. It does not claim to resolve any open problem in black hole thermodynamics; it does not predict new observational signatures; and it does not displace any existing program. It articulates a hypothetical framework with the property that, if the two hypotheses turn out to be useful, the program produces a quantitative correspondence with established Tangherlini results from a minimal set of geometric and combinatorial assumptions. Whether the framework is useful is a question that the companion papers attempt to answer; this paper only frames the question.

The remainder is organized as follows. Section 2 reviews the established results in black hole thermodynamics that the program will eventually need to recover or contact. Section 3 states the two hypotheses precisely. Section 4 outlines the verification program, identifying which paper contributes to which hypothesis and at what level of rigor. Section 5 anticipates the principal critiques the program is likely to face and addresses them in advance. Section 6 concludes.


## §2. Review of Existing Black Hole Thermodynamics

This section reviews the established results in black hole thermodynamics that the verification program will eventually need to either recover or contact. The review is intentionally compact and neutral. It is organized along the principal lines of historical development, and selected to highlight precisely those results that bear on the hypotheses introduced in §3.

### §2.1 The Schwarzschild metric and the classical description

The simplest non-trivial vacuum solution of the Einstein field equations is the Schwarzschild solution, derived in 1916. In standard $(t, r, \theta, \phi)$ coordinates,
$$ds^2 = -f(r)\, dt^2 + f(r)^{-1}\, dr^2 + r^2 \, d\Omega_2^2, \qquad f(r) = 1 - \frac{2GM}{r}.$$
The function $f(r)$ vanishes at $r_s = 2GM$, defining the Schwarzschild horizon. This horizon plays no special role in the classical mechanics of test particles — geodesics traverse it smoothly in the appropriate coordinate frames — but it acquires central importance in the thermodynamic and quantum-mechanical descriptions reviewed below.

### §2.2 Hawking radiation: thermal evaporation

In 1975, Hawking applied the formalism of quantum field theory in curved spacetime to the Schwarzschild background and derived a striking result: a Schwarzschild black hole emits thermal radiation at temperature
$$T_H = \frac{\hbar \kappa}{2\pi k_B c} = \frac{\hbar c^3}{8\pi G M k_B},$$
where $\kappa$ is the surface gravity at the horizon. The radiation has a Planckian spectrum and carries energy away from the black hole, leading to a finite evaporation lifetime $\tau \sim G^2 M^3 / (\hbar c^4)$.

The temperature is inversely proportional to mass: smaller black holes are hotter and evaporate faster. The endpoint of evaporation — the fate of a Planck-mass black hole — remains an open problem.

### §2.3 Bekenstein entropy

In 1973, Bekenstein argued on thermodynamic grounds that a black hole must possess an entropy proportional to the area of its horizon. Combining his proposal with Hawking's temperature gives the Bekenstein–Hawking entropy
$$S_{BH} = \frac{k_B c^3 A}{4 G \hbar} = \frac{A}{4 \ell_P^2},$$
where $A = 4\pi r_s^2$ is the area of the Schwarzschild horizon and $\ell_P = (G\hbar/c^3)^{1/2}$ is the Planck length. The entropy scales as the square of the mass.

The microscopic origin of $S_{BH}$ — what counts as a black hole microstate — is the principal open problem in the field. The number of microstates required is enormous: a solar-mass black hole has $S_{BH}/k_B \sim 10^{77}$.

### §2.4 The four laws of black hole thermodynamics

Bardeen, Carter, and Hawking (1973) established a complete thermodynamic analogy:

- **Zeroth law.** Surface gravity $\kappa$ is constant on a stationary horizon (analogue of temperature uniformity in equilibrium).
- **First law.** $dM = (\kappa / 8\pi G)\, dA + \Omega\, dJ + \Phi\, dQ$ (analogue of $dU = T\, dS - p\, dV$).
- **Second law.** Horizon area is non-decreasing in classical processes; in quantum processes, the generalized entropy $S_{\mathrm{gen}} = S_{BH} + S_{\mathrm{matter}}$ is non-decreasing.
- **Third law.** A horizon with $\kappa = 0$ (extremal) cannot be reached by any finite physical process.

These laws are presently understood as a robust thermodynamic structure, valid in classical and semi-classical settings.

### §2.5 The holographic principle

't Hooft (1993) and Susskind (1995) proposed the holographic principle: the number of independent quantum-mechanical degrees of freedom inside a region of space is bounded by the area of the boundary in Planck units, not by the volume. The principle is motivated by the Bekenstein–Hawking entropy: if a black hole's entropy is governed by its horizon area, then the maximum entropy that can fit inside any region is also area-bounded.

The principle suggests that a fundamental description of gravitational physics in $D$ spacetime dimensions has degrees of freedom living on a $(D-1)$-dimensional surface. It is the conceptual root of the AdS/CFT correspondence reviewed below, but it is logically more general.

### §2.6 The AdS/CFT correspondence

Maldacena (1997) proposed that Type IIB string theory on $\mathrm{AdS}_5 \times S^5$ is equivalent to $\mathcal{N} = 4$ super Yang–Mills theory on the boundary of $\mathrm{AdS}_5$. The correspondence is a concrete realization of the holographic principle: bulk gravitational physics is encoded in a boundary gauge theory.

In the context of black hole thermodynamics, AdS/CFT gives explicit microstate counts (in the BPS sector) and dual descriptions of black hole horizons in terms of thermal CFT states. The correspondence is the most successful current implementation of the holographic principle.

### §2.7 Higher-dimensional black holes

Tangherlini (1963) generalized the Schwarzschild solution to arbitrary spacetime dimension $D$. In $D = d + 1$ dimensions, the Schwarzschild–Tangherlini metric is
$$ds^2 = -f_T(r)\, dt^2 + f_T(r)^{-1}\, dr^2 + r^2\, d\Omega_{d-1}^2, \qquad f_T(r) = 1 - \frac{\mu}{r^{d-2}},$$
with $\mu$ proportional to mass. The horizon is at $r_h = \mu^{1/(d-2)}$, the surface area of the horizon is the volume of $S^{d-1}(r_h)$, and the Bekenstein–Hawking entropy formula generalizes to $S_{BH} = A_{d-1}/(4\ell_P^{d-1})$.

For our purposes the relevant case is $D = 5$ ($d = 4$), where
$$f_T(r) = 1 - \frac{8GM}{3\pi r^2}, \qquad r_h^2 = \frac{8GM}{3\pi}, \qquad T_H = \frac{1}{2\pi r_h}, \qquad S_{BH} = \frac{2\pi^2 r_h^3}{4 G_5}.$$
This is the natural comparison object for the program developed here.

Myers and Perry (1986) generalized the solution to include rotation; the Myers–Perry solutions are the higher-dimensional analogues of Kerr.

### §2.8 Quantum gravity programs

Several programs treat spacetime as fundamentally discrete at the Planck scale.

*Loop quantum gravity* (Rovelli, Smolin, and others, 1990s) constructs a Hilbert space for canonical quantum gravity on which area and volume operators have discrete spectra with a minimum unit at the Planck scale. The Bekenstein–Hawking entropy is recovered in this framework via the counting of horizon spin-network states.

*Causal dynamical triangulations* (Ambjørn, Jurkiewicz, Loll, 1998 onwards) implements a non-perturbative path integral over discrete causal structures, recovering an emergent four-dimensional spacetime in continuum simulations.

*Causal set theory* (Bombelli, Lee, Meyer, Sorkin, 1987) replaces the smooth spacetime manifold with a discrete partial order representing causal relations between Planck-scale events.

These programs share Premise C of §3.1 (Planck-scale discreteness) but differ widely in their specific implementations. The verification program of this paper does not commit to any particular implementation; it requires only the abstract premise that discreteness is meaningful at the relevant scale.

### §2.9 Synthesis

The established body of black hole thermodynamics has the following structural features:

- An entropy proportional to a geometric quantity (Premise A).
- A semi-classical evaporation mechanism whose origin is local to the horizon (Premise B).
- A range of programs invoking Planck-scale discreteness as a fundamental ingredient (Premise C).

These three premises are common ground; the differences among the programs concern the specific implementation of each. The verification program of this paper extracts these common premises and pursues a particular geometric and combinatorial implementation — central projection from a higher-dimensional structure, and integer-lattice packing of unit cubes — that has not, to our knowledge, been pursued in the literature.


## §3. Two Hypotheses

We now state the two hypotheses around which the verification program is organized. Both are formulated as conditional propositions; neither is claimed to be true on physical grounds at this stage.

### §3.1 Premises

We start from three premises which, we contend, are shared by a substantial fraction of the contemporary literature on black hole thermodynamics and quantum gravity:

**Premise A** *(Geometric origin of entropy).* The entropy of a black hole is determined by a geometric quantity associated with its horizon. In the standard formulation this is the horizon area in Planck units; we leave the precise geometric quantity unspecified at this stage.

**Premise B** *(Horizon-scale origin of evaporation).* Hawking radiation reflects dynamics at the scale of the horizon, not at infinity. Whatever microscopic mechanism underlies evaporation, it is local to the near-horizon region.

**Premise C** *(Planck-scale discreteness).* Spacetime may be effectively discrete at the Planck scale. We do not specify the form of this discreteness; we only require that, at sufficiently fine resolution, geometric quantities admit a description in terms of integer-valued combinatorial data.

Premise A is the consensus reading of the Bekenstein–Hawking formula. Premise B is the operational content of Hawking's 1975 derivation: the relevant quantum field theoretic computation is performed in the near-horizon region, and the resulting thermal spectrum at infinity reflects the redshift of these near-horizon modes. Premise C is the common point of departure for loop quantum gravity, causal dynamical triangulations, and causal set theory; we adopt it without committing to any particular implementation.

### §3.2 Hypothesis 1: The Projection Hypothesis

**Hypothesis 1.** *A black hole, as it appears within an effective four-dimensional spacetime, is the projected image of a geometric structure embedded in a higher-dimensional spacetime. The thermodynamic quantities associated with the black hole — its entropy, temperature, and evaporation rate — are properties of the projection rather than of an autonomous four-dimensional object.*

The hypothesis is consistent with Premise A. If the four-dimensional horizon is the projection of a higher-dimensional geometric object, then "horizon area" in four dimensions corresponds to some geometric quantity of that higher-dimensional object, and the Bekenstein-style relation is naturally interpreted as a relation between four-dimensional area and the higher-dimensional geometric quantity.

The hypothesis does **not** specify the higher-dimensional structure. In the verification program, we work with a particular concrete realization: gnomonic projection from a four-sphere $S^4 \subset \mathbb{R}^5$ of radius $R$ onto a four-dimensional subjective coordinate chart. This is one possible realization; the hypothesis as stated is more general.

The hypothesis is **not** equivalent to the holographic principle. Holography places the boundary at lower dimension; the projection hypothesis places the bulk at *higher* dimension and views four-dimensional physics as the projected image. The two are logically independent: one could entertain holography without projection, projection without holography, or both, or neither.

### §3.3 Hypothesis 2: The Discrete Drift Hypothesis

**Hypothesis 2.** *If spacetime is composed of discrete minimal units, and if a black hole is the projection of a finite higher-dimensional region, then Hawking evaporation corresponds to the drift of discrete units out of that region. Specifically, the difference between the continuous volume of the higher-dimensional region and the integer count of discrete units that fit inside it determines the entropy and the evaporation dynamics observable in the four-dimensional projection.*

This is consistent with Premises B and C. Premise B requires that the dynamics be local to the horizon: the discrete drift hypothesis localizes the relevant dynamics to the boundary of the higher-dimensional region. Premise C requires discreteness at the Planck scale: the hypothesis takes this discreteness seriously and asks what counting problem it gives rise to.

The hypothesis is **not** equivalent to any form of fuzzball or microstate-counting program. It does not posit that black holes are made of strings or branes; it does not invoke AdS/CFT to count boundary microstates. The "discrete units" are taken to be unit cubes in an integer lattice, and the counting problem is a four-dimensional packing problem — well-defined classically and combinatorially — whose connection to thermodynamics is the content of the program.

### §3.4 Status of the hypotheses

Both hypotheses are stated as conditional propositions. They make no claim about the physical reality of the higher-dimensional structure or the discrete units; they ask only whether, *if* one assumes such structures, one can derive the established thermodynamic quantities and at what level of agreement.

In the companion papers we show that the answer is partially affirmative:

- The Schwarzschild-like behaviour of the induced metric under central projection (Paper 2) supports Hypothesis 1 at the level of the metric structure.
- The asymptotic packing analysis (Paper 3) gives the volume deficit $\Delta(R) \sim c\,R^3$ with $c = 8/(3\pi)$, providing a concrete quantity to be matched.
- The correspondence with the four-plus-one-dimensional Tangherlini Bekenstein–Hawking entropy (Paper 4) realizes the qualitative scaling of Hypothesis 2 quantitatively.
- The dynamical self-consistency analysis (Paper 5) recovers the standard Tangherlini scalings $R \propto M^{1/2}$, $T_H \propto 1/\sqrt{M}$, and lifetime $\propto M^2$, with one quantitative discrepancy in the evaporation-rate exponent that we mark as an open problem.

The verification is partial; the discrepancies are explicit. The program does not claim to have solved the problem of black hole thermodynamics. It claims only that the two hypotheses, taken together with the minimal assumptions of central projection and discreteness, produce a calculation that contacts the established results in a non-trivial way.


## §4. The Verification Program: Overview of the Six Papers

This section sketches the verification program, identifying which paper contributes to which hypothesis and at what level of rigor. Each companion paper is self-contained and may be read independently of the others. The collective structure outlined here is descriptive: it organizes the program but is not a logical prerequisite for any individual paper.

### §4.1 Mapping of contributions

The two hypotheses of §3 are addressed across the program as follows.

**Hypothesis 1 (Projection Hypothesis):**

| Paper | Type of contribution | Level of rigor |
|---|---|---|
| 2 | Mathematical foundation: gnomonic projection metric, Schwarzschild-like form | Rigorous, theorem-proof |
| 4 | Quantitative correspondence: $\Delta(R)$ vs. four-plus-one-dimensional Tangherlini entropy | Quantitative, with explicit numerical match |
| 5 | Dynamical consistency: scaling of $R$ with mass | Quantitative, with one open exponent |
| 6 | Synthesis or open problems | (Character to be determined) |

**Hypothesis 2 (Discrete Drift Hypothesis):**

| Paper | Type of contribution | Level of rigor |
|---|---|---|
| 3 | Mathematical foundation: lattice packing, asymptotic constant $c = 8/(3\pi)$ | Rigorous, theorem-proof |
| 4 | Identification: $\Delta(R)$ as a candidate for the entropy quantity | Conjectural |
| 5 | Dynamical role: drift mechanism in the curvature radius evolution | Heuristic |
| 6 | Synthesis or open problems | (Character to be determined) |

### §4.2 Logical dependency

The companion papers are designed to be readable in isolation. Paper 2 is purely about the metric structure under central projection; Paper 3 is purely combinatorial; Paper 4 connects the two via dimensional analysis and a specific numerical correspondence; Paper 5 introduces dynamics. The sequence

$$\text{Papers 2, 3} \quad\rightarrow\quad \text{Paper 4} \quad\rightarrow\quad \text{Paper 5} \quad\rightarrow\quad \text{Paper 6}$$

is the natural reading order, but each paper is intelligible without its predecessors.

### §4.3 Domain of validity for each hypothesis

**Hypothesis 1 (Projection):** is supported at the level of the metric — the central projection construction produces a metric with Schwarzschild-like behaviour and a horizon structure. This is a mathematical statement about a particular geometric construction; it does not by itself establish that this construction describes physical black holes.

**Hypothesis 2 (Discrete Drift):** is supported at the level of scaling — the volume deficit $\Delta(R)$ has the right $R^3$ scaling to match four-plus-one-dimensional Tangherlini entropy. The interpretation of $\Delta(R)$ as the "discrete drift entropy" is a hypothesis about how the mathematical quantity should be understood physically. This interpretation is partially supported by the consistency of the dynamical scaling laws derived in Paper 5, but the precise mechanism by which discrete units drift across the boundary is not derived.

### §4.4 The evaporation rate exponent

The principal residual discrepancy in the program is the evaporation rate exponent. The simplest dynamical assumption — Stefan–Boltzmann-like radiation from the horizon area at the Hawking temperature — gives
$$\frac{dM}{dt} \propto -M^{n}, \qquad n = -1/2 \text{ or } 3/2 \text{ (depending on convention)},$$
whereas the standard semi-classical Hawking computation gives $n = -1$. The discrepancy is an integer offset in a specific exponent; it is not a sign error or a missing factor of $\pi$.

We mark this discrepancy as an open problem in the program. Paper 5 discusses possible resolutions, including:

- A more careful treatment of the back-reaction of the radiation on the geometry.
- Modifications to the discrete drift mechanism beyond the simplest scaling argument.
- A re-interpretation of what is meant by "evaporation" in the projected four-dimensional theory.

The program does not claim to have a preferred resolution at this stage. Paper 6 may either propose a resolution or leave the discrepancy as a precisely identified open problem.

### §4.5 What the program does and does not establish

By the end of the verification program, the following will have been established:

- The mathematical consistency of the central projection construction (Paper 2).
- The exact value of the asymptotic constant $c = 8/(3\pi)$ in the volume deficit (Paper 3).
- The correspondence in $R^3$ scaling between $\Delta(R)$ and four-plus-one-dimensional Bekenstein–Hawking entropy (Paper 4).
- The recovery of the Tangherlini scaling laws $R \propto M^{1/2}$, $T_H \propto 1/\sqrt{M}$, lifetime $\propto M^2$ (Paper 5).

The following will **not** be established:

- A first-principles derivation of the Bekenstein–Hawking entropy from the discrete drift mechanism (only the scaling and leading constant are matched).
- A complete agreement with the standard Hawking computation (the evaporation-rate exponent is left as an open problem).
- A resolution of the relation to the holographic principle, AdS/CFT, or the standard quantum gravity programs.
- A description of physical 3+1-dimensional black holes (the natural object is 4+1-dimensional).

The program is best understood as a *partial verification of a particular geometric and combinatorial framework* against the well-established Tangherlini results. The successes are taken to indicate that the framework is internally consistent and merits further development; the residual discrepancies are taken to indicate where the framework, as currently formulated, falls short.


## §5. Anticipated Limitations and Critiques

This section addresses the principal critiques that we expect the program to encounter. We state each critique in the form a careful reader might raise, and offer the program's response.

### §5.1 Why a four-plus-one-dimensional formulation?

**Critique.** The observed universe has three spatial dimensions and one temporal dimension. Black holes that have been theoretically studied in detail and whose thermodynamic properties are anchored in the literature are 3+1-dimensional. Why does this program take the natural comparison object to be a 4+1-dimensional Tangherlini black hole?

**Response.** The 4+1-dimensional formulation is a consequence of the geometric construction adopted in Paper 2: gnomonic projection from a four-sphere $S^4 \subset \mathbb{R}^5$ onto a four-dimensional subjective chart. The induced metric on the four-dimensional chart, equipped with the curvature radius $R$ as an additional independent geometric parameter, gives rise naturally to a five-dimensional structure $(x, y, z, t, R)$. The black hole metric that arises in this construction is, on dimensional grounds, the Tangherlini analogue rather than the Schwarzschild form.

We do not claim that 4+1 is the correct dimensionality of physical spacetime. We claim that, *given the geometric construction of Paper 2*, the natural comparison is with the 4+1 case, and that the agreement with Tangherlini thermodynamics serves as a non-trivial check on the construction itself. The relation to observed 3+1-dimensional physics — whether through a Randall–Sundrum-style brane localization, an ADD-style large extra dimension, or some other mechanism — is left for future work and is not within the scope of this program.

This positioning is shared with a substantial fraction of the higher-dimensional gravity literature: Tangherlini, Myers–Perry, Randall–Sundrum, and ADD all study generalizations whose physical realization in 3+1 dimensions is a separate question from their internal consistency.

### §5.2 Why not assume the holographic principle?

**Critique.** The holographic principle and AdS/CFT have produced a great deal of progress in black hole thermodynamics, including detailed entropy counts and explicit dynamical descriptions. Why does this program not make use of them?

**Response.** The program is constructed as an independent path of inquiry. Its starting assumptions — central projection from a higher-dimensional spacetime and discreteness at the Planck scale — are logically independent of the holographic principle. Whether the resulting picture is compatible with, supports, contradicts, or is orthogonal to holography is a question that can be examined only after the independent path has been worked out. To assume holography at the outset would obscure that question.

We emphasize that nothing in the program relies on the failure of the holographic principle. If holography is correct, the agreement with Tangherlini thermodynamics derived in Papers 4 and 5 must, by some mechanism, be consistent with the holographic counting; the precise correspondence is an open problem that the program does not attempt to solve.

### §5.3 Why not derive the results fully and quantitatively?

**Critique.** Paper 1 only states hypotheses; the quantitative derivations are deferred to companion papers. The companion papers themselves leave the evaporation-rate exponent ($n = 3/2$ in the simplest treatment, $n = -1$ in the standard Hawking computation) unreconciled. Is the program offering anything beyond a sketch?

**Response.** Yes. The program offers (i) a precise mathematical framework — the central projection construction and the integer-lattice packing problem — that admits rigorous analysis in its own right; (ii) the explicit asymptotic constant $c = 8/(3\pi)$ for the volume deficit, derived analytically and confirmed numerically; (iii) the matching of the $R^3$ scaling and the leading numerical coefficient with the 4+1-dimensional Bekenstein–Hawking entropy; and (iv) the recovery of the standard Tangherlini scaling laws $R \propto M^{1/2}$, $T_H \propto 1/\sqrt{M}$, and lifetime $\propto M^2$.

The unreconciled evaporation-rate exponent is a genuine open problem. The program is explicit that the exponent obtained from the simplest dynamical assumption ($n = 3/2$) does not match the standard semi-classical result ($n = -1$). We mark this as an open problem rather than concealing it, and we conjecture that a more careful treatment of the back-reaction or of the discrete drift mechanism is required to resolve it. This honesty about discrepancies is part of the methodology, not a defect.

### §5.4 Why is the program testable?

**Critique.** The hypotheses are stated as conditional propositions and the program does not predict new observations. In what sense is it scientific?

**Response.** The program is testable in three distinct senses:

1. *Mathematical consistency.* Papers 2 and 3 are entirely mathematical: any error in the central projection metric, the asymptotic constant $c$, or the Jacobi-style closed-form representations of $N(k)$ would falsify the foundational claim. The numerical computation is reproducible from the definition; the analytic derivation is checkable line by line.

2. *Quantitative correspondence.* Papers 4 and 5 produce specific numerical relations — the constant $c = 8/(3\pi)$, the ratio $\Delta/S_{BH} \to 32/(3\pi) \approx 3.40$, the Tangherlini scalings $R \propto M^{1/2}$ — that are either matched or not. They are matched, with the noted exception of the evaporation-rate exponent. Failure of these correspondences would falsify the program's main claim of contact with established thermodynamics.

3. *Open-problem identification.* Paper 6 formalizes the residual discrepancies. The program identifies the evaporation-rate exponent as the principal locus of disagreement with the standard semi-classical computation. This is itself a useful research output: it specifies what would have to be modified, in either the program or the standard computation, to achieve full consistency.

The program does not predict new observable signatures. We do not claim it as an empirical contribution to physics; we claim it as a mathematical and conceptual exploration of what follows from a particular set of geometric assumptions. Whether the exploration turns out to have empirical content is a question for future work that lies outside the scope of the present program.

### §5.5 What is *not* claimed

To prevent misreading, we list explicitly what the program does **not** claim:

- It does not claim that physical spacetime is four-plus-one-dimensional.
- It does not claim that the central projection construction is the unique or correct geometric description of gravity.
- It does not claim to have derived the Bekenstein–Hawking entropy from first principles.
- It does not claim that the evaporation-rate exponent is $3/2$; it identifies the discrepancy with the standard $n = -1$ as an open problem.
- It does not claim to compete with, displace, or refute the holographic, loop, causal-set, or string-theoretic programs.

What it claims is the conjunction of two conditional propositions (Hypotheses 1 and 2), the mathematical structures that follow from them, and the partial quantitative correspondence with established Tangherlini results.


## §6. Conclusion

We have proposed a hypothesis-driven research program that explores black hole thermodynamics from the standpoint of two conditional propositions. The first is that black holes, as seen in an effective four-dimensional spacetime, may be the projected images of higher-dimensional geometric structures. The second is that, if spacetime is discrete at the Planck scale, then Hawking evaporation may correspond to the drift of discrete units across the boundary of a higher-dimensional region.

Both hypotheses follow from three premises that are widely shared in contemporary work on black hole thermodynamics and quantum gravity: that black hole entropy is geometric, that evaporation is local to the horizon, and that spacetime may be discrete at the Planck scale. Neither hypothesis displaces or contradicts any established result. They specify a particular path of inquiry — central projection from a higher-dimensional structure together with integer-lattice combinatorics — that has not, to our knowledge, been pursued in the existing literature.

The verification program comprises five companion papers, each self-contained and each making a single narrow claim. By the program's end, the consistency of the central projection construction (Paper 2), the precise value of the asymptotic constant $c = 8/(3\pi)$ for the volume deficit (Paper 3), the correspondence with four-plus-one-dimensional Tangherlini entropy (Paper 4), and the recovery of the Tangherlini scaling laws (Paper 5) are established. The principal residual discrepancy is the evaporation rate exponent, which we leave as a precisely identified open problem.

The program does not solve any open problem in black hole thermodynamics. It establishes that, from minimal geometric and combinatorial assumptions, a non-trivial quantitative correspondence with established Tangherlini results can be derived. Whether this correspondence indicates something physical, or is merely a mathematical curiosity, is a question that the program does not attempt to answer.

The program's value, if any, is in its self-completeness and in the precision of its open problems. The mathematical content of Papers 2 and 3 is fully derivable from the stated assumptions. The numerical correspondences of Papers 4 and 5 are explicit and either match or fail to match. The residual discrepancy is identified at a specific quantitative level. None of these features depends on holography, AdS/CFT, or any specific quantum gravity program. The program stands or falls on its own internal consistency and on its agreement with the well-established results to which it directly compares.

We close with three remarks.

First, the program's positioning at four-plus-one dimensions, rather than three-plus-one, is a deliberate strategic choice motivated by the geometric construction of Paper 2. We make no claim about the physical dimensionality of spacetime; we work in four-plus-one because that is the natural setting for the construction we adopt.

Second, the absence of any holographic input is also deliberate. The program is designed to be informative independently of whether holography is correct. If holography is correct, the agreement with Tangherlini results derived here must be reconcilable with holographic counts; that reconciliation is not within the program's scope.

Third, the methodology — stating premises explicitly, deriving consequences rigorously, identifying discrepancies precisely — is intended to make the program useful regardless of whether one accepts its hypotheses. A reader who rejects the hypotheses still obtains a self-contained mathematical and combinatorial framework whose outputs may have interest in their own right.

The development of the program in the companion papers begins, in Paper 2, with the geometric construction.

---

## References

1. K. Schwarzschild, "Über das Gravitationsfeld eines Massenpunktes nach der Einsteinschen Theorie," *Sitzungsber. Preuss. Akad. Wiss.*, 189–196 (1916).
2. J. D. Bekenstein, "Black Holes and Entropy," *Phys. Rev. D* **7**, 2333 (1973).
3. J. M. Bardeen, B. Carter, S. W. Hawking, "The Four Laws of Black Hole Mechanics," *Comm. Math. Phys.* **31**, 161 (1973).
4. S. W. Hawking, "Particle Creation by Black Holes," *Comm. Math. Phys.* **43**, 199 (1975).
5. G. 't Hooft, "Dimensional Reduction in Quantum Gravity," arXiv:gr-qc/9310026 (1993).
6. L. Susskind, "The World as a Hologram," *J. Math. Phys.* **36**, 6377 (1995).
7. J. Maldacena, "The Large-$N$ Limit of Superconformal Field Theories and Supergravity," *Adv. Theor. Math. Phys.* **2**, 231 (1998).
8. F. R. Tangherlini, "Schwarzschild Field in $n$ Dimensions and the Dimensionality of Space Problem," *Nuovo Cimento* **27**, 636 (1963).
9. R. C. Myers, M. J. Perry, "Black Holes in Higher Dimensional Spacetimes," *Annals Phys.* **172**, 304 (1986).
10. C. Rovelli, *Quantum Gravity*, Cambridge University Press (2004).
11. J. Ambjørn, J. Jurkiewicz, R. Loll, "Emergence of a 4D World from Causal Quantum Gravity," *Phys. Rev. Lett.* **93**, 131301 (2004).
12. L. Bombelli, J. Lee, D. Meyer, R. D. Sorkin, "Space-time as a Causal Set," *Phys. Rev. Lett.* **59**, 521 (1987).
13. J.-L. Lagrange, "Démonstration d'un théorème d'arithmétique," *Nouv. Mém. Acad. Roy. Sci. Belles-Lettres Berlin*, 123–133 (1770).
14. C. G. J. Jacobi, *Fundamenta Nova Theoriae Functionum Ellipticarum*, Königsberg (1829).
15. L. Randall, R. Sundrum, "A Large Mass Hierarchy from a Small Extra Dimension," *Phys. Rev. Lett.* **83**, 3370 (1999).
16. N. Arkani-Hamed, S. Dimopoulos, G. Dvali, "The Hierarchy Problem and New Dimensions at a Millimeter," *Phys. Lett. B* **429**, 263 (1998).
