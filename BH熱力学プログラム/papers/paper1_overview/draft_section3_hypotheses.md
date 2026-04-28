# 論文1 §3 仮説の提示 試案

---

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

---

**Style notes**:
- 仮説を「conditional proposition」として提示
- ホログラフィー、ファズボール等の既存プログラムとの違いを明示
- 物理的真理を主張せず、条件付きの探求として記述
- 後続論文の貢献を整理した形で予告
