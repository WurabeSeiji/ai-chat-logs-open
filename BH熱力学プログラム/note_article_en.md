# Black Holes as Projections from a Higher Spacetime Structure — A Six-Paper Research Programme

## Purpose of This Article

This article introduces a six-paper research programme that reconstructs black hole thermodynamics — Bekenstein–Hawking entropy, Hawking radiation, evaporation lifetime — from two minimal assumptions only: **central projection and Planck-scale discreteness**.

We do **not** assume the holographic principle ('t Hooft, Susskind 1993–95) or the AdS/CFT correspondence (Maldacena 1997). The programme is an independent thought experiment: "If spacetime is the projected image of a higher-dimensional structure and is discrete at the Planck scale, what follows?"

The result: the principal thermodynamic results of a four-plus-one-dimensional Schwarzschild–Tangherlini black hole — the $R^3$ scaling of entropy, the temperature $T_H \propto M^{-1/2}$, the lifetime $\tau \propto M^{3/2}$ — are reproduced almost exactly.

All six papers are published as preprints on Zenodo (PDF, LaTeX, and Markdown in both English and Japanese, CC BY 4.0).

GitHub: https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/BH熱力学プログラム

---

## The Six Papers

| # | Topic | DOI |
|---|-------|-----|
| 1 | Overview and Hypotheses | [10.5281/zenodo.19837588](https://doi.org/10.5281/zenodo.19837588) |
| 2 | Mathematical Foundations of Central Projection | [10.5281/zenodo.19837590](https://doi.org/10.5281/zenodo.19837590) |
| 3 | Unit-Cube Packing in a 4D Ball | [10.5281/zenodo.19837592](https://doi.org/10.5281/zenodo.19837592) |
| 4 | Entropy Correspondence | [10.5281/zenodo.19837594](https://doi.org/10.5281/zenodo.19837594) |
| 5 | Dynamical Connection and Tangherlini Evaporation | [10.5281/zenodo.19837596](https://doi.org/10.5281/zenodo.19837596) |
| 6 | Synthesis and the 4+1 → 3+1 Reduction Problem | [10.5281/zenodo.19837598](https://doi.org/10.5281/zenodo.19837598) |

---

## Starting Point: Three "Field-Shared" Premises

The programme begins from three premises that are **widely shared** in the contemporary literature on black hole thermodynamics and quantum gravity:

- **Premise A**: Black hole entropy is related to a geometric quantity at the horizon (Bekenstein 1973).
- **Premise B**: Black hole evaporation is local to the horizon region (Hawking 1975).
- **Premise C**: Spacetime may be discrete at the Planck scale (the common premise of LQG, CDT, causal set theory).

From these we extract two **conditional propositions**:

> **Hypothesis 1 (Projection Hypothesis)**: A black hole, as it appears in an effective four-dimensional spacetime, may be the projected image of a higher-dimensional geometric structure (specifically four-plus-one-dimensional).

> **Hypothesis 2 (Discrete Drift Hypothesis)**: If spacetime is composed of discrete minimal units, then Hawking evaporation may correspond to the drift of these units across the boundary of a higher-dimensional region.

We make no claim about the physical truth of either hypothesis. The programme is a **conditional investigation**: what follows if both hypotheses are true?

---

## The Mathematical Heart: Unit-Cube Packing in a 4D Ball

Paper 3 contains the core mathematics.

### Problem Setting

Consider a four-dimensional ball
$$B(R) = \{x \in \mathbb{R}^4 : \|x\| \le R\}$$
of radius $R = 2k+1$. How many unit cubes, centred at integer lattice points $c \in \mathbb{Z}^4$, fit inside $B(R)$ without protrusion?

### The Packing Condition

The condition that all vertices of the cube lie in $B(R)$ reduces to:

$$\sum_{i=1}^{4}\left(|c_i|+\tfrac{1}{2}\right)^2 \le (2k+1)^2.$$

We denote the count by $N(k)$. Direct computation gives:

| $k$ | Radius $R$ | $N(k)$ |
|---:|---:|---:|
| 0 | 1 | 1 |
| 1 | 3 | 137 |
| 2 | 5 | 1,545 |
| 5 | 11 | 53,161 |
| 10 | 21 | 818,145 |
| 60 | 121 | 1,028,515,513 |

### Determination of the Asymptotic Constant

The volume deficit $\Delta(R) = V_4(R) - N(k)$ (with $V_4(R) = (\pi^2/2) R^4$) admits, by inclusion–exclusion on the auxiliary body, the asymptotic expansion:

$$\Delta(R) = \frac{16\pi}{3}R^3 - 6\pi R^2 + O(R).$$

The normalized leading constant is:

$$c = \frac{8}{3\pi} \approx 0.84883.$$

Numerical verification: a polynomial fit to data through $k = 30$ confirms this value to within **0.024%**.

### The Lagrange–Jacobi Connection

Under the change of variables $p_i = 2|x_i|+1$, the packing condition becomes a constraint on a sum of four positive odd squares:

$$p_1^2 + p_2^2 + p_3^2 + p_4^2 \le (4k+2)^2.$$

This places the problem squarely within the Lagrange–Jacobi four-square framework, with Jacobi's representation count $r_4(N) = 8\sigma(N)$ (for $N$ not divisible by 4) appearing as the natural generating object. The **special status of four dimensions** thus emerges as a number-theoretic feature.

---

## Connection with Physics

### Entropy Correspondence (Paper 4)

The Bekenstein–Hawking entropy of a four-plus-one-dimensional Schwarzschild–Tangherlini black hole (Tangherlini 1963's higher-dimensional generalization of the Schwarzschild solution) is:

$$S_{BH} = \frac{2\pi^2 r_h^3}{4 G_5} = \frac{\pi^2}{2}r_h^3 \quad (\ell_P = 1).$$

Identifying the horizon radius $r_h$ with the curvature radius $R$ of the subjective space:

$$\frac{\Delta(R)}{S_{BH}} \to \frac{32}{3\pi} \approx 3.397.$$

This is a **constant ratio independent of any free parameter**. Both quantities scale as $R^3$, automatically and from geometry.

### Dynamical Connection (Paper 5)

A Friedmann-like argument gives, for a black hole of mass $M$, $E_{\mathrm{total}} \propto R^2$, and therefore:

- $R \propto M^{1/2}$
- $T_H \propto M^{-1/2}$
- Lifetime $\tau \propto M^{3/2}$
- Evaporation rate $dM/dt \propto -M^{-1/2}$

**Exact agreement with four-plus-one-dimensional Tangherlini.**

### A Significant Finding

A noteworthy clarification emerged during the writing. The programme's original outline flagged the discrepancy "$n = 3/2$ vs $n = -1$" as an open problem. A careful analysis in Paper 5 §5.1 reveals **this is a confusion between 3+1 and 4+1 conventions**:

- **3+1 Schwarzschild**: $dM/dt \propto -M^{-2}$, lifetime $\propto M^3$
- **4+1 Tangherlini**: $dM/dt \propto -M^{-1/2}$, lifetime $\propto M^{3/2}$
- **This programme (4+1)**: $dM/dt \propto -M^{-1/2}$, lifetime $\propto M^{3/2}$

→ **In 4+1 dimensions, exact agreement holds.** The originally flagged "discrepancy" does not exist; instead, the **genuine open problem shifts to "the 4+1 → 3+1 reduction".**

---

## The Remaining Open Problem

Paper 6 formulates the open problem precisely:

> **Open Problem**: How does the 4+1-dimensional black hole thermodynamics established by the programme (lifetime $\propto M^{3/2}$) reduce to the observed 3+1-dimensional thermodynamics (lifetime $\propto M^3$)?

Three candidate mechanisms:

| Candidate | Concept | Evaluation |
|---|---|---|
| Kaluza–Klein compactification | Compactify one of the four spatial axes | Insufficient alone |
| Randall–Sundrum / ADD brane localization | Matter on brane, gravity in bulk | Viable, requires quantification |
| dS-CFT-style holographic projection | 4+1 bulk to 3+1 boundary | Speculative |

The detailed development of these mechanisms is left for future work.

---

## What We Claim and What We Do Not

**Claims**:

- From minimal assumptions (central projection and discreteness), the 4+1-dimensional Tangherlini results are quantitatively reproduced.
- The asymptotic constant $c = 8/(3\pi)$ is derived in closed form by inclusion–exclusion and verified numerically.
- A number-theoretic connection with the Lagrange–Jacobi four-square framework holds.

**Non-claims**:

- We do not claim that physical spacetime is four-plus-one-dimensional.
- We do not claim that higher-dimensional spacetime "physically exists."
- We do not claim a first-principles derivation of the Bekenstein–Hawking entropy.
- We do not claim to displace holographic, loop, or string-theoretic programmes.
- We do not predict any new observable signatures.

All cosmological language is used as **a vocabulary for describing mathematical structures intuitively**, not as physical hypotheses, claims, or predictions about the actual universe.

---

## Author

**Author**: Noriaki Kihara
**Affiliation**: WF System Co., Ltd. / Osaka University, Faculty of Engineering Science (B.Eng., 1983)
**ORCID**: [0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020)
**License**: CC BY 4.0

---

## Related Programmes

This programme is presented as an independent path of inquiry, but it sits alongside the following parallel programmes by the same author:

- **Foundation framework (10 papers)**: Geometric formulation of four-dimensional spacetime via central projection.
- **Phase Equation series (W1–W11)**: Combinatorial derivation of the 62 Standard Model states and the unification of four forces from hypercube structures.
- **Note overview (all 36 papers)**: [https://note.com/kiharanoriaki/n/n63d3c20e6b20](https://note.com/kiharanoriaki/n/n63d3c20e6b20)

---

**Zenn article (more technical exposition)**: [Black Hole Thermodynamics as Projection from a Higher Spacetime Structure](https://zenn.dev/noriaki_kihara/articles/bh-thermodynamics-projection)

---

#BlackHoleThermodynamics #CentralProjection #GnomonicProjection #Tangherlini #Holography #BekensteinHawkingEntropy #HawkingRadiation #HigherDimensionalGravity #NumberTheory #LagrangeFourSquareTheorem #JacobiFourSquareTheorem #deSitter #FourDimensionalBall #LatticePointProblem #InclusionExclusion #MathematicalPhysics #TheoreticalPhysics #QuantumGravity #Preprint #Zenodo
