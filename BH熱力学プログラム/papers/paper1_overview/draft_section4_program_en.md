# Paper 1 §4 Verification Program — Draft (English)

---

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

---

**Style notes**:
- §4.1 でテーブル形式により仮説への寄与を明示
- §4.4 で蒸発レート指数の不一致を率直に提示
- §4.5 で「establishedすること」と「establishedしないこと」を明示
