# Geometric Classification of Standard Model Particles Derived from the Orientation Structure of the 6-Dimensional Hypercube

## — Discrete Geometric Formulation of Color Charge, Mass, and Coupling Strength —

**Author:** Noriaki Kihara (WF System Co., Ltd.)

**Date:** April 2026

**DOI:** [10.5281/zenodo.19646652](https://doi.org/10.5281/zenodo.19646652)

---

## 0. Purpose

In the author's preceding work [7], six spin quantum numbers $s = 0, \frac{1}{2}, 1, \frac{3}{2}, 2, \frac{5}{2}$ were derived as bivectors from the 3-cell orientation structure of the 5-dimensional hypercube. However, the 5-dimensional framework had the following limitations:

1. **Absence of color charge**: No degree of freedom to distinguish the three colors of quarks.
2. **No quark–lepton distinction**: The same orientation corresponds to both quarks and leptons.
3. **Insufficient gauge bosons**: Only 4 spin-1 boson orientations, falling short of the 12 Standard Model states ($\gamma, g_1$–$g_8, W^\pm, Z^0$).
4. **Graviton placement**: The spin-2 boson requires $n = 2$, but with only two non-spatial axes, only one combination was available.

This paper resolves all of these problems by adding a sixth axis to the 5-dimensional framework. The additional axis is a **discrete label axis** qualitatively distinct from space, time, and curvature, which geometrically introduces internal degrees of freedom corresponding to color charge.

All results in this paper are stated as propositions in discrete geometry and number theory.

---

## 1. Basic Structure of the 6-Dimensional Framework

### 1.1 Definition

The 6-dimensional structure treated in this paper consists of the 5-dimensional hypercube (vertices $(\pm 1, \pm 1, \pm 1, \pm 1, \pm 1)$) augmented by a 4-valued discrete label axis $Q$.

The five geometric axes $(x_1, x_2, x_3, x_4, x_5)$ take coordinate values $\pm 1$ and form a hypercube structure. The sixth axis $Q$ is not a geometric coordinate but a discrete label value $Q \in \{0, r, g, b\}$ attached to each state.

### 1.2 Element Counts of the 5-Dimensional Hypercube (Restated)

| $k$ | Element | Count |
|:---:|:-------:|:-----:|
| 0 | Vertices | 32 |
| 1 | Edges | 80 |
| 2 | Square faces | 80 |
| 3 | Cubic cells (3-cells) | 40 |
| 4 | Hypercubic facets (4-facets) | 10 |

This is the same geometric foundation as in the preceding work [7]. The $Q$ axis does not alter the geometric structure of the hypercube but provides an additional quantum number to each state.

### 1.3 Six-Axis Notation

In this paper, each state is denoted by a six-component symbol $(x, y, z, t, R, Q)$. The first five components represent the position in the hypercube (fixed axes as $\pm 1$, free axes as $0$), and the sixth component represents the color label.

---

## 2. Axis Classification: Four Qualitatively Distinct Degrees of Freedom

### 2.1 Extended Assumption

Assumption A (physical classification of five axes) from the preceding works [6, 7] is extended to six axes.

**Assumption B (Physical classification of six axes).** The six axes are classified as follows:

- **Spatial axes** ($x_1, x_2, x_3$): Three axes corresponding to 3-dimensional space directly perceivable by the observer. Sign reversal is not permitted. The three axes are mutually symmetric under the permutation group $S_3$.
- **Temporal axis** ($x_4 = t$): The axis corresponding to time. Sign reversal $\pm 1$ is possible. The range is finite (2 values).
- **Curvature axis** ($x_5 = R$): The axis corresponding to the curvature radius of space. Takes discrete values $R = \pm n^2$ ($n \in \mathbb{N}$). The range is **unbounded**.
- **Color axis** ($x_6 = Q$): The axis corresponding to internal degrees of freedom. Takes discrete label values $Q \in \{0, r, g, b\}$ (4 values).

### 2.2 Four Qualitative Differences

The six axes are classified into four qualitatively distinct types.

| Type | Axis | Values | Meaning | Bound |
|:----:|:----:|:------:|:-------:|:-----:|
| Position | $x, y, z$ | Continuous | "Where is it?" | $S_3$ constrained |
| Direction | $t$ | $\pm 1$ | "Which way?" | 2 values |
| Scale | $R$ | $\pm n^2$ | "How large?" | **Unbounded** |
| Label | $Q$ | $0, r, g, b$ | "What kind?" | 4 values |

**Proposition 2.1.** Among the six axes, only the $R$ axis has an unbounded range. The **type** of a particle (spin, color charge, generation) is determined by the finite discrete structure of $(x, y, z, t, Q)$, while the **mass** of a particle is determined by the sole unbounded degree of freedom ($R$).

### 2.3 Consequences of Q Being a Label

The $Q$ axis represents a label (which kind), not a quantity (how much). The following consequences follow immediately from this property:

1. **Charge quantization**: Since $Q$ takes only discrete values, charge is intrinsically quantized — not as a consequence of continuous U(1) symmetry, but as an inherent property of the axis discreteness.
2. **Physical determination of label count**: The number of values of $Q$ is not constrained by hypercube geometry. The number of values is determined by what the physics requires — namely, the number of color degrees of freedom.
3. **Coupling strength and channel count**: The coupling strength of a force is proportional to the number of interaction channels provided by $Q$ labels (detailed in Section 7).

---

## 3. Stability of the Spin Classification

### 3.1 Spin Derivation (Restated)

As in the preceding work [7], spin is determined by the number $n$ of non-spatial axes contained in the normal bivector.

Since a bivector is composed of two axes, the number of non-spatial axes is limited to $n = 0, 1, 2$. This holds regardless of whether there are two or three non-spatial axes.

| $n$ | Structure of $C$ | Boson spin | Fermion spin |
|:---:|:----------------:|:----------:|:------------:|
| 0 | No non-spatial reversal | 0 | 1/2 |
| 1 | One non-spatial reversal | 1 | 3/2 |
| 2 | Two non-spatial reversals | 2 | 5/2 |

**Proposition 3.1 (Spin stability).** The set of spin quantum numbers $s \in \{0, \frac{1}{2}, 1, \frac{3}{2}, 2, \frac{5}{2}\}$ remains unchanged when the number of non-spatial axes increases from two (5-dimensional) to three (6-dimensional).

*Proof.* A bivector is composed of two axes. Since $n$ is the number of non-spatial axes contained in the bivector, $n \leq 2$. The maximum value of $n$ cannot exceed 2 regardless of how many non-spatial axes exist. $\square$

### 3.2 Increase in the Number of Types

While the set of spin values is invariant, the number of types belonging to each spin increases.

| Spin | $n$ | 5D types | 6D types | Source of increase |
|:----:|:---:|:--------:|:--------:|:------------------:|
| $s = 0, \frac{1}{2}$ | 0 | 3 | 3 | None |
| $s = 1, \frac{3}{2}$ | 1 | 6 | **9** | 3 types $e_i \wedge e_Q$ |
| $s = 2, \frac{5}{2}$ | 2 | 1 | **3** | 2 types $e_t \wedge e_Q$, $e_R \wedge e_Q$ |
| **Total** | | **10** | **15** | $+5$ |

---

## 4. Geometric Derivation of Color Charge

### 4.1 Four Values of the Q Axis and Color

The four values $Q \in \{0, r, g, b\}$ have the following physical correspondence:

| $Q$ value | Color charge | Particle type |
|:---------:|:------------:|:-------------:|
| $0$ | Colorless (color-neutral) | Lepton |
| $r$ | Red | Quark |
| $g$ | Green | Quark |
| $b$ | Blue | Quark |

**Proposition 4.1 (Quark–lepton theorem).** States with $Q = 0$ correspond to leptons (no color charge), while states with $Q \in \{r, g, b\}$ correspond to quarks (with color charge). For a given geometric orientation $(x, y, z, t, R)$, the four values of $Q$ cause each generation of fermions to split into one lepton type and three quark colors.

### 4.2 Derivation of Eight Gluon States

Gluons are spin-1 gauge bosons that mediate color charge transitions. For the three nonzero values $\{r, g, b\}$ of the $Q$ axis, gluons are represented as color–anticolor pairs.

**Proposition 4.2 (Gluon theorem).** The number of independent gauge boson states mediating transitions among color labels $\{r, g, b\}$ is

$$3 \times 3 - 1 = 8$$

The $-1$ arises because the color singlet (the state proportional to $r\bar{r} + g\bar{g} + b\bar{b}$) is not included among the physical gluons.

*Proof.* Color transitions are described by pairs $Q_\text{initial} \to Q_\text{final}$. There are $3^2 = 9$ transitions from three colors to three colors. However, one linear combination of the diagonal components (the color singlet) acts identically on all colors and does not mediate color transitions. Therefore, the number of independent states is $9 - 1 = 8$. $\square$

**Remark.** In the Standard Model, the eight gluon states are derived as the 8-dimensional adjoint representation of the Lie algebra $\mathfrak{su}(3)$ of SU(3). In this paper, the same number 8 is derived purely combinatorially from the 4-label structure of the $Q$ axis. No continuous group structure of SU(3) needs to be assumed.

### 4.3 Complete Classification of Gauge Bosons

According to the type of non-spatial axis, spin-1 gauge bosons are classified as follows:

| Fixed axis | $n=1$ type | Boson orientation | SM correspondence | Count |
|:----------:|:----------:|:-----------------:|:-----------------:|:-----:|
| $t$ axis | Temporal | $(0,0,0,\pm,0,Q)$ | $W^+, W^-$ | 2 |
| $R$ axis | Curvature | $(0,0,0,0,\pm,Q)$ | $\gamma, Z^0$ | 2 |
| $Q$ axis | Color | $(0,0,0,0,0,Q_i \to Q_j)$ | Gluons $g_1$–$g_8$ | 8 |
| **Total** | | | | **12** |

The two $W$ boson states ($t = +1, -1$) correspond to positive and negative charges. $\gamma$ and $Z^0$ correspond to the sign of the $R$ axis, with the distinction between massless ($M(\gamma) = 0$) and massive ($M(Z) \neq 0$). The eight gluon states follow from Proposition 4.2.

The total of 12 agrees with the number of gauge boson states in the Standard Model ($\gamma + g \times 8 + W^+ + W^- + Z^0 = 12$).

---

## 5. Geometric Placement of the Graviton and Spin-2 Bosons

### 5.1 Resolution by Three Non-Spatial Axes

In five dimensions (two non-spatial axes), the only $n = 2$ bivector was $e_4 \wedge e_5$ ($t \wedge R$) — **one combination**. In six dimensions (three non-spatial axes: $t, R, Q$), $n = 2$ increases to **three combinations**:

| Normal bivector | Fixed axes | $C$ matrix |
|:--------------:|:----------:|:----------:|
| $e_t \wedge e_R$ | $t, R$ | $\mathrm{diag}(+,+,+,-,-,Q)$ |
| $e_t \wedge e_Q$ | $t, Q$ | $\mathrm{diag}(+,+,+,-,R,-)$ |
| $e_R \wedge e_Q$ | $R, Q$ | $\mathrm{diag}(+,+,+,t,-,-)$ |

**Proposition 5.1.** In the 6-dimensional framework, spin-2 bosons come in three types.

| Type | Fixed axes | Force property | SM correspondence |
|:----:|:----------:|:--------------:|:-----------------:|
| $tR$ type | $t, R$ | Attractive only ($n=2$, even) | Graviton |
| $tQ$ type | $t, Q$ | Attractive only | **Prediction** |
| $RQ$ type | $R, Q$ | Attractive only | **Prediction** |

The graviton naturally corresponds to the $tR$ type. The $tQ$ and $RQ$ types are predictions of spin-2 bosons with no counterpart in the Standard Model.

### 5.2 Resolution of the 5-Dimensional Problem

The problem in five dimensions — that "the graviton does not fit into a 4-facet center" — was due to dimensional insufficiency. In six dimensions, with three non-spatial axes, multiple $n = 2$ bivectors exist, providing the graviton with a natural geometric position.

---

## 6. Geometric Derivation of Mass via the Fifth Axis (R Axis)

### 6.1 Uniqueness of the R Axis: The Sole Unbounded Degree of Freedom

$R = \pm n^2$ ($n \in \mathbb{N}$) takes unbounded discrete values. The $R$ axis is the only one among the six axes with an unbounded range.

| Axis | Number of possible values |
|:----:|:------------------------:|
| $x, y, z$ | $S_3$ constrained (3 generations) |
| $t$ | 2 |
| $Q$ | 4 |
| $R$ | $\infty$ |

### 6.2 Mass and Area

Due to the structure $R = n^2$, the $R$ axis has the dimension of $L^2$ (area). Among the 80 square faces of the 5-dimensional hypercube, the area of faces containing the $R$ axis is proportional to $n^2$.

**Definition 6.1.** The **signed area** $M(\sigma)$ of an orientation structure $\sigma$ is defined as the sum of the products of orientation signs and areas over all square faces containing the $R$ axis:

$$M(\sigma) = \sum_{\text{square face } f \text{ containing } R} \mathrm{sgn}(f) \cdot \mathrm{Area}(f)$$

### 6.3 Classification of R-Containing Faces

Among the 80 square faces, those containing the $R$ ($= x_5$) axis are classified as follows:

| Face type | Bivector | Type count | Face count | Area scale |
|:---------:|:--------:|:----------:|:----------:|:----------:|
| Spatial–$R$ | $e_i \wedge e_5$ | 3 | $3 \times 8 = 24$ | $L \cdot n^2$ |
| $t$–$R$ | $e_4 \wedge e_5$ | 1 | $8$ | $T \cdot n^2$ |
| **Total** | | **4** | **32** | |

Faces containing the $R$ axis account for 32 out of 80 total faces, or $\frac{2}{5}$ of the total.

### 6.4 Masslessness Condition

**Proposition 6.1.** When the signed area $M(\sigma) = 0$ (complete cancellation of signs among all $R$-containing faces), the orientation corresponds to a massless particle.

The masslessness of the photon and gluons is explained as a consequence of the vanishing signed area of the corresponding orientations.

### 6.5 Mass and Frequency

$R$ is the curvature radius of space and relates to spatial extent (wavelength $\lambda$). Wavelength and frequency are connected by $\lambda \nu = c$, and from the mass–energy equivalence $E = mc^2$ and the quantum relation $E = h\nu$:

$$m = \frac{h\nu}{c^2}$$

**Proposition 6.2.** For a given orientation structure $(x, y, z, t, Q)$, different values of $R = n^2$ yield different mass states. This is analogous to the fundamental tone and overtones of a vibrating string: different values of $n$ form a mass spectrum of particles with identical quantum numbers (spin, color charge, generation).

**Proposition 6.3.** The experimentally measured mass of a particle is the representative value of the most stable vibrational mode (fundamental frequency) for each orientation structure. The "resonances" observed in accelerator experiments may correspond to higher-order vibrational modes (overtones) of the same orientation structure.

---

## 7. Determination of Coupling Strength by Channel Count

### 7.1 Channel Count and Coupling Constants

**Proposition 7.1.** The coupling strength of a force is proportional to the number of $Q$ channels available for the interaction.

| Force | Mediator | $Q$ channels | Relative coupling |
|:-----:|:--------:|:------------:|:-----------------:|
| Strong | $g_1$–$g_8$ | 8 ($3 \times 3 - 1$) | Strongest |
| Electromagnetic | $\gamma$ | Couples to charged ($Q \neq 0$) | Intermediate |
| Weak | $W^\pm, Z^0$ | $t$ axis reversal only | Weak |
| Gravity | $G$ | $Q$-independent ($n = 2$) | Weakest |

*Rationale.* The greater the number of channels, the larger the total transition probability in scattering processes, and hence the larger the effective coupling constant. The strong force utilizes transitions among $Q$ labels (8 channels) and is therefore far stronger than gravity ($n = 2$, spacetime structure only), which is $Q$-independent.

**Remark.** In the Standard Model, the hierarchy of coupling constants $\alpha_s \approx 1, \alpha \approx 1/137, \alpha_W \approx 10^{-6}, \alpha_G \approx 10^{-39}$ cannot be quantitatively explained by channel count differences alone. Channel count provides the qualitative origin of the hierarchy; precise derivation of the values remains a task for future work.

---

## 8. Complete Correspondence Table with Standard Model Particles

### Table A: Confirmed Particles ↔ 6-Dimensional Orientations

#### Bosons

| Particle | Symbol | $s$ | Orientation $(x,y,z,t,R,Q)$ | $n$ | Mass (MeV) | Note |
|----------|--------|-----|-----------------------------|-----|------------|------|
| Higgs | $H^0$ | 0 | $(\pm,0,0,0,0,0)$ etc. | 0 | 125,250 | $Q = 0$, color-neutral |
| W boson | $W^+$ | 1 | $(0,0,0,+,0,0)$ | 1 | 80,379 | $t$ type |
| W boson | $W^-$ | 1 | $(0,0,0,-,0,0)$ | 1 | 80,379 | $t$ type |
| Photon | $\gamma$ | 1 | $(0,0,0,0,+,0)$ | 1 | 0 | $R$ type, $M = 0$ |
| Z boson | $Z^0$ | 1 | $(0,0,0,0,-,0)$ | 1 | 91,188 | $R$ type, $M \neq 0$ |
| Gluon | $g_1$–$g_8$ | 1 | $(0,0,0,0,0,Q_i{\to}Q_j)$ | 1 | 0 | $Q$ type, 8 color transitions |
| Graviton | $G$ | 2 | $(0,0,0,\pm,\pm,0)$ | 2 | 0 | $tR$ type |

#### Fermions: Leptons ($Q = 0$, color-neutral)

| Particle | Symbol | $s$ | Gen. | Orientation $(x,y,z,t,R,Q)$ | Mass (MeV) |
|----------|--------|-----|------|------------------------------|------------|
| Electron | $e^-$ | 1/2 | 1 | $(+,+,0,0,0,0)$ | 0.511 |
| Electron $\nu$ | $\nu_e$ | 1/2 | 1 | $(+,-,0,0,0,0)$ | $< 10^{-6}$ |
| Muon | $\mu^-$ | 1/2 | 2 | $(+,0,+,0,0,0)$ | 105.66 |
| Muon $\nu$ | $\nu_\mu$ | 1/2 | 2 | $(+,0,-,0,0,0)$ | $< 0.17$ |
| Tau | $\tau^-$ | 1/2 | 3 | $(0,+,+,0,0,0)$ | 1,776.9 |
| Tau $\nu$ | $\nu_\tau$ | 1/2 | 3 | $(0,+,-,0,0,0)$ | $< 15.5$ |

Antileptons correspond to orientations with reversed spatial axis signs.

#### Fermions: Quarks ($Q \in \{r, g, b\}$, color-charged)

| Particle | Symbol | $s$ | Gen. | Orientation $(x,y,z,t,R,Q)$ | Mass (MeV) |
|----------|--------|-----|------|------------------------------|------------|
| Up | $u$ | 1/2 | 1 | $(+,+,0,0,0,r/g/b)$ | 2.16 |
| Down | $d$ | 1/2 | 1 | $(+,-,0,0,0,r/g/b)$ | 4.67 |
| Charm | $c$ | 1/2 | 2 | $(+,0,+,0,0,r/g/b)$ | 1,270 |
| Strange | $s_q$ | 1/2 | 2 | $(+,0,-,0,0,r/g/b)$ | 93.4 |
| Top | $t_q$ | 1/2 | 3 | $(0,+,+,0,0,r/g/b)$ | 172,760 |
| Bottom | $b$ | 1/2 | 3 | $(0,+,-,0,0,r/g/b)$ | 4,180 |

Each quark has three color states $Q = r, g, b$. Antiquarks correspond to reversed spatial signs and anticolors.

#### Summary

| Category | $s$ | $Q$ | Types | States | SM states |
|:--------:|:---:|:---:|:-----:|:------:|:---------:|
| Higgs | 0 | 0 | 1 | 1 | 1 |
| Gauge boson ($t$ type) | 1 | 0 | 2 | 2 | 2 ($W^\pm$) |
| Gauge boson ($R$ type) | 1 | 0 | 2 | 2 | 2 ($\gamma, Z$) |
| Gluon ($Q$ type) | 1 | transition | 8 | 8 | 8 |
| Graviton ($tR$ type) | 2 | 0 | 1 | 1 | 1 |
| Lepton | 1/2 | 0 | 6 | 6 | 6 |
| Antilepton | 1/2 | 0 | 6 | 6 | 6 |
| Quark | 1/2 | r,g,b | 6 | 18 | 18 |
| Antiquark | 1/2 | $\bar{r},\bar{g},\bar{b}$ | 6 | 18 | 18 |
| **Total** | | | | **62** | **62** |

### Table B: Predictions of Undiscovered Particles

| Prediction | $s$ | B/F | $n$ | Orientation feature | Gen. | States |
|------------|-----|-----|-----|---------------------|------|--------|
| $tQ$ type boson | 2 | B | 2 | $t, Q$ fixed | — | Predicted |
| $RQ$ type boson | 2 | B | 2 | $R, Q$ fixed | — | Predicted |
| $t$ type fermion | 3/2 | F | 1 | Spatial–$t$ normal | 3 gen. | 12 × 4($Q$) |
| $R$ type fermion | 3/2 | F | 1 | Spatial–$R$ normal | 3 gen. | 12 × 4($Q$) |
| $Q$ type fermion | 3/2 | F | 1 | Spatial–$Q$ normal | 3 gen. | 12 × 4($Q$) |
| $tR$ type fermion | 5/2 | F | 2 | $t$–$R$ normal | None | 4 × 4($Q$) |
| $tQ$ type fermion | 5/2 | F | 2 | $t$–$Q$ normal | None | 4 × 4($Q$) |
| $RQ$ type fermion | 5/2 | F | 2 | $R$–$Q$ normal | None | 4 × 4($Q$) |

**Remark.** The $s = 3/2$ fermions split into three families ($t$ type, $R$ type, $Q$ type), each possessing a three-generation structure. Compared with the gravitino (one species) of supersymmetry theory, this predicts a far richer structure.

### Summary of Improvements from Five Dimensions

| Problem | 5D (W7) | 6D (this paper) | Status |
|:-------:|:-------:|:----------------:|:------:|
| Color charge | Absent | $Q \in \{0,r,g,b\}$ | **Resolved** |
| Quark–lepton distinction | Impossible | $Q = 0$ vs $Q \neq 0$ | **Resolved** |
| Gluon state count | Insufficient (4 vs 12) | $3 \times 3 - 1 = 8$ | **Resolved** |
| Gauge boson total | 4 | 12 ($W^{\pm} + \gamma + Z + g \times 8$) | **Resolved** |
| Graviton | Position unclear | $tR$ type, $n = 2$ | **Resolved** |
| Spin types | 6 | 6 (stable) | **Confirmed** |

---

## 9. Open Problems and Future Work

### 9.1 Quantitative Computation of Signed Area

Explicitly computing the signed area $M(\sigma)$ of Definition 6.1 for each orientation structure and comparing with known particle mass ratios is the most direct test of this theory. Particularly important are:

- $m_\mu / m_e \approx 206.8$: inter-generation mass ratio
- $m_W / m_Z \approx 0.881$: relation to the Weinberg angle
- $m_t / m_e \approx 3.4 \times 10^5$: maximum mass ratio

### 9.2 Derivation of the Gauge Group

There is a suggestive correspondence between the 4-label structure of the $Q$ axis and SU(3), the 2 values of the $t$ axis and SU(2), and the continuous scale of the $R$ axis and U(1). A rigorous proof of the mechanism by which continuous gauge groups emerge from discrete structures remains open.

### 9.3 CKM and PMNS Matrices

How inter-generation mixing (CKM matrix, PMNS matrix) is derived as a breaking of $S_3$ symmetry remains unresolved.

### 9.4 Higgs Mechanism

The Higgs boson corresponds to $n = 0$ (spatial axes only), but how spontaneous symmetry breaking by the Higgs field is understood within this framework is a task for future work.

### 9.5 Additional Spin-2 Boson Predictions

If $tQ$-type and $RQ$-type spin-2 bosons (other than the graviton) exist, their experimental consequences must be clarified.

---

## 10. Conclusion

From the 6-dimensional framework obtained by adding a sixth axis ($Q$ axis: 4-valued discrete label) to the 5-dimensional hypercube, the following have been derived:

1. **Stability of spin classification**: The set of spin types is the same six as in five dimensions ($s = 0, \frac{1}{2}, 1, \frac{3}{2}, 2, \frac{5}{2}$) and is stable under dimensional increase.

2. **Geometric derivation of color charge**: The four labels $\{0, r, g, b\}$ of the $Q$ axis provide color charge and distinguish $Q = 0$ (leptons) from $Q \neq 0$ (quarks).

3. **Derivation of eight gluon states**: The number is determined as $3 \times 3 - 1 = 8$ color transitions without assuming SU(3).

4. **Complete correspondence of 12 gauge boson states**: The three families — $t$ type ($W^\pm$), $R$ type ($\gamma, Z^0$), and $Q$ type ($g_1$–$g_8$) — reproduce all 12 Standard Model gauge boson states exactly.

5. **Natural placement of the graviton**: With three non-spatial axes, the $n = 2$ (spin-2 boson) naturally corresponds to the graviton as the $tR$ type. Additionally, $tQ$-type and $RQ$-type spin-2 bosons are newly predicted.

6. **Areal origin of mass**: The signed area of faces containing the $R$ axis determines mass, and the discreteness $R = n^2$ provides a frequency structure for mass.

7. **Coupling strength**: The number of $Q$ channels qualitatively explains the hierarchy of coupling strengths. Coupling constants originate from the state count of the orientation structure, not as free parameters.

8. **Reproduction of all 62 Standard Model states**: A complete correspondence has been given for 14 boson states ($H + W^\pm + \gamma + Z + g \times 8 + G$) and 48 fermion states (12 leptons + 36 quarks).

The structure in which six coordinate axes divide into four qualitatively distinct types — **position** (3 spatial axes), **direction** (1 temporal axis), **scale** (1 curvature axis), and **label** (1 color axis) — suggests that many of the free parameters of the Standard Model are geometric consequences.

---

## References

[1] Finkelstein, D. & Rubinstein, J. "Connection between Spin, Statistics, and Kinks." *J. Math. Phys.* 9, 1762 (1968).

[2] Atiyah, M.F., Manton, N.S. & Schroers, B.J. "Geometric Models of Matter." *Proc. R. Soc. A* 468, 1252 (2012). arXiv: 1108.5151.

[3] Furey, C. "Standard Model Physics from an Algebra?" arXiv: 1611.09182 (2016).

[4] Kaluza, Th. "Zum Unitätsproblem der Physik." *Sitzungsber. Preuss. Akad. Wiss. Berlin* 966-972 (1921).

[5] Coxeter, H.S.M. *Regular Polytopes.* Dover (1973).

[6] Kihara, N. "Geometric Classification of Spin Derived from the Orientation Structure of the 5-Dimensional Orthoplex." DOI: [10.5281/zenodo.19630972](https://doi.org/10.5281/zenodo.19630972) (2026).

[7] Kihara, N. "Bivectorial Classification of Spin Derived from the Orientation Structure of the 5-Dimensional Hypercube." DOI: [10.5281/zenodo.19643358](https://doi.org/10.5281/zenodo.19643358) (2026).
