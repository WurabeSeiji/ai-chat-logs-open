# Geometric Structure of Generation Mixing — Qualitative Derivation of the CKM and PMNS Matrices

## —— A Structural Answer to the Generation Mixing Problem of [4] §9.11.6 ——

**Author**: Noriaki Kihara (WF System Co., Ltd.)

**ORCID**: https://orcid.org/0009-0002-7952-9855

**DOI**: https://doi.org/10.5281/zenodo.19731602

**Date**: April 2026

---

## §0 Purpose of This Paper

In [4] §9.11.6, we declared that "how generation mixing (the CKM matrix and the PMNS matrix) is derived as the breaking of $S_3$ symmetry of the $xyz$ axes remains an open problem."

The present paper responds to this declaration and demonstrates the following.

1. The diagonally dominant structure of the CKM matrix follows directly from the $\mathrm{SO}(3) \to \mathrm{SO}(2)$ symmetry breaking caused by the Higgs (standing wave) axis selection (§1).
2. The large mixing angle structure of the PMNS matrix follows directly from the (spatial, spatial)-type 2-face structure of neutrinos (§2).
3. The structural difference between the CKM and PMNS matrices reduces to the geometric difference between (spatial, $t$)-type and (spatial, spatial)-type (§3).

The claims of this paper are limited to the derivation of the **qualitative structure** of the mixing matrices (diagonal dominance vs. large mixing angles). The derivation of specific numerical values of the mixing angles is an open problem sharing the same root as the quantitative derivation of mass ratios ([4] §9.11.4) and lies outside the scope of this paper.

**References**:
- [4] Set Structure and Combinatorial Properties of the 6-Dimensional Hyperrectangle (W8)
- [5] Interactions of Shape-Invariant Waves (W11)

---

## §1 The CKM Matrix: Generation Mixing of the (Spatial, $t$)-Type

### 1.1 Correspondence Between Generations and Spatial Axes

By [4] §9.1, charged fermions (quarks and charged leptons) have 2-faces of the (spatial 1-axis, $t$)-type. The three spatial axes $x, y, z$ correspond to the three generations:

| Generation | Spatial axis | Up-type quark | Down-type quark |
|:---:|:---:|:---:|:---:|
| 1 | $x$ | $u = (+, 0, 0, +, 0, Q)$ | $d = (+, 0, 0, -, 0, Q)$ |
| 2 | $y$ | $c = (0, +, 0, +, 0, Q)$ | $s = (0, +, 0, -, 0, Q)$ |
| 3 | $z$ | $t = (0, 0, +, +, 0, Q)$ | $b = (0, 0, +, -, 0, Q)$ |

### 1.2 Axis Selection by the Higgs and Symmetry Breaking

By [4] §9.11.7, the Higgs boson (standing wave, $\kappa = 1$) takes a nonzero component along one of the three $xyz$ axes. In the initial state, the $xyz$ axes are completely equivalent under $\mathrm{SO}(3)$ symmetry, but the symmetry is spontaneously broken when the Higgs selects one axis. We denote the selected axis as $z$.

As a result:
- $z$ axis: the nonzero axis of the Higgs → direct coupling → largest mass (third generation)
- $x, y$ axes: do not share the nonzero axis of the Higgs → only indirect coupling → lighter masses (first and second generations)

The symmetry is broken as $\mathrm{SO}(3) \to \mathrm{SO}(2)$, where the $z$ axis is separated while an approximate $\mathrm{SO}(2)$ symmetry remains between the $x$ and $y$ axes.

### 1.3 Diagonally Dominant Structure of the CKM Matrix

**Proposition 1.1**. In generation mixing of (spatial, $t$)-type fermions, the Higgs $z$-axis selection leads to the following hierarchical structure:

(i) Mixing between the third generation and the first/second generations is small (separation of the $z$ axis).

(ii) Mixing between the first and second generations is larger than (i) (approximate symmetry of the $x, y$ axes).

(iii) The CKM matrix is diagonally dominant.

*Argument*. The weak interaction mediated by $W^{\pm}$ is an isospin transition ([5] §6), which flips the $t$-axis sign of the fermion and converts up-type quarks into down-type quarks. Since each generation of quarks occupies only one spatial axis, inter-generation transitions involve changing the spatial axis.

Due to the Higgs $z$-axis selection, fermions on the $z$ axis (third generation) have a mass scale significantly different from those on other axes. The deviation between mass eigenstates and flavor eigenstates (mixing angle) depends on the coupling strength between axes, which is inversely proportional to the mass difference.

(i) Since $m_3 \gg m_{1,2}$, the mixing angles between the third generation and the first/second generations are small.

(ii) Since the $x$ and $y$ axes are approximately equivalent under $\mathrm{SO}(2)$ symmetry, $m_1 \approx m_2$, and the mixing angle between the first and second generations is larger than in (i).

(iii) From (i) and (ii), the diagonal elements of the CKM matrix are close to 1 and the off-diagonal elements are small (diagonal dominance). $\square$

**Comparison with experimental values**:

$$V_{\text{CKM}} \approx \begin{pmatrix} 0.974 & 0.225 & 0.004 \\ 0.225 & 0.973 & 0.041 \\ 0.009 & 0.040 & 0.999 \end{pmatrix}$$

- $|V_{ub}| = 0.004$, $|V_{td}| = 0.009$: mixing between the 3rd and 1st generations is smallest ✓
- $|V_{cb}| = 0.041$, $|V_{ts}| = 0.040$: mixing between the 3rd and 2nd generations is small ✓
- $|V_{us}| = 0.225$: mixing between the 1st and 2nd generations is largest ✓

---

## §2 The PMNS Matrix: Generation Mixing of the (Spatial, Spatial)-Type

### 2.1 The 2-Face Structure of Neutrinos

By [4] §9.9 Table A, neutrinos have 2-faces of the (spatial, spatial)-type. The $\binom{3}{2} = 3$ combinations of spatial axes correspond to the three flavors:

| Flavor | Nonzero axes | Set |
|:---:|:---:|:---|
| $\nu_e$ | $(x, y)$ | $(+, +, 0, 0, 0, 0)$ |
| $\nu_\mu$ | $(y, z)$ | $(0, +, +, 0, 0, 0)$ |
| $\nu_\tau$ | $(x, z)$ | $(+, 0, +, 0, 0, 0)$ |

### 2.2 Axis-Sharing Structure

**Proposition 2.1**. Any two neutrino flavors share exactly one spatial axis.

*Proof*. Among the three combinations of choosing 2 elements from the 3-element set $\{x, y, z\}$ — namely $\{x,y\}, \{y,z\}, \{x,z\}$ — the intersection of any two combinations consists of exactly one element:

| Pair | Shared axis | Non-shared axes |
|:---:|:---:|:---:|
| $\nu_e, \nu_\mu$ | $y$ | $x \leftrightarrow z$ |
| $\nu_\mu, \nu_\tau$ | $z$ | $y \leftrightarrow x$ |
| $\nu_e, \nu_\tau$ | $x$ | $y \leftrightarrow z$ |

This is a combinatorial property of $\binom{3}{2}$ and holds in general for the case of 3 axes with 2 chosen. $\square$

### 2.3 Large Mixing Angle Structure of the PMNS Matrix

**Proposition 2.2**. In generation mixing of (spatial, spatial)-type fermions, the Higgs $z$-axis selection leads to the following structure:

(i) Both $\nu_\mu$ and $\nu_\tau$ contain the $z$ axis (Higgs axis), and their non-shared axes $(y, x)$ are unaffected by the Higgs, so $\nu_\mu \leftrightarrow \nu_\tau$ mixing is nearly maximal.

(ii) Since $\nu_e$ does not contain the $z$ axis, the mixing $\nu_e \leftrightarrow \nu_\mu$ and $\nu_e \leftrightarrow \nu_\tau$ is smaller than (i).

(iii) The PMNS matrix has large off-diagonal elements (in contrast to the CKM matrix).

*Argument*.

(i) $\nu_\mu = (y, z)$ and $\nu_\tau = (x, z)$ share the axis $z$, with non-shared axes $y$ and $x$, respectively. Since $x$ and $y$ are unaffected by the Higgs and the $\mathrm{SO}(2)$ symmetry remains, the coupling strengths of $\nu_\mu$ and $\nu_\tau$ to the Higgs are nearly equal. The mixing angle of two states with equal coupling strengths approaches the maximum ($\theta \approx 45°$).

(ii) $\nu_e = (x, y)$ does not contain the $z$ axis. The coupling of $\nu_e$ to the Higgs is limited to indirect coupling through the $x, y$ axes, and its coupling structure differs from that of $\nu_\mu, \nu_\tau$, which contain the $z$ axis. This asymmetry causes the mixing angle of $\nu_e$ with $\nu_{\mu,\tau}$ to deviate from the maximum.

(iii) The coexistence of the maximal mixing of (i) and the moderate mixing of (ii) results in large off-diagonal elements in the PMNS matrix overall. $\square$

**Comparison with experimental values**:

- $\theta_{23} \approx 49°$ (atmospheric angle, $\nu_\mu \leftrightarrow \nu_\tau$): nearly maximal ✓ — Proposition (i)
- $\theta_{12} \approx 33°$ (solar angle, $\nu_e \leftrightarrow \nu_\mu$): large but not maximal ✓ — Proposition (ii)
- $\theta_{13} \approx 8.6°$ (reactor angle, $\nu_e \leftrightarrow \nu_\tau$): smallest ✓ — Proposition (ii)

---

## §3 Origin of the Structural Difference

### 3.1 (Spatial, $t$)-Type vs. (Spatial, Spatial)-Type

**Proposition 3.1**. The structural difference between the CKM matrix (diagonally dominant) and the PMNS matrix (large mixing angles) reduces to the difference in axis composition of the 2-faces.

| Property | (Spatial, $t$)-type | (Spatial, spatial)-type |
|:---|:---|:---|
| Number of occupied spatial axes | **1** | **2** |
| Generations containing the Higgs axis ($z$) | Third generation only | **Two**: $\nu_\mu$ and $\nu_\tau$ |
| Axis sharing between two flavors | None (each generation has an independent axis) | Always share 1 axis (Proposition 2.1) |
| Separation by the Higgs | Split into 1 vs. 2 ($z$ vs. $x, y$) | Weak separation (2 containing $z$ vs. 1 not containing $z$) |
| Consequence for mixing angles | Small (diagonally dominant) | Large |

*Argument*. In the (spatial, $t$)-type, each generation occupies one independent spatial axis, with no axis sharing between generations. The Higgs $z$-axis selection separates only one generation, and the coupling between the remaining two generations is also indirect, so the mixing angles are small overall.

In the (spatial, spatial)-type, each flavor occupies two spatial axes, and any two flavors share one axis (Proposition 2.1). The Higgs $z$ axis is contained in two out of three flavors, so the separation effect is weak, and the direct coupling through shared axes makes the mixing angles large. $\square$

### 3.2 Relation to Supplementary Lecture 1

As shown in Supplementary Lecture 1 §2.2, (spatial, spatial)-type neutrinos have $\text{sign}(M_F) = +1$ (positive orientation), which differs in orientation from (spatial, $t$)-type charged fermions. By Supplementary Lecture 1 Proposition 2.2, antineutrinos also have $M_F^{\text{anti}} = M_F$ (sign-invariant).

The structural difference in generation mixing shares the same geometric origin as this difference in orientation — the difference in 2-face structure between the (spatial, $t$)-type and the (spatial, spatial)-type.

---

## §4 Conclusions

In this paper, we have provided the following structural answers to the generation mixing problem declared open in [4] §9.11.6.

1. **The diagonally dominant structure of the CKM matrix** follows from the fact that each generation of (spatial, $t$)-type fermions occupies an independent spatial axis, and the Higgs $z$-axis selection induces $\mathrm{SO}(3) \to \mathrm{SO}(2)$ symmetry breaking that separates the third generation (Proposition 1.1).

2. **The large mixing angle structure of the PMNS matrix** follows from the fact that each flavor of (spatial, spatial)-type neutrinos occupies two spatial axes, and any two flavors share one axis (Proposition 2.1). In particular, since both $\nu_\mu$ and $\nu_\tau$ contain the Higgs axis $z$, $\theta_{23}$ is nearly maximal (Proposition 2.2).

3. **The structural difference between the CKM and PMNS matrices** reduces to the difference in 2-face structure between the (spatial, $t$)-type and the (spatial, spatial)-type (Proposition 3.1).

The derivation of specific numerical values of the mixing angles remains as an open problem sharing the same root as the quantitative derivation of mass ratios in [4] §9.11.4.

---

## References

- [4] Noriaki Kihara, "Set Structure and Combinatorial Properties of the 6-Dimensional Hyperrectangle," Zenodo preprint, DOI: 10.5281/zenodo.19721125 (2026).
- [5] Noriaki Kihara, "Interactions of Shape-Invariant Waves — Axial Displacement Transfer, Wave Packet Deformation, and Retroactive Constitution of Causality," Zenodo preprint, DOI: 10.5281/zenodo.19721128 (2026).
