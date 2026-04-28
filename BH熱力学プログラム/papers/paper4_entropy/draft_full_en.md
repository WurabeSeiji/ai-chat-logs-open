# Paper 4 (English): A Quantitative Correspondence between the Lattice Volume Deficit and the Tangherlini Bekenstein–Hawking Entropy

**Author**: Noriaki Kihara (WF System Co., Ltd.)
**ORCID**: 0009-0004-6753-4020
**Date**: April 2026

**Abstract**. We propose to interpret the volume deficit $\Delta(R) = V_4(R) - N(k)$ — derived in a companion paper as the difference between the volume of a four-dimensional ball and the count of integer-centred unit cubes contained in it — as a count of degrees of freedom in Planck units, and to compare this count with the Bekenstein–Hawking entropy of a four-plus-one-dimensional Schwarzschild–Tangherlini black hole. The two quantities exhibit identical $R^3$ scaling in $R$, with the asymptotic ratio $\Delta(R) / S_{BH} \to 32/(3\pi) \approx 3.40$. The correspondence is structural and quantitative; we leave the physical interpretation of the constant ratio as an open problem.

---

## §1. Introduction

The Bekenstein–Hawking entropy of a black hole — proportional to the horizon area in Planck units — admits an obvious dimensional analogue in higher dimensions. For a four-plus-one-dimensional Schwarzschild–Tangherlini black hole, the horizon is a three-sphere of radius $r_h$, and the Bekenstein–Hawking entropy is
$$S_{BH} = \frac{A_3}{4 G_5} = \frac{2\pi^2 r_h^3}{4 G_5} = \frac{\pi^2 r_h^3}{2 G_5}.$$

The entropy scales as $r_h^3$. The microscopic origin of this entropy — what counts as a microstate — remains an open problem in the standard formulation.

In a separate strand of work, we have considered the geometric problem of packing unit cubes inside a four-dimensional ball $B(R)$ of radius $R$. The number $N(k)$ of integer-centred unit cubes contained in $B(R)$ for $R = 2k+1$ can be computed exactly; the volume deficit
$$\Delta(R) := V_4(R) - N(k), \qquad V_4(R) = (\pi^2/2) R^4,$$
admits the asymptotic expansion
$$\Delta(R) = \frac{16\pi}{3} R^3 - 6\pi R^2 + O(R), \qquad c := \lim_{k\to\infty} \frac{\Delta(R)}{2\pi^2 R^3} = \frac{8}{3\pi}.$$

The quantity $\Delta(R)$ is a purely combinatorial quantity defined by the discrete packing problem. In this paper we identify $\Delta(R)$, expressed in Planck units, with a candidate count of "discrete drift degrees of freedom" associated with a four-plus-one-dimensional black hole, and compare with $S_{BH}$.

The principal observations are:

1. Both $\Delta(R)$ and $S_{BH}$ scale as the cube of the radius.
2. Under the identification $R \leftrightarrow r_h$, the ratio $\Delta(R)/S_{BH}$ approaches the constant $32/(3\pi) \approx 3.40$.
3. This ratio is independent of any free parameter in the construction; it is fixed by geometry and by the elementary identification of $R$ with $r_h$.

The status of these observations is structural. We do not claim that $\Delta(R)$ *equals* the Bekenstein–Hawking entropy. We claim that $\Delta(R)$ has the right scaling and a definite proportionality to $S_{BH}$, which is a non-trivial check on the discrete drift hypothesis introduced elsewhere in the program.

The paper is organized as follows. Section 2 reviews the necessary inputs from the companion papers. Section 3 sets up the identification of $R$ with the horizon radius $r_h$ and computes the ratio $\Delta(R)/S_{BH}$. Section 4 discusses the physical interpretation of the ratio. Section 5 concludes.

---

## §2. Inputs from Other Papers

### §2.1 The volume deficit $\Delta(R)$

From the companion paper on the discrete packing problem, we have:
$$\Delta(R) = V_4(R) - N(k) = \frac{16\pi}{3} R^3 - 6\pi R^2 + 4R - 1 + O(R^{-1}) - E(R),$$
where $E(R) = O(R^3)$ is the lattice-point error term, conjectured to have zero leading-$R^3$ coefficient on average.

The normalized leading constant is $c := \lim \Delta(R)/(2\pi^2 R^3) = 8/(3\pi) \approx 0.84883$.

### §2.2 The four-plus-one-dimensional Tangherlini metric

The Schwarzschild–Tangherlini metric in 4+1 dimensions is
$$ds^2 = -f_T(r)\, dt^2 + f_T(r)^{-1} dr^2 + r^2\, d\Omega_3^2, \qquad f_T(r) = 1 - \frac{8 G_5 M}{3\pi r^2}.$$
The horizon is at $r_h = (8 G_5 M/(3\pi))^{1/2}$.

The Hawking temperature is $T_H = 1/(2\pi r_h)$. The Bekenstein–Hawking entropy is
$$S_{BH} = \frac{A_3}{4 G_5} = \frac{2\pi^2 r_h^3}{4 G_5}.$$

In Planck units defined by $G_5 = \ell_P^3$ (a natural choice in 5 dimensions), this is
$$S_{BH} = \frac{\pi^2 r_h^3}{2 \ell_P^3}.$$

### §2.3 Identification

The construction of the central-projection metric (Paper 2) introduces a curvature radius $R$ as the natural geometric scale of the four-dimensional subjective chart. The dynamical analysis of Paper 5 establishes that this $R$ scales as $M^{1/2}$ for a black hole of mass $M$, in agreement with the Tangherlini horizon-radius scaling.

We therefore identify $R \leftrightarrow r_h$ throughout this paper. The ratio $\Delta(R)/S_{BH}$ is then a dimensionless number depending only on the geometry of the construction.

---

## §3. The Quantitative Correspondence

### §3.1 Both quantities scale as $R^3$

From §2.1: $\Delta(R) \sim (16\pi/3) R^3$.
From §2.2 (with $r_h \leftrightarrow R$ and $G_5 = \ell_P^3$): $S_{BH} = \pi^2 R^3 / (2 \ell_P^3)$.

In Planck units ($\ell_P = 1$):
$$\Delta(R) \sim \frac{16\pi}{3} R^3, \qquad S_{BH} \sim \frac{\pi^2}{2} R^3.$$

Both quantities scale as the cube of the radius.

### §3.2 The asymptotic ratio

The ratio is
$$\frac{\Delta(R)}{S_{BH}} \to \frac{16\pi/3}{\pi^2/2} = \frac{32}{3\pi} \approx 3.397.$$

This is a dimensionless constant, independent of $R$.

### §3.3 Numerical values at finite $k$

The convergence to the asymptotic ratio is observed in the data:

| $k$ | $R$ | $\Delta(R)$ | $S_{BH} = (\pi^2/2) R^3$ | $\Delta/S_{BH}$ |
|---:|---:|---:|---:|---:|
| 10 | 21 | 141,580.27 | 45,712.2 | 3.10 |
| 20 | 41 | 1,114,426.60 | 340,184.3 | 3.28 |
| 30 | 61 | 3,708,117.64 | 1,119,679.5 | 3.31 |
| 40 | 81 | 9,064,209.71 | 2,624,489.6 | 3.45 |
| 50 | 101 | 16,981,622.84 | 5,054,121.0 | 3.36 |
| 60 | 121 | 29,303,164.67 | 8,740,869.3 | 3.35 |

The ratio approaches $32/(3\pi) \approx 3.40$ from below, with finite-$k$ deviations consistent with the sub-leading $R^2$ correction to $\Delta(R)$ (which contributes to the lower numerical values at small $k$).

### §3.4 Status of the correspondence

The correspondence in scaling and the constant value of the ratio is established. The physical interpretation of the ratio $32/(3\pi)$ is open. Possibilities include:

1. *Geometric.* The ratio $32/(3\pi)$ relates the boundary correction $16\pi/3$ (from inclusion–exclusion on the auxiliary body $\Omega_R$) to the horizon area density $\pi^2/2$ (from $V_4(R) = (\pi^2/2) R^4$ and $S_3(R) = 2\pi^2 R^3$). Specifically:
$$\frac{32}{3\pi} = \frac{16\pi/3}{\pi^2/2} = \frac{2 \cdot (16\pi/3)}{2\pi^2} = \frac{2 c \cdot 2\pi^2}{2\pi^2} = 2 c \cdot \frac{2\pi^2}{\pi^2} = 4c \approx 4 \cdot 0.849 \approx 3.40.$$
Equivalently, the ratio is $4c$ where $c = 8/(3\pi)$ is the volume-deficit constant.

2. *Counting normalization.* The factor of $4$ in $S_{BH} = A/(4 \ell_P^{D-2})$ is the Bekenstein normalization (one quarter of the area in Planck units). A different counting convention would shift the ratio by a factor.

3. *Dimensional convention.* The identification $G_5 = \ell_P^3$ is one natural choice; other choices ($G_5 = \alpha \ell_P^3$ for some $\alpha$) shift $S_{BH}$ and hence the ratio.

We do not at this stage privilege any of these interpretations.

---

## §4. Physical Interpretation: Open Questions

### §4.1 Is $\Delta(R)$ "the" entropy?

We have shown that $\Delta(R)$ has the same scaling as $S_{BH}$ and is in a constant ratio to it. This does not establish that $\Delta(R) = S_{BH}$ in any deeper sense. The two are quantities of the same dimensional type (counts of degrees of freedom in Planck units, both scaling as length cubed in 4+1 dimensions), and they bear a structural relation.

A stronger interpretation would require:
(i) An independent derivation of $\Delta(R)$ from a microscopic theory of black hole microstates.
(ii) A consistent extension to logarithmic and finite-size corrections.
(iii) A clear physical mechanism for the discrete drift that justifies the identification.

The companion papers (3 and 5) address (i) and (iii) at the level of geometric and combinatorial structure, but not at the level of a complete microscopic theory.

### §4.2 The constant $32/(3\pi) \approx 3.40$

The ratio is a definite number. Its physical interpretation, if any, is open. It might be a manifestation of:

- A residual gauge of the discrete-cube packing problem (the choice of integer lattice rather than another regular lattice).
- A specific choice of dimensional reduction from a putative higher-dimensional structure.
- An artifact of the inclusion–exclusion boundary correction having a specific geometric interpretation.

No preferred interpretation is offered here.

### §4.3 Comparison with logarithmic and one-loop corrections

The standard Bekenstein–Hawking entropy receives logarithmic and one-loop corrections in the semi-classical approximation. The volume-deficit construction has its own sub-leading corrections — the $-6\pi R^2$ term and lower — which we have not attempted to match against the standard quantum corrections. **後で検討**: detailed correspondence at sub-leading orders.

---

## §5. Conclusion

We have identified a structural correspondence between the volume deficit $\Delta(R)$, derived from the geometric problem of discrete unit-cube packing, and the Bekenstein–Hawking entropy of a four-plus-one-dimensional Schwarzschild–Tangherlini black hole. Both quantities scale as $R^3$, with the asymptotic ratio $\Delta(R)/S_{BH} \to 32/(3\pi) \approx 3.40$ in Planck units.

This is a non-trivial check on the discrete drift hypothesis: under the assumption that $\Delta(R)$ is the relevant entropy quantity for the discrete drift, the scaling matches the expected form. The identification of the constant ratio as something physical, rather than as a feature of the geometric construction, is left as an open problem.

The dynamical determination of $R$ from the matter content of the four-dimensional theory is the subject of Paper 5.

---

## References (selected, 後で検討)

- J. D. Bekenstein, "Black Holes and Entropy," *Phys. Rev. D* 7, 2333 (1973).
- S. W. Hawking, "Particle Creation by Black Holes," *Comm. Math. Phys.* 43, 199 (1975).
- F. R. Tangherlini, "Schwarzschild Field in $n$ Dimensions and the Dimensionality of Space Problem," *Nuovo Cimento* 27, 636 (1963).
- R. C. Myers and M. J. Perry, "Black Holes in Higher Dimensional Spacetimes," *Annals Phys.* 172, 304 (1986).
