# A Geometric Identity for the Fine-Structure Constant: From the 4D Unit Ball Volume and its Cube-Packing Deficit (Paper 7: α-Identity)

**Author**: Noriaki Kihara (WF System Co., Ltd.)
**ORCID**: [0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020)
**Date**: April 2026
**Concept DOI (always resolves to the latest version)**: [10.5281/zenodo.19869266](https://doi.org/10.5281/zenodo.19869266)
**v1 Version DOI**: [10.5281/zenodo.19869267](https://doi.org/10.5281/zenodo.19869267)
**v2 Version DOI (adds §6.5 duality outlook)**: [10.5281/zenodo.19875074](https://doi.org/10.5281/zenodo.19875074)
**v3 Version DOI (latest, precision corrected to 8.7 ppb)**: [10.5281/zenodo.19876200](https://doi.org/10.5281/zenodo.19876200)

---

## Abstract

Let $N(1)$ denote the count of integer-centred unit cubes fully contained in the 4-dimensional ball $B(R)$ of radius $R = 3$. Direct enumeration gives $N(1) = 137$. The volume of the 4D unit ball is $V_4(1) = \pi^2/2$. The combination
$$ \alpha^{-1} \stackrel{?}{=} N(1) + V_4(1) \cdot \alpha $$
is a self-consistent quadratic equation
$$ \tfrac{\pi^2}{2}\alpha^2 + 137 \alpha - 1 = 0 $$
whose positive root is $\alpha = 7.29735194\ldots \times 10^{-3}$, agreeing with the observed value $\alpha_{\text{obs}} = 7.2973525693 \times 10^{-3}$ (CODATA 2018) to a relative accuracy of $8.7 \times 10^{-8}$ (about **8.7 ppb**, with $\alpha^{-1}$ matching to seven significant digits). In this paper we (i) re-derive $N(1) = 137$ and $V_4(1) = \pi^2/2$ by elementary calculation, (ii) verify the numerical validity of the identity, and (iii) propose a structural reading as "tree-level coupling on the inside (packed cubes) plus self-energy correction on the outside (boundary gap)". Whether the identity reflects geometric necessity or numerical coincidence is left as an open problem.

---

## §1. Introduction

The fine-structure constant $\alpha = e^2/(4\pi\varepsilon_0 \hbar c) \approx 1/137.036$ is among the most precisely measured dimensionless physical constants and governs the strength of the electromagnetic interaction. That its value lies "near 137" has been known since Eddington (1929), and number-theoretic explanations have been attempted for a century. These attempts—from Eddington's naive integer-fitting to Wyler's (1969) multi-element formula—have nearly all been rejected, either by improved measurement precision or by the absence of independent verification.

This paper observes that a problem in 4-dimensional geometry, motivated independently, satisfies a self-consistent equation that reproduces both the integer part 137 and the fractional part 0.036 of $\alpha^{-1}$ using only **parameters that cannot be externally tuned**.

Specifically, from two elementary facts—(i) the classical problem of packing integer-centred unit cubes in a 4D ball (a natural generalisation of the Gauss circle problem to 4D), and (ii) the volume formula for the 4D ball—an algebraic identity emerges that predicts $\alpha$ to **ppb (parts per billion) accuracy**.

This paper does not propose a number-theoretic "single integer fit" but expresses $\alpha$ as the **solution of a self-consistent equation**. Since the correction term $V_4(1) \cdot \alpha$ contains $\alpha$ itself, the identity has the structure of a perturbation expansion
$$ 1 = 137 \alpha + V_4(1) \alpha^2, $$
analogous to the renormalisation structure in QFT. This is the decisive difference from naive number-theoretic predecessors.

**Structure of the paper.** §2 re-derives $N(1) = 137$ and $V_4(1) = \pi^2/2$ independently. §3 performs the numerical verification of the identity. §4 presents a physical reading (inside/outside decomposition). §5 enumerates candidate mechanisms and §6 discusses open problems and the residual. §7 concludes.

The paper is written self-contained: every numerical value is reproducible by elementary calculation in the main text, and references to external papers serve only for context, not for the verification of the main claims.

---

## §2. Re-derivation of the Two Geometric Quantities

### §2.1 $N(1) = 137$: Direct Enumeration

Consider the 4-dimensional ball $B(3) = \{x \in \mathbb{R}^4 : \|x\|_2 \le 3\}$. An integer-centred cube $C(c) = \prod_i [c_i - 1/2, c_i + 1/2]$ with $c \in \mathbb{Z}^4$ is **fully contained** in $B(3)$ iff every vertex of $C(c)$ lies in $B(3)$, equivalently
$$ \sum_{i=1}^{4} (|c_i| + 1/2)^2 \le 9. $$

Setting $s_i := |c_i| \in \mathbb{Z}_{\ge 0}$, the inequality becomes
$$ \sum_i s_i^2 + \sum_i s_i + 1 \le 9 \quad\iff\quad \sum_i s_i^2 + \sum_i s_i \le 8. $$

Enumerating $(s_1, s_2, s_3, s_4)$ in increasing order:

| Pattern $(s_1,s_2,s_3,s_4)$ | $\sum s_i^2 + \sum s_i$ | Permutations | Sign DOF | Cube count |
|---|---|---|---|---|
| $(0,0,0,0)$ | $0$ | $1$ | $1$ | $1$ |
| $(1,0,0,0)$ | $2$ | $4$ | $2^1$ | $8$ |
| $(1,1,0,0)$ | $4$ | $6$ | $2^2$ | $24$ |
| $(1,1,1,0)$ | $6$ | $4$ | $2^3$ | $32$ |
| $(1,1,1,1)$ | $8$ | $1$ | $2^4$ | $16$ |
| $(2,0,0,0)$ | $6$ | $4$ | $2^1$ | $8$ |
| $(2,1,0,0)$ | $8$ | $12$ | $2^2$ | $48$ |
| $(2,1,1,0)$ | $10$ | — | — | excluded |
| $(2,2,0,0)$ | $12$ | — | — | excluded |
| $(3,0,0,0)$ | $12$ | — | — | excluded |

**Permutations** count placement choices for non-zero $s_i$ among 4 slots; **sign DOF** is $\pm$ for each non-zero $s_i$. Total:
$$ N(1) = 1 + 8 + 24 + 32 + 16 + 8 + 48 = \boxed{137}. $$

This calculation is fully elementary and rigorous, requiring no external reference.

### §2.2 $V_4(R) = (\pi^2/2) R^4$ and $V_4(1) = \pi^2/2$

The $n$-dimensional ball volume formula:
$$ V_n(R) = \frac{\pi^{n/2}}{\Gamma(n/2 + 1)} R^n. $$

For $n = 4$, $\Gamma(3) = 2$ gives
$$ V_4(R) = \frac{\pi^2}{2} R^4. $$

In particular $V_4(1) = \pi^2/2 \approx 4.9348022$.

### §2.3 Unit Value of the 4+1D Schwarzschild–Tangherlini BH Entropy

For reference, the horizon area of the 4+1D Tangherlini (1963) solution is the 3-dimensional measure of a 3-sphere of radius $r_h$:
$$ A_3 = 2\pi^2 r_h^3, $$
and the Bekenstein–Hawking entropy in Planck units ($G_5 = \ell_P^3 = 1$):
$$ S_{BH}(r_h) = \frac{A_3}{4} = \frac{\pi^2}{2} r_h^3. $$

At $r_h = 1$:
$$ S_{BH}(1) = \frac{\pi^2}{2} = V_4(1). $$

That is, the unit value of the 4+1D BH entropy coincides with the 4D unit ball volume. This follows directly from standard 4D and 4+1D geometric formulae and is not strictly required for the main result, but is used in the physical reading in §4.

---

## §3. The Identity and its Numerical Verification

### §3.1 Main Result

Combining the two quantities $N(1) = 137$ and $V_4(1) = \pi^2/2$ from §2, consider the algebraic identity:

$$ \boxed{\; \alpha^{-1} = N(1) + V_4(1) \cdot \alpha \;} \tag{1} $$

This is equivalent to a self-consistent quadratic in $\alpha$:
$$ V_4(1) \cdot \alpha^2 + N(1) \cdot \alpha - 1 = 0 \tag{2} $$
$$ \frac{\pi^2}{2}\alpha^2 + 137 \alpha - 1 = 0. $$

### §3.2 Solution and Numerical Verification

Positive root of (2):
$$ \alpha = \frac{-137 + \sqrt{137^2 + 2\pi^2}}{\pi^2}. $$

50-digit high-precision evaluation (mpmath):
$$ \sqrt{18788.7392088022\ldots} = 137.07202197677751\ldots $$
$$ \alpha = \frac{0.07202197677751\ldots}{9.86960440108936\ldots} = 7.29735193535903\ldots \times 10^{-3}. $$

Reciprocal:
$$ \alpha^{-1} = 137.03601098838876\ldots $$

Observed value (CODATA 2018):
$$ \alpha_{\text{obs}} = 7.2973525693(11) \times 10^{-3}, \qquad \alpha_{\text{obs}}^{-1} = 137.035999084(21). $$

| Quantity | Predicted (high precision) | Observed | Relative error |
|---|---|---|---|
| $\alpha$ | $7.29735194 \times 10^{-3}$ | $7.29735257 \times 10^{-3}$ | $\mathbf{8.7 \times 10^{-8}}$ |
| $\alpha^{-1}$ | $137.03601099$ | $137.03599908$ | $8.7 \times 10^{-8}$ |
| $\alpha^{-1} - 137$ | $0.03601099$ | $0.03599908$ | $3.3 \times 10^{-4}$ |

The integer 137 matches exactly. As $\alpha$ itself, the prediction matches observation to **8.7 ppb** (parts per billion, $\approx 8.7 \times 10^{-8}$); $\alpha^{-1}$ matches to four decimal places (137.0360).

**Note**: Previous paper versions v1/v2 reported "relative error $1.4 \times 10^{-4}$ (0.02%)", which was an under-estimation due to coarse 4-digit rounding ($\alpha \approx 7.2984 \times 10^{-3}$). v3 corrects this following high-precision external verification by Grok (xAI).

### §3.3 Equivalent Forms

(1) admits the equivalent forms
$$ 1 = 137 \alpha + \frac{\pi^2}{2} \alpha^2 \tag{3} $$
which is a power expansion in $\alpha$: tree-level $137\alpha$ + first-order correction $V_4(1)\alpha^2$ summing to unity. Physically reminiscent of perturbative expansions in QFT.

Alternative form:
$$ 137 = \alpha^{-1} - \frac{\pi^2}{2} \alpha = \alpha^{-1}\left(1 - \frac{\pi^2}{2}\alpha^2\right). \tag{4} $$
"Integer 137 = observed $\alpha^{-1}$ minus self-energy term = bare coupling."

---

## §4. Physical Reading: Inside/Outside Decomposition

The identity (1) was established as a mathematical fact in §3. This section attempts a physical reading. **Emphasis**: the following interpretation is a structural metaphor, not a first-principles derivation.

### §4.1 Geometry of the R=3 System

As seen in §2.1, the 4D ball of radius $R=3$ contains exactly **137 fully inscribed unit cubes**. The ball volume is $V_4(3) = 81 \cdot \pi^2/2 \approx 399.72$. The deficit:
$$ \Delta(3) = V_4(3) - N(1) = 81 \cdot \frac{\pi^2}{2} - 137 \approx 262.72. $$
This is the volume of the ball not occupied by full cubes—**65.7\%** of the total ball volume. In the R=3 system the boundary dominates over the interior.

| Region | Volume | Structure |
|---|---|---|
| Interior (137 cubes) | 137 | discrete, rigid |
| Boundary (gap) | 262.72 (66\%) | continuous, surplus |
| Total (ball) | 399.72 | $V_4(3)$ |

### §4.2 Naturalness of the Inside/Outside Decomposition

Restating (3):
$$ 1 = 137 \alpha + \frac{\pi^2}{2} \alpha^2. $$

We re-read each term:

| Term | Physical reading | Origin |
|---|---|---|
| $137 \alpha$ | Inside tree-level coupling: 137 geometric contact points contributing at $\alpha$ | $N(1) = 137$ (packing integer) |
| $\frac{\pi^2}{2} \alpha^2$ | Outside self-energy correction: boundary degrees of freedom contributing at $\alpha^2$ | $V_4(1) = \pi^2/2$ (unit ball volume) |
| Total = 1 | Normalisation of total coupling | — |

Physical picture: an electron (or analogous solitonic structure) organises itself at scale $R=3$. Its 137 internal geometric units carry the tree-level electromagnetic coupling, and its boundary layer (which dominates 66\% of the volume) generates a self-correlated renormalisation correction at order $\alpha^2$.

### §4.3 Why $V_4(1)$? — Unit Influence-Region Reading

The correction term contains $V_4(R=1) = \pi^2/2$, not $V_4(R=3) = 81\pi^2/2$. This suggests that the self-energy correction is a **local process at the unit-cube scale**:

- Each cube correlates with adjacent unit influence regions (4D balls of radius 1, volume $V_4(1) = \pi^2/2$)
- All 137 cubes share the same unit influence region (**intensive**)
- The total correction collapses to a single unit influence-region weight

An equivalent reading: since the unit 4+1D BH entropy $S_{BH}(r_h=1) = \pi^2/2$ (§2.3) coincides with $V_4(1)$, the correction term can also be read as "**unit-boundary entropy degrees of freedom multiplied by coupling**". These are mathematically the same value with different physical metaphors.

---

## §5. Candidate Mechanisms

To physically justify (1) requires an independent derivation of the correction $V_4(1) \cdot \alpha$. We list candidates.

### §5.1 Candidate A: Vacuum Polarisation (Standard QED Analogue)

QED vacuum polarisation dresses the photon propagator with virtual fermion loops, shifting the low-energy coupling from the bare value:
$$ \frac{1}{\alpha(0)} - \frac{1}{\alpha(\Lambda)} = \frac{1}{3\pi}\sum_\ell \theta(\Lambda - m_\ell)\left[\ln\frac{\Lambda^2}{m_\ell^2} - \frac{5}{3}\right]. $$

In a lattice model, $\Lambda$ corresponds to a natural cutoff at the lattice spacing. Restricting to the electron loop and demanding the observed shift $\Delta(1/\alpha) = 0.036$, the implied cutoff is $\Lambda \sim 0.6$ MeV—just above $m_e$. The physical interpretation requires further study, but is consistent with kinematic thresholds at the lightest lepton scale.

### §5.2 Candidate B: Boundary Self-Correlation

From the inside/outside decomposition in §4.2, where the boundary occupies 66\%, a direct mechanism is that the self-correlation energy of boundary modes generates the correction:
$$ \delta\alpha^{-1} \sim \int_{\partial} \langle \mathcal{O}(x) \mathcal{O}(y) \rangle \, dx\, dy \stackrel{?}{=} V_4(1) \cdot \alpha. $$

Here $\mathcal{O}$ denotes a boundary field and $\partial$ the boundary layer of the R=3 ball. Concrete computation is left for future work.

### §5.3 Candidate C: Beat Coupling with a Longitudinal Scalar Field

Independent thought-experiment work suggests identifying the Higgs (or analogous scalar longitudinal mode) with the longitudinal vibrations of a minimal-phase lattice. **Beat interference** between transverse-wave solitons (e.g. electrons) and the longitudinal Higgs vibration enables small energy transfer between modes that are nominally orthogonal. This mechanism activates when curvature (background geometry) breaks the orthogonality, yielding a correction to the observed $\alpha$. Details deferred to a separate paper.

### §5.4 Candidate Comparison

| Candidate | Mechanism | Verifiability | Estimated contribution |
|---|---|---|---|
| A Vacuum polarisation | Virtual lepton pair production | Semi-quantitative via existing QED | $\sim 10^{-3}$ to $\sim 10^{-1}$ (cutoff-dependent) |
| B Boundary self-correlation | Direct implementation of inside/outside | Lattice-computable | not yet estimated |
| C Scalar longitudinal beat | Thought-experiment based | Conjectural mechanism | not yet estimated |

These are not exclusive; they may represent different facets of the same phenomenon.

---

## §6. Open Problems and the Residual

### §6.1 The 8.7 ppb Residual

As shown in §3.2, identity (1) predicts $\alpha$ to **$8.7 \times 10^{-8}$ (8.7 ppb)** accuracy. Since experimental precision is $\sim 10^{-10}$, the residual is three orders of magnitude above experimental noise and clearly significant. At the same time, this is approximately **1/3000** of the Eddington-style integer-fitting deviation ($\alpha^{-1} \to 137$ exactly, deviation $2.6 \times 10^{-4}$): a strikingly precise agreement.

Candidate origins of the residual:

1. **Higher-order corrections**: extension to $1 = 137\alpha + V_4(1)\alpha^2 + c_3 \alpha^3 + \cdots$. Reverse-fitting from observed $\alpha$ gives residual $\sim 6 \times 10^{-10}$, divided by $\alpha^3 \approx 3.9 \times 10^{-7}$, yielding $c_3 \approx 1.6 \times 10^{-3}$. This is far smaller than $V_4(1) \approx 4.93$ and does not directly correspond to an obvious geometric quantity such as $\pi^2/2$. **Geometric identification of $c_3$ remains open**.

2. **Corrections to $V_4(1)$**: the basic quantity may itself need higher-order $\pi$-corrections (e.g. $\pi^2/2 \to \pi^2/2 + \delta$).

3. **Modification of $N(1) = 137$**: discrete count, so does not admit continuous deformation. Unlikely.

4. **Different mechanism**: identity (1) is approximate; the true equation has a different form.

Quantitatively, the residual is $|\alpha_{\text{pred}} - \alpha_{\text{obs}}| \approx 6.3 \times 10^{-10}$. This is three orders of magnitude smaller than $V_4(1) \cdot \alpha^3 \approx 1.9 \times 10^{-6}$, but larger than $V_4(1) \cdot \alpha^4 \approx 1.4 \times 10^{-8}$. This suggests **an $\alpha^3$-order geometric correction with a coefficient much smaller than $V_4(1)$**.

### §6.2 The Privilege of $R = 3$

The identity combines a $R = 1$ quantity ($V_4(1)$) and a $R = 3$ quantity ($N(1) = 137$). Why these two scales pair up is unresolved:

- $R = 1$: unit ball (smallest scale)
- $R = 3$: diameter $6 = 2 \times \sqrt{4}$ = (4D unit hypercube diagonal) $\times 3$. Smallest non-trivial complete packing?

In §2.1 the packing problem yields $N(0) = 1$ (trivial), $N(1) = 137$, $N(2) = 1545$, $N(3) = 7281$, …Why the electron corresponds to $N(1)$ (with possibly the muon at $N(2)$ and tau at $N(3)$ in a generation hierarchy) is for separate study.

### §6.3 Generation Independence

Electron, muon, and tau share the same charge, so $\alpha$ is generation-universal. The identity predicts a generation-universal $\alpha$, consistent with this. The mass hierarchy $(m_e, m_\mu, m_\tau)$ is not addressed and lies outside this paper's scope.

### §6.4 Relation to the Eddington Trap

As mentioned in §1, Eddington-style number theory failed by attempting to fit $\alpha^{-1}$ with **integers alone**. The present identity:

1. Has the structure **integer 137 + correction containing $\alpha$**, hence is a self-consistent equation
2. Both the integer 137 and the correction coefficient $\pi^2/2$ are derived **directly and independently from 4D geometry** (§2)
3. The correction has the structural form of a perturbative QFT expansion
4. Candidate physical mechanisms are proposed from a thought experiment (§5)

These are decisive differences from Eddington. Yet whether the identity itself reflects geometric necessity or numerical coincidence remains an open question. This paper does not adjudicate, only presents the observed structure.

### §6.5 Topological Duality and Node Interpretation (Future Outlook)

The identity (1) combines a discrete count $N(1)$ and a continuous volume $V_4(1)$. At first glance these have different dimensions: $N(1)$ is a dimensionless count, while $V_4(1)$ formally carries dimension $[L^4]$. This section sketches an alternative reading that addresses this mismatch through the Schläfli duality between the 4-dimensional hypercube (tesseract) and the 4-dimensional cross-polytope (16-cell).

**Re-reading via duality.** Each cell (3-volume element) of the tesseract corresponds to a vertex of its dual ($f_k$(tesseract) $= f_{4-k}$(16-cell)). The 137 unit cubes inscribed in the radius-3 4D ball can therefore be re-read on the dual lattice as **137 vertices (interaction nodes)**. Under this interpretation $N(1)$ becomes a purely dimensionless topological quantity, and the dimension-mismatch concern dissolves.

**Correspondence with QED vertex factors.** In quantum field theory, Feynman rules attach a factor of $\sqrt{\alpha}$ to each interaction vertex of an amplitude (and $\alpha$ to each vertex of a probability). Reading the perturbative form $1 = 137\alpha + V_4(1)\alpha^2$ in this light yields a structural correspondence:

| Term | Standard QFT reading | Dual-model reading |
|---|---|---|
| $137\alpha$ | tree-level vertex contribution (vertex × $\alpha$) | coupling at 137 nodes of the dual lattice |
| $V_4(1)\alpha^2$ | one-loop correction (two-vertex path) | boundary self-correlation, two-node loop |
| Total = 1 | unitarity / normalisation | normalisation of the interaction graph |

If $V_4(1)$ is likewise re-interpreted as a discrete information measure (boundary-lattice node density, or a count of spin-network states in the LQG sense), the identity becomes a fully dimensionless "normalisation condition for an interaction graph".

**Caveats.** This duality reading is a **motivational argument**, not a rigorous derivation of QED Feynman rules from the 4D packing problem. Concretely unresolved:

- What physical interaction process the 137 dual vertices represent
- The structure of the "lines" (propagators) connecting these nodes
- The relationship between the implied UV cutoff $\Lambda$ (estimated $\sim 0.6$ MeV from §5.1 Candidate A) and the dual-lattice spacing

Even so, the fact that each term of the identity parallels the canonical structure of perturbative QFT (Feynman vertex rules) provides indirect circumstantial support for a geometric origin of $\alpha$. Affinity with spin networks (LQG) and lattice gauge theory may elevate the identity from a numerical coincidence to a starting point of a programme that **reconstructs interaction theory from discrete structure**.

---

## §7. Conclusion

From two elementary geometric quantities—$N(1) = 137$, the count of integer-centred unit cubes fully contained in a 4D ball of radius 3, and $V_4(1) = \pi^2/2$, the volume of the 4D unit ball—a self-consistent algebraic equation
$$ \alpha^{-1} = N(1) + V_4(1) \cdot \alpha $$
is constructed whose solution predicts the observed fine-structure constant to **8.7 ppb ($\sim 8.7 \times 10^{-8}$)** accuracy. This is approximately 1/3000 the deviation of Eddington-style integer-fitting: a strikingly precise agreement.

This identity is:

1. **Mathematically**: independently verifiable via the elementary calculations in §2 and the numerical check in §3.
2. **Structurally**: of perturbative form $1 = 137\alpha + V_4(1)\alpha^2$, analogous to renormalisation in QFT.
3. **Physically**: readable as inside (137 packed cubes) / outside (66\% boundary) decomposition (§4), with multiple candidate mechanisms (§5).

Whether the identity reflects geometric necessity or a mechanism-specific coincidence is unresolved. The contribution of this paper is not a verdict but the presentation of a self-consistent equation, built from 4D geometric quantities not borrowed from elsewhere, that agrees with observation to 8.7 ppb.

---

## References

1. Bekenstein, J. D. "Black Holes and Entropy." *Phys. Rev. D* **7**, 2333 (1973).
2. Hawking, S. W. "Particle Creation by Black Holes." *Commun. Math. Phys.* **43**, 199 (1975).
3. Tangherlini, F. R. "Schwarzschild Field in $n$ Dimensions and the Dimensionality of Space Problem." *Nuovo Cimento* **27**, 636 (1963).
4. Eddington, A. S. "The Charge of an Electron." *Proc. Roy. Soc. London* **A122**, 358 (1929).
5. Wyler, A. "Les groupes des potentiels de Coulomb et de Yukawa." *C. R. Acad. Sci. Paris* **A269**, 743 (1969).
6. Mohr, P. J., Newell, D. B., Taylor, B. N. "CODATA Recommended Values of the Fundamental Physical Constants: 2018." *Rev. Mod. Phys.* **93**, 025010 (2021).

---

## Appendix: Numerical Reproduction

### A.1 Calculator-level (4-digit rounding)

1. $V_4(1) = \pi^2/2 \approx 4.9348$
2. $N(1) = 137$
3. Quadratic $4.9348\, \alpha^2 + 137\, \alpha - 1 = 0$
4. Discriminant $D = 137^2 + 4 \times 4.9348 \approx 18788.7$
5. $\sqrt{D} \approx 137.072$
6. $\alpha \approx 0.0072974$
7. $\alpha^{-1} \approx 137.036$ (matches observation to 4 digits)

### A.2 High precision (Python / mpmath, 50 digits)

```python
from mpmath import mp, mpf, pi, sqrt
mp.dps = 50

V4 = pi**2 / 2          # 4.9348022005446793094172454999380755676568497036204
N1 = mpf(137)
disc = N1**2 + 4 * V4   # 18788.739208802178717237668981999752302270627398814
sqrt_disc = sqrt(disc)  # 137.07202197677751743057523907144419105450522908974
alpha = (-N1 + sqrt_disc) / (2 * V4)
                        # 0.0072973519353590299272569788872011757597089129229091
alpha_inv = 1 / alpha   # 137.03601098838875871528761953572209552725261454479
```

CODATA 2018: $\alpha = 7.2973525693 \times 10^{-3}$, $\alpha^{-1} = 137.035999084$

| Quantity | High-precision predicted | CODATA 2018 | Relative error |
|---|---|---|---|
| $\alpha$ | $7.29735194 \times 10^{-3}$ | $7.29735257 \times 10^{-3}$ | $\mathbf{8.7 \times 10^{-8}}$ (8.7 ppb) |
| $\alpha^{-1}$ | 137.03601099 | 137.03599908 | $8.7 \times 10^{-8}$ |

The main claim is reproducible at calculator level (approximate) and trivially extensible to high precision in any computational environment.
