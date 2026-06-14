# Paper 11: The Necessity of Dimension

## Both-Endpoint Theorem, Survival Interval $[0,4]$, Mod 8 Selection Theorem, Developmental Selection, and Four-Legs Theorem

Author: Noriaki Kihara  
Affiliation: WF System Co., Ltd.  
ORCID: 0009-0004-6753-4020  
Version: v0.3 (content and results unchanged from v0.2; nomenclature unified to the displacement-record definition of Paper 0.5 — ledger/register sorted into displacement record / conserved quantity / invariant / accounting; integer ledger and the rigid 1D ledger renamed to displacement record)  
Date: June 2026  
DOI (this version): 10.5281/zenodo.20689962  
Concept DOI: 10.5281/zenodo.20640466  
License: CC BY 4.0  

* * *

## Abstract

In the preceding papers of the reciprocal dual model, the spatial dimension $d=4$ (the lattice $\mathbb{Z}^4$) was an input. This paper organizes, as several mutually independent selection mechanisms, the extent to which $d=4$ is **necessary** under the three assumptions ($\nu\lambda=1$, the zero point $1/2$, and the asymmetric one bit).

First, the **both-endpoint theorem**: the condition for the half-wavelength censorship that protects composites (Paper 9) is, via the closed form of the maximal anharmonic shift of the container spectrum, $\Delta\nu_1(d)=(\sqrt d-1)^2/2$, equivalent to

$$
\Delta\nu_1(d)\le\tfrac12\iff d\in[0,4]
$$

and **saturation (equality) holds exactly at the two endpoints $d\in\{0,4\}$ only**. For $d\ge5$ the shift of the lowest mode exceeds half a quantum, a continuous decay channel opens, and no stable composite can exist. The interior ($0<d<4$) is stable with genuine margin, and this is the existence condition for a developmental path that "grows from zero with 4 as the ceiling."

Second, the **mod 8 selection theorem**: since the dressed label of a $d$-axis lattice satisfies $4s=\sum(2|k_i|+1)^2\equiv d\pmod 8$, the label is an integer only for $d\equiv0\pmod4$, and **it is an odd integer (the foundation of the classification theorem for single states, Paper 5) only for $d=4$** ($d=8$ yields even labels and single states disappear). Sectors that can carry records (an integer-valued displacement record) are restricted to 4 axes — a second selection mechanism independent of the censorship ceiling.

Third, **fixing the normalization**: hierarchical relativity forces the censorship criterion to be level-local (the ladder unit of each composite itself), and the comparison ladder $\ell+\tfrac{d-1}2$ is the unique canonical equally spaced ladder by the exact identity of the spherical spectrum $\omega_\ell^2=L^2-(\tfrac{d-1}2)^2$. The remaining normalization freedom reduces to the single condition of the curvature coupling $\xi$ of the wave operator, and we take minimal coupling $\xi=0$ as the default (we also make explicit that with conformal coupling the dimension selection disappears).

Fourth, the **four-legs theorem (conditional necessity)**: if a stable composite exists, then its state space must be discrete (the discrete spectrum of a compact container), bounded without boundary, of positive curvature, and with $d\le4$. Flat infinite and negatively curved spaces have continuous spectra, hence no discrete ladder for the censorship to protect, and cannot support stable composites.

Fifth, **developmental selection and the selector from below**: since the capacity at fixed budget $S$ grows rapidly with dimension as $N^{(d)}(S)\approx V_d S^{d/2}$, dimensional growth is driven by the same counting principle as splitting. On the other hand, the actualization of an axis carries a zero-point cost, giving the budget bound $d\le4S$, and the only dimension in which the minimal cell (pure zero point) coincides with the dual fixed point $s=1$ is $d=4$ (the self-dual seed). The ceiling from above (censorship) and the selector from below (self-duality) coincide at $d=4$.

These constitute a **convergence** of several independent selection mechanisms; this paper makes explicit the preconditions of each mechanism and identifies the remaining tasks for removing dimension from the inventory of axioms.

* * *

## Keywords

spacetime dimension, dimension selection, Weyl asymptotics, anharmonicity, mod 8, developmental selection, curvature sign, self-duality, normalization, emergence

* * *

## 1. Introduction

The question "why is space three-dimensional (spacetime four-dimensional)?" has a long history. Ehrenfest [5] argued from the relation between the inverse-square law and stable orbits, Whitrow [6] from the possibility of life, and Tegmark [7] discussed the special character of $3{+}1$ from the viewpoint of predictability. On the dynamical side, causal dynamical triangulations (CDT) have numerically exhibited the emergence of four-dimensional macroscopic structure (Ambjørn, Jurkiewicz & Loll [8]). Algebraically, Hurwitz's theorem [9,10], stating that the dimensions of associative division algebras are limited to $1,2,4$, is a representative example of exceptional structure specific to $d=4$.

The setting of this paper follows the lattice model of Paper 4 [1] and Paper 5 [2]. The approach of this paper is independent of the above, and shows that the dimension is constrained by **two internal requirements: the stability of composite states (censorship) and the arithmetic of records (mod 8)**. The relation to prior approaches (in particular the four-dimensional emergence of CDT and the intersection with the dimension sequence of division algebras) is discussed as structural correspondence in §8.

* * *

## 2. Both-Endpoint Theorem

### 2.1 Closed form

We examine the condition for the half-wavelength censorship that protects composites (Paper 9 [4]). The maximal deviation between the spectrum $\ell(\ell+d-1)$ of the $S^d$ container and the canonical ladder $\ell+\tfrac{d-1}2$ (whose canonicity is shown in §4) is, by the reduction with the substitution $t=\sqrt d$,

$$
\Delta\nu_1(d)=\frac{\bigl(\tfrac{d-1}2\bigr)^2}{\tfrac{d+1}2+\sqrt d}=\boxed{\ \frac{(\sqrt d-1)^2}{2}\ }
$$

($d=1:0$, $d=2:0.086$, $d=3:0.268$, $d=4:\tfrac12$, $d=5:0.764$, $d=9:2$).

### 2.2 Theorem

The censorship condition (distortion shift ≤ half the resolution $=\tfrac12$):

$$
\frac{(\sqrt d-1)^2}{2}\le\frac12\iff|\sqrt d-1|\le1\iff d\in[0,4]
$$

> **Both-endpoint theorem**: The survival interval is the closed interval $[0,4]$, and saturation (equality) holds exactly at the two endpoints $d=0$ and $d=4$ only. This is a uniqueness statement not restricted to integers, and is not a "coincidence of perfect squares." For $d\ge5$ a continuous decay channel opens, and stable composites, persistent records, and hierarchy cannot stand.

$d=1$ (the circle) has zero curvature and exactly zero shift — the first step of growth is cost-free. The shift then increases monotonically and is exactly used up at the ceiling $d=4$. **The margin at low dimensions is not a reason for elimination but a passability, and criticality is not a selection criterion but a consequence of reaching the ceiling.**

![Figure 1. Both-endpoint theorem: survival interval [0,4]](figure_paper11_1_survival_interval.png)

**Figure 1. Both-endpoint theorem: $\Delta\nu_1(d)=(\sqrt d-1)^2/2$ and the censorship limit $1/2$. The survival interval is $[0,4]$, with saturation only at the two endpoints.**


* * *

## 3. Mod 8 Selection Theorem

Since odd squares are $\equiv1\pmod 8$, the dressed label of a $d$-axis lattice satisfies

$$
4s=\sum_{i=1}^d(2|k_i|+1)^2\equiv d\pmod 8.
$$

Therefore:

| $d$ | 1 | 2 | 3 | **4** | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| Label integer | ✗ | ✗ | ✗ | **✓** | ✗ | ✗ | ✗ | ✓ |
| Parity when integer | — | — | — | **odd** | — | — | — | even |

> **Theorem**: Integer labels occur only for $d\equiv0\pmod4$. Odd integer labels — the entire foundation of the classification theorem, the splitting theorem, $Z_2$, the tiling, and the genealogy (Paper 5 [2], Paper 6 [3], Paper 7) — occur **only for $d=4$** ($d=8$ yields even labels and single states disappear).

The censorship ceiling (§2) selects $d=4$ from "stability," and the mod 8 theorem selects it from "the arithmetic of records," **independently of each other**.

* * *

## 4. Fixing the Normalization

The dimension-selection claims could depend on three normalization freedoms. This paper decomposes them, settles two, and reduces one to a named condition.

**(A) Choice of unit**: What the censorship protects is the internal coherence of a composite, so the criterion must be level-local (the ladder unit of each composite itself) (hierarchical relativity, Paper 6 [3] §5). This makes the marginal saturation (exactly $1/2$ at $\ell=1$) universal across all hierarchy levels.

**(B) Canonicity of the ladder**: The spherical spectrum satisfies the exact identity $\omega_\ell^2=L^2-a^2$ ($L=\ell+a$, $a=\tfrac{d-1}2$), and the equally spaced ladder whose asymptotic deviation vanishes is $c=a$ only. The comparison ladder is neither an approximation nor a convenience but **canonical**.

**(C) Curvature coupling**: Among the one-parameter family of wave operators $-\Delta+\xi R_{\mathrm{scal}}$ on the container, with minimal coupling ($\xi=0$) the defect $a^2$ is dimension-dependent and the entire chain of §2 holds. With conformal coupling ($\xi_c$) the defect becomes universally $\tfrac14$ in all dimensions, and **the dimension selection evaporates** (this contrast is itself diagnostic; for the context of curvature coupling see Birrell & Davies [11]). By a principle of parsimony (adding a curvature term = adding structure) we take $\xi=0$ as the default, and state the condition of the $d\in[0,4]$ theorems explicitly as "under $\xi=0$." An in-model derivation of $\xi$ (the curvature response of the zero-crossing set of logical waves) remains an open task.

* * *

## 5. Four-Legs Theorem: Trichotomy of Curvature Signs

The stability mechanism (censorship) presupposes the existence of "a discrete ladder to protect." Inspecting the three curvature signs:

- **Flat infinite $E^d$**: The spectrum is continuous. Every decay channel is continuous, and there is no discrete structure for the censorship to protect. Moreover, the anharmonicity is identically zero, so the very selection mechanism of the dimension ceiling disappears.
- **Negative curvature $H^d$**: Non-compact, with a continuous spectral band. Likewise no stable composites.
- **Positive curvature $S^d$**: Compact ⟹ discrete spectrum ⟹ there exists an object for the censorship to protect, the anharmonicity realizes the ceiling, and condensation by shared curvature (Paper 8) is possible.

> **Four-legs theorem (conditional necessity)**: If a stable composite exists, then its state space must be **discrete, bounded without boundary, of positive curvature, and with $d\le4$ (critically saturated at the ceiling)**. Positive curvature is not an aesthetic choice; it is "the only curvature sign under which the half-wavelength censorship has something to protect."

Finiteness follows from the self-consistency condition $R_c^2=S$ (finite content forces a finite positive curvature radius), and the absence of boundary follows from hierarchical relativity and relationalism (the exclusion of boundaries = absolute markers).

* * *

## 6. Developmental Selection: Division of Labor between Push-Up and Ceiling

### 6.1 Counting drive

The capacity at budget $S$ on the shifted lattice $(\mathbb{Z}+\tfrac12)^d$ is, for large $S$, $N^{(d)}(S)\approx V_d S^{d/2}$ ($V_d$ the unit-ball volume), and the gain per step is $N^{(d+1)}/N^{(d)}\approx(V_{d+1}/V_d)\sqrt S>1$. **Dimensional growth is driven by the same microcanonical counting principle as splitting and branching ratios; no new weighting principle is needed.** However, the optimal dimension by counting alone is $d^\ast\sim2\pi S\gg4$, so what stops it at 4 is not the drive but the **brake (censorship)** — the division of labor is settled: push-up = counting, ceiling = stability.

### 6.2 Reduction of dimension events and the budget bound

$\mathbb{Z}^d$ is an unexcited slice of $\mathbb{Z}^{d+1}$, so "adding a dimension" is reduced not to a new type of move but to **the first excitation of a new axis = an ordinary configuration transition**. An actualized axis pays the zero-point cost $\tfrac14$, so an object with budget $S$ can actualize at most $d\le4S$ axes — **the finiteness of dimension itself is a consequence of the zero point**.

### 6.3 The self-dual seed (selector from below)

The representative value of the minimal cell (all axes unexcited = pure zero point) is $s_{\min}(d)=d/4$. Since the self-dual point of $\nu\lambda=1$ is $s=1$,

$$
s_{\min}(d)=1\iff d=4
$$

**The only dimension in which the minimal cell sits exactly at the fixed point of the duality is 4.** If we adopt as a principle that "the undifferentiated seed should sit at the fixed point of the dual transformation" (the asymmetric bit cannot choose a side without a record), then $d=4$ comes out strictly and uniquely from below — and **meets the ceiling from above (§2) in a pincer**. The coincidence is not demanded: it is the confluence of two independent computations. Note that the selector from below does not depend on $\xi$ (its freedom from the condition of §4(C) is an important point).

![Figure 2. Three independent dimension selectors](figure_paper11_2_selectors.png)

**Figure 2. Three independent dimension selectors: (a) mod 8 arithmetic (odd integer labels only at $d=4$), (b) the self-dual seed $s_{\min}(d)=d/4=1\iff d=4$ (from below), (c) the counting drive (capacity gain at every step, push-up).**


* * *

## 7. The Emergence Sequence (Synthesis)

Combining the above, the emergence from a high-frequency undifferentiated single state can be read as the following sequence (the individual links are theorems; the unproven junctions are made explicit in §9):

$$
\text{seed}(d=0)
\;\xrightarrow{\text{counting drive}}\;
\text{dimensional growth}\ d\le4
\;\xrightarrow{\text{critical at the ceiling}}\;
d=4
\;\xrightarrow{\text{stability limit exceeded}}\;
\text{splitting cascade}
\;\xrightarrow{\text{terminal}}\;
\text{stable species}\{1,3,5\}
$$

The fragment number $n$ marks the order of appearance of gauge-invariant quantities: at $n=1$ nothing (cloud), at $n=2$ length, at $n=3$ angle, curvature, and spin (the birth point of gauge structure, Paper 10), and at the terminal stage permanent labels. **The rigid one-dimensional displacement record (the clock) exists before splitting, while space-like relational quantities (two-dimensional and higher) are born only together with splitting — in this sequence, time-like structure precedes space-like structure.**

* * *

## 8. Discussion: Convergence of Independent Selectors

None of the selection mechanisms of this paper (censorship ceiling, mod 8, self-dual seed) directly uses the exceptional structures of the lattice. On the other hand, $d=4$ possesses an independent series of exceptional structures — the maximal dimension of associative division algebras [9,10], the self-duality of $\mathbb{Z}^4$/the 24-cell, the non-simplicity of SO(4) — and CDT dynamically exhibits the emergence of four-dimensional macroscopic structure [8]. The convergence of multiple independent routes to the same value is circumstantial evidence of a stronger kind than any single argument alone. However, the proof of a common origin (a single structure generating all of these) has not been undertaken, and the convergence is recorded as a conjecture.

* * *

## 9. Scope of Claims of This Paper

What this paper claims: (1) the both-endpoint theorem (exact under $\xi=0$); (2) the mod 8 selection theorem (unconditional); (3) the three-way decomposition of the normalization and the settlement of (A)(B); (4) the four-legs theorem (conditional necessity: existence of a stable composite ⟹ discrete, bounded without boundary, positive curvature, $d\le4$); (5) the arithmetic of the counting drive, the budget bound, and the self-dual seed.

What this paper does not claim: (1) an explanation of the physical spacetime dimension ("dimension" here is the number of axes of the visible sector of the model); (2) an in-model derivation of $\xi=0$ (open task); (3) the promotion of the self-dual seed to a "principle" (it remains a candidate); (4) the determination of the interleaving order of dimensional growth drive and splitting; (5) a common origin with the Hurwitz family and CDT (conjecture).

* * *

## 10. Conclusion

In this model, $d=4$ is selected by three independent mechanisms: (i) stability (the censorship ceiling, from above), (ii) the arithmetic of records (mod 8, unconditional), and (iii) duality (the self-dual seed, from below). Including, further, the curvature sign and finiteness, the single condition "a stable composite exists" entails as necessary conditions that the state space be discrete, bounded without boundary, of positive curvature, and four-dimensional (critical). What remains for removing dimension entirely from the inventory of axioms are two points: the internal derivation of $\xi$ and the promotion of the self-dual principle.

* * *

## Appendix: Key Points of Verification

The algebra of the reduction $(\sqrt d-1)^2/2$ (the substitution $t=\sqrt d$), the mod 8 table ($d=1$–$8$), the exact flatness of $d=1$, the formal saturation of $d=0$, the numerical table of capacity gains $V_{d+1}/V_d$, and the arithmetic of $s_{\min}(d)=d/4$. Scripts and research notes are provided in the public repository (https://github.com/WurabeSeiji/ai-chat-logs-open).

* * *

## References

[1] Noriaki Kihara, Paper 4, v1.0, 2026. Concept DOI: 10.5281/zenodo.20638962.

[2] Noriaki Kihara, Paper 5, v0.2, 2026. Concept DOI: 10.5281/zenodo.20640454.

[3] Noriaki Kihara, Paper 6, v0.2, 2026. Concept DOI: 10.5281/zenodo.20640456.

[4] Noriaki Kihara, Paper 9, v0.2, 2026. Concept DOI: 10.5281/zenodo.20640462.

[5] P. Ehrenfest, "In what way does it become manifest in the fundamental laws of physics that space has three dimensions?," Proceedings of the Amsterdam Academy 20 (1917), 200–209.

[6] G. J. Whitrow, "Why physical space has three dimensions," British Journal for the Philosophy of Science 6 (1955), 13–31.

[7] M. Tegmark, "On the dimensionality of spacetime," Classical and Quantum Gravity 14 (1997), L69–L75.

[8] J. Ambjørn, J. Jurkiewicz, and R. Loll, "Emergence of a 4D world from causal quantum gravity," Physical Review Letters 93 (2004), 131301.

[9] A. Hurwitz, "Über die Composition der quadratischen Formen von beliebig vielen Variablen," Nachr. Königl. Ges. Wiss. Göttingen (1898), 309–316.

[10] J. C. Baez, "The octonions," Bulletin of the American Mathematical Society 39 (2002), 145–205.

[11] N. D. Birrell and P. C. W. Davies, *Quantum Fields in Curved Space*, Cambridge University Press, 1982.

* * *

## License

This paper is published under CC BY 4.0.  
Reuse, adaptation, translation, and citation are permitted, provided that the author name, version, publication date, and source are clearly indicated.
