# QCD Expressed in Octonions
## ― Writing the Lagrangian of Quantum Chromodynamics As It Stands in Octonions, as a Reference for Checking Answers ―

**Author:** Noriaki Kihara (WF System Co., Ltd.)  ORCID 0009-0004-6753-4020
**Date:** 2026-11-19    **Version:** v1.0
**Version DOI:** 10.5281/zenodo.24369466    **Concept DOI:** 10.5281/zenodo.24369465
**License:** CC BY 4.0

---

## Abstract

The author has been working toward the goal of reading out spacetime, particles, and interactions from the relations among the states of a finite complex Hilbert space, without placing a background spacetime. Quantum chromodynamics (QCD), however, had not been taken up in that work. The author attended a lecture by Professor Tetsuya Onogi [X1] in a public lecture series of the Graduate School of Science, The University of Osaka. The lecture was easy to follow for a general audience, and yet it explained the difficult concepts of QCD carefully, without blurring them. The general form of the QCD Lagrangian left a particularly strong impression. The author therefore first considered whether this formula can be expressed in octonions as it stands.

When one imaginary unit of the octonions is fixed, the octonions split into one complex number and three complex numbers, and the automorphisms of the octonions that do not move that unit form SU(3). This is the known construction of Günaydin and Gürsey [X2]. This paper re-checks that construction with the author's multiplication table, and then (1) writes the QCD Lagrangian in the language of octonions, (2) expands the gluon-only term in octonion components and writes it down as seven conditions, (3) gives the correspondence with the reading of the octonions that the author has used so far, $(t, xyz)(R, Q_1Q_2Q_3)$, and (4) sorts what came out of the octonion side in this rewriting from what was still placed from outside.

This paper contains no new theorem. It does not derive QCD from the octonions. Its purpose is to prepare a reference for checking answers when, in future work, the author examines whether spacetime, particles, and interactions can be read out from the relations among states expressed in octonions without placing a background. All formulas in the text were checked numerically with the attached program.

**Keywords:** quantum chromodynamics, Lagrangian, octonions, SU(3), $G_2$, lattice gauge theory, reference for checking answers

---

## Introduction: The Lecture That Prompted This Paper

This paper grew out of attending the lecture "How Far Can the World of Elementary Particles Be Reproduced on a Computer?" [X1], given by Professor Tetsuya Onogi (Department of Physics, particle theory) on November 18, 2026, in "Science Night 2026," a public lecture series of the Graduate School of Science, The University of Osaka.

In the lecture, it was shown that quantum chromodynamics (QCD), the theory of the strong force, is written in the following single line.

$$
\mathcal{L}_{\rm QCD} = -\frac{1}{4g^2}F^a_{\mu\nu}F^{a\,\mu\nu} + \sum_q \bar q\,(i\gamma^\mu D_\mu - m_q)\,q
$$

The only quantities that can be adjusted in this formula are the strength of the force $g$ and the quark masses $m_q$. As a method for solving this formula, the lecture explained lattice QCD [X4], in which spacetime is approximated by a lattice, quarks are placed on the lattice sites, and gluons are placed on the links that connect the sites.

The lecture was delivered in a way that was easy to follow for a general audience, and yet the difficult concepts of QCD were not blurred; each was explained carefully, one by one. Thanks to this, it was very easy to understand even for the author, an independent researcher and beginner without professional training in the field. Above all, the general form of the QCD Lagrangian left a strong impression. The author therefore considered whether this formula can be expressed in octonions as it stands. This paper is the result.

---

## 1. Motivation and Position of This Paper

The author's research plan aims to read out spacetime, particles, and interactions solely from a set of states and their relations, without placing a background spacetime, a ruler, or the kinds of particles from outside. The Standard Model and general relativity are not things to be replaced; they are the reference for checking answers, which must be reproduced. So far, the author has written a map from a set of complex numbers to spacetime (1+3) with quaternions [S4], and has examined the reading (1+3)(1+3) with octonions made of two quaternions [S2][S3].

On the other hand, the author had not worked on QCD, the theory of the strong interaction. After hearing the lecture [X1] described above, the author came to think that QCD might be approachable even without professional training.

As a first step, the author considered whether the general form of the Lagrangian shown in the lecture can be expressed in octonions without changing its form. This paper is the record of that attempt.

The position of this paper is as follows.

- **What was done:** The standard QCD Lagrangian was rewritten with the color degrees of freedom placed inside the octonions. The content of the formula is the same as standard QCD.
- **Purpose:** To make a reference for checking answers in future work. When checking whether something corresponding to QCD has been read out from a formulation that places no background, the target of the comparison is written down in advance in the language of octonions.
- **What was not done:** The formulation without a background spacetime itself. In this paper, spacetime $x^\mu$ and spinors are used exactly as in standard QCD (§9).

The classification of the claims in each section is summarized in "Classification" at the end.

---

## 2. Starting Point: The QCD Lagrangian

The starting point is the following formula shown in the lecture [X1].

$$
\mathcal{L}_{\rm QCD} = -\frac{1}{4g^2}F^a_{\mu\nu}F^{a\,\mu\nu} + \sum_q \bar q\,(i\gamma^\mu D_\mu - m_q)\,q \tag{2.1}
$$

The first term is the gluon-only term, and the second term is the term for quarks and their interaction. $g$ is the coupling strength, $m_q$ are the quark masses, and the sum runs over the kinds of quarks.

In what follows, the gauge field is written as $A_\mu = -iA^a_\mu\lambda^a/2$ ($\lambda^a$ are the Gell-Mann matrices), the covariant derivative as $D_\mu = \partial_\mu + A_\mu$, and the field strength as $F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu + [A_\mu, A_\nu] = -iF^a_{\mu\nu}\lambda^a/2$.

---

## 3. Splitting the Octonions into "One Complex Number and Three Complex Numbers"

### 3.1 Multiplication Table

Let the basis of the octonions $\mathbb{O}$ be $1, e_1, \dots, e_7$, and define the product by

$$
e_i e_j = -\delta_{ij} + f_{ijk}\,e_k \tag{3.1}
$$

Here $f_{ijk}$ is totally antisymmetric and equals $+1$ for $(ijk) = (123)$, $(145)$, $(176)$, $(246)$, $(257)$, $(347)$, $(365)$. This product satisfies the multiplicativity of the norm, $|xy| = |x||y|$, and does not satisfy the associative law.

### 3.2 Fixing One Imaginary Unit

Choose one imaginary unit and fix it. In this paper $e_7$ is chosen, and multiplication by $e_7$ from the left is regarded as the complex number $i$. In the table (3.1), $e_7e_1 = -e_6$, $e_7e_2 = e_5$, and $e_7e_3 = e_4$, so the octonions split as

$$
\mathbb{O} = \mathbb{C} \oplus \mathbb{C}^3,\qquad
\mathbb{C} = \langle 1, e_7\rangle,\qquad
\mathbb{C}^3 = \langle e_1, e_2, e_3\rangle_{\mathbb C} \tag{3.2}
$$

The real six-dimensional part consists of three planes, $(e_1, e_6)$, $(e_2, e_5)$, and $(e_3, e_4)$. An element of $\mathbb{C}^3$ is written as $\psi = \psi^1e_1 + \psi^2e_2 + \psi^3e_3$, where the complex coefficients $\psi^a$ multiply from the left. In what follows, an octonion is written as $x = a + \psi$ ($a \in \mathbb{C}$, $\psi \in \mathbb{C}^3$), and its two components are denoted by $[x]_{\mathbb C} = a$ and $[x]_{\mathbb C^3} = \psi$.

### 3.3 The Product in Components

The octonion product of $x = a + \psi$ and $y = b + \chi$ takes the following form.

$$
(a + \psi)(b + \chi) = \underbrace{\left(ab - \chi^\dagger\psi\right)}_{\mathbb{C}}
\;+\; \underbrace{\left(a\chi + \bar b\,\psi + \bar\psi\times\bar\chi\right)}_{\mathbb{C}^3} \tag{3.3}
$$

Here $\chi^\dagger\psi = \sum_a \bar\chi^a\psi^a$ is the Hermitian inner product on $\mathbb{C}^3$, and $(\bar\psi\times\bar\chi)^c = \epsilon_{abc}\,\bar\psi^a\bar\chi^b$. Both the inner product and $\epsilon_{abc}$ appear inside the product of two octonions. The inner product can be taken out as the $\mathbb{C}$ component of an octonion product,

$$
\psi^\dagger\chi = \left[\chi\,\psi^*\right]_{\mathbb C} \tag{3.4}
$$

where $\psi^*$ is the conjugate as an octonion.

---

## 4. The Automorphisms That Do Not Move $e_7$ Form SU(3)

The transformations that preserve the octonion product (the automorphisms) form the 14-dimensional group $G_2$. Among them, those that do not move $e_7$ are

$$
\{\,g \in G_2 \mid g(e_7) = e_7\,\} \cong SU(3) \tag{4.1}
$$

which is 8-dimensional. These transformations do nothing on the $\mathbb{C}$ part and act on the $\mathbb{C}^3$ part as $3\times3$ unitary matrices (with determinant 1). Conversely, the map $D^a$ made from the Gell-Mann matrices by

$$
D^a\big|_{\mathbb{C}^3} = -\tfrac{i}{2}\lambda^a,\qquad D^a\big|_{\mathbb{C}} = 0 \tag{4.2}
$$

is an infinitesimal transformation that preserves the octonion product.

Under this SU(3), the eight components of the octonions split as $1 + 1 + 3 + \bar 3$. The eight gluons are not the eight components of the octonions themselves; they are on the side of the transformations that act on the octonions.

The above is the construction shown by Günaydin and Gürsey [X2], and the general properties of the octonions are summarized in [X3].

---

## 5. Writing the Lagrangian in Octonions

Using §3 and §4, the fields of QCD can be re-read as follows.

| | Standard way of writing | Reading in octonions |
|---|---|---|
| Quark $q$ | A field with three colors | A field whose components each take values in $\mathbb{C}^3\subset\mathbb{O}$ |
| Gluon $A_\mu$ | The SU(3) gauge field | A field taking values in the transformations of the octonions that do not move $e_7$ |
| Sum over colors | $\psi^\dagger\chi$ | $[\chi\psi^*]_{\mathbb C}$ |

Then (2.1) can be written as

$$
\mathcal{L}_{\rm QCD} = \frac{1}{4g^2}\,{\rm Tr}_{\mathbb O}\!\left(F_{\mu\nu}F^{\mu\nu}\right)
+ \sum_q \sum_\alpha \left[\,\left(\gamma^0 (i\gamma^\mu D_\mu - m_q)\,q\right)_\alpha \left(q_{\alpha}\right)^{*}\,\right]_{\mathbb C} \tag{5.1}
$$

Here $\alpha$ is the spinor index, and each component $q_\alpha$ is an element of $\mathbb{C}^3$. By (3.4), the second term is equal to $q^\dagger\gamma^0(i\gamma^\mu D_\mu - m_q)q = \bar q(i\gamma^\mu D_\mu - m_q)q$. The fields are treated here as classical fields. ${\rm Tr}_{\mathbb O}$ is the trace as an $8\times8$ matrix acting on the eight components of the octonions, and

$$
{\rm Tr}_{\mathbb O}(F_{\mu\nu}F^{\mu\nu}) = -F^a_{\mu\nu}F^{a\,\mu\nu} \tag{5.2}
$$

so (5.1) is the same as (2.1). The form of the formula does not change; the rewriting is such that SU(3) appears as a property of the octonions.

The gluon commutator $[A_\mu, A_\nu]$ is a commutator of matrices (transformations) and is not affected by the fact that the octonion product does not satisfy the associative law.

---

## 6. Expanding the First Term in Octonion Components

The theory with the first term alone (the theory of gluons only) is a complete theory in itself. This term is written down in octonion components.

The field strength $F_{\mu\nu}$ is represented as an antisymmetric real $6\times6$ matrix $F^{ij}_{\mu\nu}$ acting on $e_1, \dots, e_6$. $F^{ij}$ is the component that represents a rotation within the plane spanned by $e_i$ and $e_j$. The first term then becomes

$$
\mathcal{L}_1 = -\frac{1}{4g^2}\sum_{i,j=1}^{6} F^{ij}_{\mu\nu}F^{ij\,\mu\nu} \tag{6.1}
$$

However, $F^{ij}$ is not a free $6\times6$ antisymmetric matrix (15 components); it satisfies the seven conditions determined by the multiplication table,

$$
f_{ijk}F^{jk} = 0\quad(i = 1, \dots, 7),\qquad F^{7j} = 0 \tag{6.2}
$$

Written out with the table (3.1), they are

$$
\begin{aligned}
&F^{23} + F^{45} = 0,\quad F^{13} = F^{46},\quad F^{12} = F^{56},\quad F^{15} + F^{26} = 0,\\
&F^{14} + F^{36} = 0,\quad F^{24} = F^{35},\quad F^{16} = F^{25} + F^{34}
\end{aligned}
$$

15 − 7 = 8, which agrees with the eight gluons.

These seven conditions are the table of signs of the product copied as it is. For example, the products that make $e_1$ are the three pairs $e_2e_3 = +e_1$, $e_4e_5 = +e_1$, and $e_6e_7 = -e_1$, and this corresponds to the first condition $F^{23} + F^{45} - F^{67} = 0$ (with $F^{67} = 0$).

Geometrically, this can be read as follows.

- The six dimensions split into three planes, $(e_1, e_6)$, $(e_2, e_5)$, and $(e_3, e_4)$. The operation of multiplying by $e_7$ rotates each of the three planes by 90°.
- The eight gluons consist of two rotations within the planes and six rotations that turn two planes in a linked way. There are three rotations within the planes, but by the seventh condition only two of them are independent.

---

## 7. The Form on the Lattice

The form on the lattice explained in the lecture [X1] (lattice gauge theory [X4]) can also be written in octonions.

The link variable $U$ is an element of SU(3), and by §4 it is an automorphism of the octonions that does not move $e_7$. Since $U$ is the identity on $\mathbb{C}$ and a $3\times3$ unitary matrix on $\mathbb{C}^3$, we have ${\rm Tr}_{\mathbb O}\,U = 2 + 2\,{\rm Re}\,{\rm tr}_3 U$, and the Wilson action ($\beta = 6/g^2$) can be written as

$$
S_G = \beta\sum_P\left(1 - \tfrac13\,{\rm Re}\,{\rm tr}_3 U_P\right) = \frac{1}{g^2}\sum_P {\rm Tr}_{\mathbb O}\!\left(1 - U_P\right) \tag{7.1}
$$

Here $U_P$ is the product of the links around the smallest square of the lattice.

---

## 8. Correspondence with the Author's Earlier Reading

The author has so far read the octonions as a pair of two quaternions, $o = q + p\,\ell$ [S3][S4][S5]. In this reading, $q = t + x\mathbf{i} + y\mathbf{j} + z\mathbf{k}$ is the external (spacetime) side and $p = R + Q_1\mathbf{i} + Q_2\mathbf{j} + Q_3\mathbf{k}$ is the internal side.

Taking $\mathbf{i},\mathbf{j},\mathbf{k} = e_1,e_2,e_3$ and $\ell = e_7$, and multiplying $p$ by $\ell$ from the right, the components in the splitting of §3 become

$$
[o]_{\mathbb C} = t + iR,\qquad [o]_{\mathbb C^3} = \mathbf{x} - i\mathbf{Q} \tag{8.1}
$$

The component without color is the pair $(t, R)$, and the component with color is the pair $(\mathbf{x}, \mathbf{Q})$.

In the splitting of (8.1), the SU(3) transformations rotate $\mathbf{x}$ on the external side and $\mathbf{Q}$ on the internal side together as a pair. This, however, is how things look when one chooses the splitting that places $\ell$ on the internal side. For example, the center of SU(3) (the transformations of order 3) does not mix the external side and the internal side in a splitting that takes a quaternion subalgebra containing $\ell$ as the external side. This point has been checked numerically in the author's working memo [S5]. Which splitting is selected is not determined.

When more is fixed inside the octonions, the way the eight components split and the automorphisms that remain change as follows [S1][S5][X3].

| What is fixed | Splitting of the 8 components | Remaining automorphisms |
|---|---|---|
| Nothing | 1 + 7 | $G_2$ (14-dimensional) |
| One imaginary unit | 2 + 6 | SU(3) (8-dimensional) |
| One quaternion subalgebra (as a whole) | 4 + 4 | SO(4) (6-dimensional) |

§3 to §7 of this paper apply the case of the second row to the formula of QCD.

---

## 9. What Came from the Octonion Side and What Was Still Placed from Outside

In the rewriting of this paper, we sort what appeared as a property of the octonions from what was placed from outside exactly as in standard QCD.

| What appeared from the octonion side | What was still placed from outside |
|---|---|
| The three components of color ($\mathbb{C}^3$) | Spacetime $x^\mu$ and the derivative $\partial_\mu$ |
| The SU(3) of color (the automorphisms that do not move $e_7$) | The signature of the metric |
| The eight gluons (15 − 7) | Spinors and $\gamma^\mu$ |
| The inner product over colors and $\epsilon_{abc}$ (components of the product) | The kinds of quarks and their masses $m_q$ |
| The component without color ($\mathbb{C}$) | The coupling strength $g$ |
| | The lattice (in the case of §7) and the way of taking the sum as a quantum theory |

Only the part concerning color appeared from the octonion side. Everything concerning spacetime is placed from outside. In addition, the choice of one imaginary unit $e_7$ is itself made by hand in this paper.

Care is also needed regarding the identification of the $i$ of quantum theory with the operation of multiplying by $e_7$ from the left. The operation of multiplying by $e_7$ from the right is antilinear ($\psi a = \bar a\psi$), so left and right are not the same.

---

## 10. As a Reference for Checking Answers

When, in future work, the author examines the readout of interactions from the relations among states expressed in octonions without placing a background, the following points will be used as the reference for comparison for the strong interaction.

1. Whether the three components corresponding to color and the component without color appear in the splitting (3.2).
2. Whether the transformations that carry the interaction appear as eight transformations satisfying the seven conditions (6.2).
3. Whether the combinations without color appear in the form of the inner product and $\epsilon_{abc}$ of (3.3).
4. Whether the gluon-only term appears in the form (6.1) or in the form (7.1).
5. Whether the adjustable quantities are only the coupling strength and the quark masses.

These are the content of standard QCD written in the language of octonions, and they are not new predictions. If a point of disagreement appears as a result of the comparison, it will be recorded.

---

## 11. What Is Not Claimed

- It is not claimed that QCD has been derived from the octonions.
- The construction of §3 to §5 is due to Günaydin and Gürsey [X2] and is not a new result. The expansion in §6 and the rewriting in §7 also write down the same content in components.
- It is not claimed that the Standard Model is replaced.
- It is not claimed that a formulation without a background spacetime has been given. This paper places spacetime from outside.
- It is not claimed that the correspondence (8.1) is the physically correct reading. It is a table of correspondence that sets the author's notation and the splitting of §3 side by side.
- Nothing is stated about the kinds of quarks, their masses, or the electroweak interaction.

---

## Acknowledgments

This paper owes a great deal to the public lecture [X1] by Professor Tetsuya Onogi. The author is deeply grateful for a lecture that was both easy to follow and sincere. The content of this paper is the author's own attempt, and the author alone is responsible for any errors.

---

## References

### Self-Citations (Concept DOI, `10.5281/zenodo.<number>`)

[S1] Noriaki Kihara, "Projection Axis Selection and Covariance of Conservation Laws ― Two Radii Carried by All Waves, Projection Asymmetry of 1-, 2-, and 4-Dimensional Subalgebras, and Reduction to Minimal Systems ―," 2026-11-16. Concept DOI 10.5281/zenodo.24273723.  
[S2] Noriaki Kihara, "Emergence of Spacetime and Elementary Particles from a Bounded Hilbert Space and Octonions ― From Photon Geometry and Gnomonic Central Projection to Open Problems in Background-Free Quaternionic and Octonionic Dynamics ―," 2026-11-16. Concept DOI 10.5281/zenodo.24267347.  
[S3] Noriaki Kihara, "Did We Understand Complex Numbers Correctly? Where the Amplitude Lives, the Imaginary Reading, and Attribution across $\mathbb{C} \to \mathbb{H} \to \mathbb{O}$," 2026-10-23. Concept DOI 10.5281/zenodo.23585436.  
[S4] Noriaki Kihara, "A Map from a Set of Complex Numbers to Spacetime (1+3): Readout by Quaternions and Partial Projection," 2026-10-18. Concept DOI 10.5281/zenodo.23451494.  
[S5] Noriaki Kihara, "Octonion Closure and the Decomposition of Relations" (in Japanese), working memo and numerical checks, 2026-10-19. GitHub: https://github.com/WurabeSeiji/ai-chat-logs-open , folder `次元の生成構造/八元数閉塞と関係の分解_20261019/` (not registered on Zenodo).

### External References

[X1] Tetsuya Onogi, "How Far Can the World of Elementary Particles Be Reproduced on a Computer?" (in Japanese), public lecture series Science Night 2026, Graduate School of Science, The University of Osaka, November 18, 2026.  
[X2] M. Günaydin and F. Gürsey, "Quark structure and octonions," *J. Math. Phys.* 14, 1651–1667 (1973).  
[X3] J. C. Baez, "The Octonions," *Bull. Amer. Math. Soc.* 39, 145–205 (2002). DOI 10.1090/S0273-0979-01-00934-X. arXiv:math/0105155  
[X4] K. G. Wilson, "Confinement of quarks," *Phys. Rev. D* 10, 2445–2459 (1974). DOI 10.1103/PhysRevD.10.2445

---

## Appendix: Verification Program

`verify_qcd_octonion_paper_v1.py` (Python 3, NumPy). It checks (3.1) to (8.1) of the text with values made from random numbers. The random seed is fixed, and all 25 items hold. The output is saved in `results/verify_qcd_octonion_paper_v1.txt`.

---

## Classification

| Category | Content |
|---|---|
| Identities / known mathematics | The multiplication table (3.1) and the multiplicativity of the norm / the splitting (3.2) / the product in components (3.3)(3.4) / (4.1)(4.2) that the automorphisms not moving $e_7$ form SU(3) [X2][X3] / the trace relation (5.2) / the seven conditions (6.2) and 15 − 7 = 8 / the rewriting of the Wilson action (7.1) / the table in §8 (subgroups of $G_2$) |
| Rewriting (the content is the same as standard QCD) | (5.1) the Lagrangian expressed in octonions / (6.1) the component expansion of the first term |
| Readings / interpretations | The geometric reading in §6 (three planes and linked rotations) / (8.1) the correspondence with the author's notation |
| Organization | The sorting in §9 / the reference for comparison in §10 |
| Open / outside the scope of this paper | Reading out spacetime, spinors, and the signature of the metric from the octonion side / what selects the imaginary unit $e_7$ and the splitting / the meaning of the difference between multiplying from the left and from the right / the kinds of quarks, their masses, and the electroweak interaction |

---

## Version History

- **v1.0 (2026-11-19):** First version. Organized into the form of a paper from a study memo written after attending the public lecture [X1] at The University of Osaka.
- Environment for organizing and checking: Claude (Claude Code) was used. The component expansion in §6 is based on an examination with claude.ai and was confirmed by the author. All formulas in the text were checked with the attached program.
