# The Geometric Identity of the Zero-Square-Sum Constraint under Scale Invariance
## Isotropic cone, projective quadric, and intrinsic quantum structure — an expository note

**Author:** Noriaki Kihara<br>
**ORCID:** 0009-0004-6753-4020<br>
**Date:** July 22, 2026<br>
**Version DOI:** 10.5281/zenodo.21495306<br>
**Concept DOI:** 10.5281/zenodo.21495305<br>
**Position:** Expository note. It contains no new theorems; its contribution is the synthesis of known mathematical results and their placement onto one object.

---

## Abstract

Consider a system defined by exactly two postulates.

1. States are complex vectors $x=(x_1,\dots,x_N)\in\mathbb C^N\setminus\{0\}$ satisfying the zero-square-sum constraint $\sum_{n=1}^N x_n^2=0$.
2. Only ratios are physical: $x$ and $\lambda x$ ($\lambda\in\mathbb C^\ast$) represent the same state.

The purpose of this note is to identify the mathematical identity of the state space these two postulates define, by citation of known results alone. The conclusions are organized in five points.

1. The solution set of postulate 1 is the complex isotropic (null) cone, and after projectivization by postulate 2 the state space is a quadric hypersurface in complex projective space — the **projective quadric** $Q_{N-2}\subset\mathbb{CP}^{N-1}$. This is a compact Kähler manifold [4]
2. By the general theory of geometric quantization, a system with compact phase space has, with no room for choice, **finite-dimensional state spaces and discrete spectra** [1,2,3]. Discreteness (quantumness) is not an additional hypothesis but a structure that follows automatically from the two postulates. The system is not merely quantizable — it is born quantized
3. Writing $x_n=q_n+ip_n$, the real part of the constraint is $\sum q_n^2=\sum p_n^2$ (equipartition) and the imaginary part is $\sum q_np_n=0$ — the condition that the generator of scale transformations (the dilatation) of classical mechanics vanishes identically. Scale invariance (postulate 2) may not be an independent requirement but already present in the imaginary part of postulate 1
4. The space of polynomial functions on the cone is exactly the space of **harmonic polynomials** [6], and harmonic polynomials give the spherical harmonics — integer-labelled irreducible multiplets of the rotation group. The structure of "harmonics" appears not as a metaphor but as a theorem of the state space
5. Line bundles on compact projective manifolds are classified by integers (Chern classes) [4,5], and the exact integrality of winding-type quantities is the integrality of cohomology itself. In particular, for $N=3$ the projectivized isotropic cone is a conic isomorphic to the Riemann sphere $\mathbb{CP}^1$, and null vectors can be written as squares of spinors (Cartan's construction [7]). The state space of the minimal nontrivial system is thus isomorphic to that of a qubit (the projective state space of spin 1/2)

---

## 1. Postulates and Notation

We restate the two postulates.

> **Postulate 1 (zero-square-sum constraint)** $x\in\mathbb C^N\setminus\{0\}$, $\displaystyle\sum_{n=1}^N x_n^2=0$.
>
> **Postulate 2 (scale invariance)** Only ratios are physical: $x\sim\lambda x$ ($\lambda\in\mathbb C^\ast$).

The solution set of postulate 1,

$$
\mathcal C=\Bigl\{x\in\mathbb C^N\setminus\{0\}\ :\ \sum_{n=1}^N x_n^2=0\Bigr\},
$$

is the **isotropic cone** (null cone) of the complex quadratic form $\sum x_n^2$. Over the reals there is no nontrivial solution (a vanishing sum of real squares forces every component to vanish). A nontrivial solution requires imaginary parts; in this sense the complex structure is built into postulate 1.

For $N=2$ the cone degenerates into the two lines $x_2=\pm ix_1$, and the projectivization is two points. Throughout we consider $N\ge3$, where smooth geometry appears.

The identification $x\sim\lambda x$ of postulate 2 is the quotient by the group $\mathbb C^\ast$. Corresponding to the decomposition $\mathbb C^\ast=\{|\lambda|\}\times\{e^{i\varphi}\}$, the unreadability of absolute scale ($|\lambda|$) and the unreadability of common phase ($e^{i\varphi}$) are the modulus part and the argument part of one and the same group.

---

## 2. The State Space Is a Projective Quadric

The cone $\mathcal C$ is invariant under $\lambda x$ ($\lambda\in\mathbb C^\ast$), so the quotient

$$
Q_{N-2}\ :=\ \mathcal C/\mathbb C^\ast\ =\ \bigl\{[x]\in\mathbb{CP}^{N-1}\ :\ \textstyle\sum_n x_n^2=0\bigr\}
$$

is well defined. It is a nondegenerate quadric hypersurface in complex projective space $\mathbb{CP}^{N-1}$ — the **projective quadric** — a smooth compact complex manifold of complex dimension $N-2$ [4].

$\mathbb{CP}^{N-1}$ is a Kähler manifold with the Fubini–Study metric, and the quadric, as a complex submanifold, inherits the Kähler structure [4]. Therefore:

> The state space defined by the two postulates is a compact Kähler manifold.

This is the identification from which the note departs. Everything below is a known consequence of this one line.

---

## 3. Compactness Supplies Quantumness

The general theory of geometric quantization (Kostant [1], Souriau [2]) teaches the following. A symplectic manifold $(M,\omega)$ admits a prequantum line bundle $L$ when $[\omega/2\pi]$ is an integral cohomology class (the prequantization condition); on a Kähler manifold, choosing the holomorphic polarization, the quantum state space is constructed as the space of holomorphic sections $H^0(M,L^{\otimes k})$. If $M$ is **compact**, this space is **finite-dimensional** for every $k$.

On the quadric the restriction of the Fubini–Study form represents the hyperplane class (an integral class), so the prequantization condition can be satisfied automatically. Therefore:

> The system of the two postulates has, with no room for choice in the quantization procedure, finite-dimensional state spaces and discrete spectra.

The classical explicit model showing that quantization on a compact phase space forces discreteness and finite-dimensionality is Berezin's example of the sphere [3] (quantization of the sphere = spin systems = finite dimensions). In ordinary quantization, discreteness is obtained by introducing $\hbar$ from outside; in this system the compactness of the (projectivized) cone supplies the discreteness. In this sense the system is not merely quantizable — it is **born quantized**.

---

## 4. Real and Imaginary Parts of the Constraint: Equipartition and Dilatation

Writing $x_n=q_n+ip_n$ ($q_n,p_n\in\mathbb R$),

$$
\sum_n x_n^2=\Bigl(\sum_n q_n^2-\sum_n p_n^2\Bigr)+2i\sum_n q_np_n,
$$

so postulate 1 decomposes into two real conditions:

$$
\text{real part:}\quad \sum_n q_n^2=\sum_n p_n^2,
\qquad
\text{imaginary part:}\quad \sum_n q_np_n=0.
$$

The real part, reading $q$ as configuration and $p$ as momentum, is the harmonic-oscillator-type **equipartition** condition. The imaginary part, $D=\sum_n q_np_n$, is the classical **generator of scale transformations (the dilatation)**. The imaginary part is thus the condition that the generator of scale transformations vanishes identically: the zero-square-sum constraint has, built in from the start, a structure that kills the scale degree of freedom identically.

From this one observation follows: **scale invariance (postulate 2) may not be an independent requirement but already present in the imaginary part of postulate 1.** This note records the observation and does not decide the logical independence or dependence of the two postulates.

---

## 5. The Function Space on the Cone Is Harmonic Polynomials

Restricting polynomials on $\mathbb C^N$ to the isotropic cone $\mathcal C$ is the same as quotienting by the ideal $(\sum x_n^2)$. By the classical direct-sum decomposition

$$
\mathcal P_k=\mathcal H_k\oplus r^2\,\mathcal P_{k-2}
\qquad(\mathcal P_k:\text{homogeneous polynomials of degree }k,\ \mathcal H_k:\text{harmonic polynomials},\ r^2=\textstyle\sum x_n^2),
$$

the space of degree-$k$ functions on the cone is exactly identified with the harmonic polynomials $\mathcal H_k$ (a consequence of the vanishing of the symbol of the Laplacian on the cone; the standard reference is Stein–Weiss [6], Chapter IV). Its dimension is

$$
\dim\mathcal H_k=\binom{N+k-1}{k}-\binom{N+k-3}{k-2},
$$

and on the sphere $\mathcal H_k$ gives the spherical harmonics — the irreducible multiplets of the rotation group with integer label $k$.

The natural function system of this state space is therefore literally **harmonics**. On this geometry, the word "harmonic" is not a metaphor but a theorem of the function space.

---

## 6. Integrality: Chern Classes and Winding

Complex line bundles on a compact projective manifold are classified by the first Chern class $c_1\in H^2(M;\mathbb Z)$ — an **integer-coefficient** cohomology class (Chern [5]; textbook treatment in [4]). On the quadric, $H^2(Q_{N-2};\mathbb Z)\cong\mathbb Z$ (with the single exception $N=4$, where $Q_2\cong\mathbb{CP}^1\times\mathbb{CP}^1$ and $H^2\cong\mathbb Z^2$) [4].

Consequently, the exact integrality of "winding"-type quantities on this state space (the number of turns of the phase along a closed loop, the label of a line bundle) is not an approximate or statistical property but the integrality of the cohomology of a compact projective manifold itself. The integrality does not depend on any detail of dynamics and is invariant under continuous deformations and scale transformations.

---

## 7. The Minimal Nontrivial System $N=3$: Conic, Spinor, Qubit

For $N=3$, the projectivized isotropic cone is a smooth conic (quadric curve) in $\mathbb{CP}^2$, and a conic is isomorphic, as a rational normal curve, to $\mathbb{CP}^1$ — the Riemann sphere [4].

Moreover, null vectors can be written explicitly as squares of spinors (Cartan's construction [7]). For a spinor $(a,b)\in\mathbb C^2$, set

$$
x_1=a^2-b^2,\qquad x_2=i\,(a^2+b^2),\qquad x_3=-2ab;
$$

then $x_1^2+x_2^2+x_3^2=0$ holds identically, and conversely every isotropic vector of $\mathbb C^3$ can be written in this form. The ratio of the spinor components, $[a:b]\in\mathbb{CP}^1$, determines the state.

The state space of the minimal nontrivial system is therefore the Riemann sphere — **isomorphic to the projective state space of a qubit (a two-level system), i.e., of spin 1/2**. $SU(2)$ acts on the spinors, $SO(3)$ acts on the isotropic vectors, and the double cover $SU(2)\to SO(3)$ is implemented by Cartan's construction itself. The structures of spin and $SU(2)$ are in a position to emerge from this geometry without additional assumptions.

We note that the program of founding physical theory on projectivized null structures and spinorization is the adjacent field that Penrose's twistor theory [8] has developed for half a century. The components of this note — null cone, spinors, holomorphic structure — are all standard tools there.

---

## 8. What This Note Does Not Claim

1. All identifications in this note are **kinematic** (at the level of the state space). No correspondence with any particular dynamics, update rule, or physical system is asserted.
2. The statement of Section 4 that "scale invariance is already present in the imaginary part of the constraint" is an observation, not a decision on the logical relation of the two postulates.
3. The statement of Section 7 that "the origin of spin is in a position to emerge" is the indication of a mathematical isomorphism, not a claim of deriving physical spin.
4. This note contains no new theorems. Every mathematical claim belongs to the cited literature; the contribution is the placement and synthesis onto the system of the two postulates.

## 9. Classification of Claims

| Claim | Classification |
|---|---|
| Projectivized isotropic cone = projective quadric, compact Kähler | Application of known mathematics [4] |
| Quantization of compact phase space = finite-dimensional, discrete | Application of known mathematics [1,2,3] |
| Real part = equipartition; imaginary part = vanishing dilatation | Elementary computation |
| Scale invariance is present in the imaginary part of the constraint | Observation (indication of a possibility) |
| Function space on the cone = harmonic polynomials = harmonics | Application of known mathematics [6] |
| Exact integrality of winding = integrality of cohomology | Application of known mathematics [4,5] |
| $N=3$: conic $\cong\mathbb{CP}^1$, squares of spinors, isomorphic to a qubit | Application of known mathematics [4,7] |
| Physical origin of spin and $SU(2)$ | Not claimed (indication of position only) |

---

## References

1. Bertram Kostant, "Quantization and unitary representations," in *Lectures in Modern Analysis and Applications III*, Lecture Notes in Mathematics 170, 87–208, Springer, 1970. DOI: [10.1007/BFb0079068](https://doi.org/10.1007/BFb0079068).
2. Jean-Marie Souriau, *Structure of Dynamical Systems: A Symplectic View of Physics*, Progress in Mathematics 149, Birkhäuser, 1997 (French original: *Structure des systèmes dynamiques*, Dunod, 1970). DOI: [10.1007/978-1-4612-0281-3](https://doi.org/10.1007/978-1-4612-0281-3).
3. Felix A. Berezin, "General concept of quantization," *Communications in Mathematical Physics*, 40, 153–174, 1975. DOI: [10.1007/BF01609397](https://doi.org/10.1007/BF01609397).
4. Phillip Griffiths and Joseph Harris, *Principles of Algebraic Geometry*, Wiley, 1978 (Wiley Classics Library, 1994). DOI: [10.1002/9781118032527](https://doi.org/10.1002/9781118032527).
5. Shiing-Shen Chern, "Characteristic classes of Hermitian manifolds," *Annals of Mathematics*, 47(1), 85–121, 1946. DOI: [10.2307/1969037](https://doi.org/10.2307/1969037).
6. Elias M. Stein and Guido Weiss, *Introduction to Fourier Analysis on Euclidean Spaces*, Princeton University Press, 1971, Chapter IV. ISBN: 978-0-691-08078-9.
7. Élie Cartan, *The Theory of Spinors*, Hermann, 1966 (French original: *Leçons sur la théorie des spineurs*, 1938; Dover reprint 1981). ISBN: 978-0-486-64070-9.
8. Roger Penrose, "Twistor algebra," *Journal of Mathematical Physics*, 8(2), 345–366, 1967. DOI: [10.1063/1.1705200](https://doi.org/10.1063/1.1705200).
