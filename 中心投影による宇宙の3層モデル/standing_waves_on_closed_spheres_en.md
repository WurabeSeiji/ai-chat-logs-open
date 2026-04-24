# Shape-Invariant Standing Waves on Closed Spheres — Dimensional Stability and the Geometric Necessity of 3+1 Spacetime

**Author**: Noriaki Kihara  
**Affiliation**: WF System Co., Ltd.  
**ORCID**: 0009-0004-6753-4020  
**DOI**: https://doi.org/10.5281/zenodo.19731594  
**Version**: v1.0  
**Date**: April 2026

---

## Abstract

We systematically analyze the phase wave equation on closed $n$-dimensional spheres $S^n$ from one to four dimensions. We show that the phase wave equation is independent of the spatial curvature $R$, and that $R$ enters only through the standing wave condition (the global boundary condition). We demonstrate that the parameter description $(\theta_0, \lambda_{\text{base}}, A, N)$ of shape-invariant standing waves is preserved across all dimensions, that the conserved quantity $\int |u|^2 \, dV$ carries the identity of the wave packet, and that Huygens' principle holds only on odd-dimensional spheres and fails on even-dimensional ones. From these facts, we show that **the minimal nontrivial sphere on which shape-invariant standing waves can stably exist is $S^3$ (three-dimensional)**. Furthermore, we add an independent argument from the discrete-algebraic side via Frobenius' theorem and argue that both the continuous and discrete sides necessitate a **3+1 spacetime structure**.

**Keywords**: phase wave equation, shape-invariant standing wave, closed sphere, Huygens' principle, Frobenius' theorem, 3+1 spacetime

---

## §1 Introduction and Scope of Claims

This paper describes, as purely geometric and number-theoretic facts, the dimensional dependence of standing waves of the phase wave equation on closed $n$-dimensional spheres $S^n$.

The claims of this paper are limited to the following three points:

1. The phase wave equation is an operation essentially orthogonal to the spatial curvature $R$, and $R$ enters only through the standing wave condition (the commensurability of the circumference of the closed space with the wavelength).
2. Shape-invariant standing waves that stably preserve their identity as wave packets exist only on odd-dimensional spheres where Huygens' principle holds. The minimal nontrivial dimension is $S^3$.
3. By Frobenius' theorem (associative division algebras are limited to $\mathbb{R}, \mathbb{C}, \mathbb{H}$), four dimensions are independently necessitated from the discrete-algebraic side, and the $1+3$ decomposition of the quaternions algebraically embeds the $3+1$ spacetime structure.

The following elements are **outside the scope of this paper**:

- Dynamical stability (eigenvalue analysis under time evolution, convexity of the energy functional)
- Rigorous formulation as soliton solutions of nonlinear wave equations
- Concrete correspondences with physical spacetime

The object called "shape-invariant standing wave" in this paper is distinguished from "soliton," which carries dynamical implications. "Shape-invariant" refers solely to the geometric property that the amplitude distribution $|\psi|$ is preserved under parallel transport, and "standing" refers solely to the property that the nodal structure is compatible with the arrangement of the family of spheres.

---

## §2 Shape-Invariant Standing Waves on a Closed String ($S^1$)

### §2.1 Parameter Description

On a closed string ($S^1$) of circumference $L$, a shape-invariant standing wave is completely described by the following four parameters:

- $\theta_0$: positional phase of the maximum wave height
- $\lambda_{\text{base}}$: fundamental mode wavelength
- $A$: amplitude of the maximum wave height
- $N$: number of modes in the Fourier expansion (truncation order)

Here $N$ is the parameter controlling the "roundness/sharpness" of the waveform:

- $N = 1$: fundamental mode only (pure sinusoidal wave)
- $N$ intermediate: trapezoidal waveform with added odd-order harmonics
- $N \to \infty$: perfect square wave (fully localized pulse)

### §2.2 Standing Wave Condition

The necessary and sufficient condition for a standing wave to exist is:

$$L = k \cdot \lambda_{\text{base}} \quad (k \in \mathbb{Z}_{>0})$$

That is, the circumference of the closed string must be a positive integer multiple of the fundamental mode wavelength. Under this condition, the permitted wavelengths are automatically discretized to $\lambda_n = L/n$ ($n = 1, 2, 3, \ldots$).

### §2.3 Unification of Boundary Conditions

The phase wave equation requires the string to be closed but does not specify the manner of closure. The difference between free-end reflection (phase shift $0$) and fixed-end reflection (phase shift $\pi$) can be treated uniformly, via the method of images, as the even-parity/odd-parity subspaces of a closed string of length $2L$. The difference in boundary conditions is absorbed into the phase selection ($0$ or $\pi$) of the fundamental mode.

---

## §3 Orthogonality of the Phase Wave Equation and Curvature $R$

### §3.1 $R$-Independence

The phase wave equation is a local differential equation, and the spatial curvature $R$ (or equivalently $L = 2\pi R$) appears nowhere in it. $R$ enters only at **the level of the global boundary condition**, which selects which wavelengths can constitute standing waves:

- $\lambda = L/n$ ($n$ integer): phase shifts by $2\pi n$ over one circuit, and a standing wave is established
- $\lambda > L$: phase shifts by less than $2\pi$ over one circuit, and no standing wave forms
- $\lambda$ incommensurate with $L$: phase never returns to its original value regardless of the number of circuits, and the trajectory fills the phase space quasi-periodically

In every case the wave equation itself is identical; whether $R$ is small or large, the wave obeys the same equation.

### §3.2 A Posteriori Definition of $R$

From this structure, $R$ is not a quantity within the wave equation but rather **a geometric quantity that can be read off a posteriori by surveying the set of standing-wave solutions**:

$$R = \frac{\lambda_{\text{base}}}{2\pi}$$

That is, $R$ is "$(2\pi)^{-1}$ times the longest wavelength admitted as a standing wave" — a derived quantity, not an a priori parameter.

### §3.3 Structural Orthogonality

**This is not a Euclidean approximation.** The phase wave equation and the spatial curvature are essentially orthogonal operations, and calculations in Euclidean coordinates are exactly correct. Curvature $R$ enters only through the global boundary condition (the commensurability of the circumference of the closed space with the wavelength) and has no effect whatsoever on the local solutions of the wave equation.

---

## §4 Extension to the Sphere $S^2$ and Antipodal Convergence

### §4.1 Preservation of Parameter Structure

Shape-invariant standing waves on the sphere $S^2$ have essentially the same parameter structure as on $S^1$. By the $\mathrm{SO}(3)$ symmetry of $S^2$, the standing wave condition takes the same form regardless of which great circle is chosen:

$$2\pi R = n \cdot \lambda_{\text{base}}$$

The parameters change only in the dimension of the positional phase and are exhausted by the four quantities $(\theta_0, \lambda_{\text{base}}, A, N)$. The position of the wave packet is a completely free parameter and is indistinguishable from internal observation.

### §4.2 Antipodal Convergence Phenomenon

On $S^2$, when the great-circle circumference is not an integer multiple of $\lambda$, a wave packet localized at a single point executes the following motion:

1. Initially: a peak localized at a single point
2. Expands as a ring-shaped ($S^1$) wavefront
3. Ring reaches maximum near the equator
4. Contracts toward the antipodal point on the opposite side
5. Converges again to a single point at the antipodal point
6. Repeats

### §4.3 Geometric Dispersion

The eigenfrequencies of $S^2$ are $\omega_l = c\sqrt{l(l+1)}/R$, which are **non-equally spaced**. As a result, the wave packet is strictly quasi-periodic rather than periodic, and the shape of the wave packet deteriorates with each successive cycle.

---

## §5 Conserved Quantity $\int |u|^2 \, dV$

### §5.1 Distinction Between Instantaneous Wave Height and Conserved Quantity

In a perfectly lossless linear wave equation, the instantaneous wave height $u(x,t)$ itself is not a "stable quantity" either locally or globally. What is conserved is the following nonlocal integral quantity:

$$\int_{S^n} |u|^2 \, dV = \sum_{l,m} |a_{lm}|^2 = \text{const}$$

The "identity" of a wave packet resides not in its instantaneous shape but in this squared integral over the entire space.

### §5.2 Structural Isomorphism with the Schrödinger Equation

This conservation law is structurally isomorphic to $\int |\psi|^2 \, dV = 1$ in the Schrödinger equation of quantum mechanics:

| | Shape-invariant standing wave on a closed sphere | Schrödinger equation |
|---|---|---|
| Field | $u(x,t)$ (real-valued wave height) | $\psi(x,t)$ (complex amplitude) |
| Conserved density | $u^2$ | $|\psi|^2$ (probability density) |
| Conserved quantity | $\int u^2 \, dA$ | $\int |\psi|^2 \, dV = 1$ |
| Expansion basis | Spherical harmonics $Y_{lm}$ | Hamiltonian eigenfunctions |

This isomorphism is not coincidental but a property arising from the common mathematical skeleton of "a linear equation expandable in an orthogonal basis"; the distinction between classical and quantum is not essential.

---

## §6 Validity of Huygens' Principle on $S^3$

### §6.1 Qualitative Difference Between Odd and Even Dimensions

The propagation properties of the wave equation differ essentially between even and odd spatial dimensions (Hadamard's result):

- **Odd dimensions** ($S^1, S^3, S^5, \ldots$): strict Huygens' principle holds. A point-source pulse propagates solely as a sharp wavefront, leaving no trailing wake.
- **Even dimensions** ($S^2, S^4, \ldots$): Huygens' principle fails. A "wake" persists after the wavefront has passed.

### §6.2 Shape-Invariant Standing Waves on $S^3$

On $S^3$ (the hyperspherical shell of a perfectly symmetric hypersphere in four-dimensional space):

- Conserved quantity: $\int |u|^2 \, dV$ ($V = 2\pi^2 R^3$)
- Standing wave condition: $2\pi R = n \cdot \lambda_{\text{base}}$ ($\mathrm{SO}(4)$ symmetry)
- Parameters: $(\theta_0, \phi_0, \psi_0, \lambda_{\text{base}}, A, N)$

By Huygens' principle, a wave packet localized at a single point **propagates solely as a sharp spherical shell ($S^2$-shaped wavefront)**, leaving nothing behind. It converges exactly to a single point at the antipodal point, where the peak is restored, and then expands again as a spherical shell. Compared to the quasi-periodic motion with afterglow on $S^2$, a **cleaner restoration** is achieved on $S^3$.

### §6.3 Failure of Huygens' Principle on $S^4$

On $S^4$, Huygens' principle fails again, and afterglow persists behind the wavefront, as on $S^2$. The clean restoration achieved on $S^3$ is lost.

### §6.4 Alternating Pattern of Dimension and Stability

| Dimension | Huygens' principle | Behavior of wave packet |
|------|----------------|----------------|
| $S^1$ | (one-dimensional) exact | Phase rotation only, perfectly periodic |
| $S^2$ | Fails | Antipodal convergence + afterglow, quasi-periodic |
| $S^3$ | **Holds** | Sharp spherical-shell propagation, clean restoration |
| $S^4$ | Fails | Higher-dimensional analogue of $S^2$, afterglow present |
| $S^5$ | **Holds** | Higher-dimensional analogue of $S^3$, clean |

---

## §7 Conditions for Stable Existence of Shape-Invariant Standing Waves

### §7.1 Loss of Wave Packet Identity in Even Dimensions

We organize the phenomena occurring in even dimensions ($S^2, S^4$):

- The conserved quantity $\int |u|^2 \, dV$ is exactly conserved
- Each individual mode amplitude $|a_{lm}|$ is also exactly conserved
- However, **the identity as a localized wave packet** is not maintained

Because the eigenfrequencies $\omega_l$ are arranged in irrational ratios, the phases of the individual modes drift independently, undergoing quasi-periodic motion in phase space. Under long-time averaging, a distribution spread over the entire sphere becomes the typical state, with localization partially reviving by occasional coincidental phase alignment (Poincaré recurrence). This is consistent with the picture of a quasi-periodic dynamical system in the sense of KAM theory.

### §7.2 Necessary and Sufficient Conditions for Stable Shape-Invariant Standing Waves

The conditions for a shape-invariant standing wave to "stably exist while preserving its identity as a wave packet" are:

1. **Huygens' principle holds**: the wavefront leaves no trailing wake
2. **Eigenfrequencies are close to integer ratios**: phase coherence is regenerated

These are simultaneously satisfied **only on odd-dimensional spheres**.

### §7.3 Minimality

$S^1$ trivially satisfies Huygens' principle, but in one dimension it exhibits only the degenerate behavior of phase rotation. The minimal dimension that supports nontrivial wave packet formation, restoration, and antipodal convergence is **$S^3$**.

**Proposition**: The minimal dimension $n$ on which shape-invariant standing waves can exist nontrivially and stably on a closed $n$-dimensional sphere is $n = 3$.

---

## §8 Independent Argument via Frobenius' Theorem

### §8.1 Associative Division Algebras

By Frobenius' theorem (1877), the finite-dimensional associative division algebras over the real numbers are limited to the following three:

- $\mathbb{R}$ (real numbers, 1-dimensional)
- $\mathbb{C}$ (complex numbers, 2-dimensional)
- $\mathbb{H}$ (quaternions, 4-dimensional)

No natural division algebra exists in three dimensions. Extending while preserving algebraic structure (addition, subtraction, multiplication, division, and norm) skips three dimensions and arrives at four.

### §8.2 The $1+3$ Decomposition of Quaternions

A quaternion $q = a + bi + cj + dk$ decomposes intrinsically within the algebraic structure into $1$ dimension (scalar part $a$) $+ 3$ dimensions (vector part $bi + cj + dk$). This decomposition is not arbitrary but a natural structure demanded by the conjugation, norm, and multiplication of quaternions.

### §8.3 Upper Bound: Non-Associativity of Octonions

Extending to the 8-dimensional octonions $\mathbb{O}$ breaks the associative law (Hurwitz's theorem). Since losing associativity makes it difficult to describe the composition and interaction of waves, four dimensions constitute the maximal dimension for an associative division algebra.

---

## §9 Convergence of the Dual Argument

### §9.1 Independent Arguments from the Continuous and Discrete Sides

| Aspect | Conclusion | Basis |
|------|------|------|
| Continuous waves (phase wave equation on $S^n$) | Three-dimensional space is preferred | Huygens' principle in odd dimensions |
| Discrete algebra (lattice structure) | Four dimensions are necessitated | Frobenius' theorem |
| Synthesis | **3 + 1 = 4** | $1+3$ decomposition of quaternions |

### §9.2 Unified Proposition

The "shape-invariant standing waves stable in odd dimensions" on the continuous side and "Frobenius' theorem" on the discrete side, though seemingly independent, are two facets of a single requirement: **the minimal spacetime in which stable shape-invariant standing waves can exist**:

1. **Spatial structure**: a space that geometrically protects the identity of wave packets → odd-dimensional spheres → the minimal nontrivial dimension is $S^3$ (3-dimensional)
2. **Algebraic structure**: interactions of wave packets close within the same category → associative division algebras → the minimal nontrivial dimension is $\mathbb{H}$ (4-dimensional)
3. **Synthesis**: via the $1+3$ decomposition of quaternions, unified as 3-dimensional space $+$ 1-dimensional time axis

**Theorem (Convergence of the dual argument)**: The minimal closed space on which shape-invariant standing waves stably exist and possess algebraically closed interactions is $S^3$, and the minimal structure describing this within an associative division algebra is the $1+3$ decomposition of the quaternions $\mathbb{H}$, namely the $3+1$ spacetime structure.

---

## §10 Orthogonality of the Phase Wave Equation and Spatial Curvature

The orthogonality between the phase wave equation and the curvature $R$, demonstrated in §3, is preserved in higher dimensions. For $S^1 \to S^2 \to S^3 \to S^4$ alike:

- The wave equation is a universal equation independent of $R$
- $R$ is a geometric quantity determined a posteriori through the standing wave condition
- Calculations in Euclidean coordinates are exactly correct (not an approximation)

This orthogonality means that shape-invariant standing waves are indifferent to spatial distortion. Waves propagate along geodesics, but whether the space is distorted or not has no effect on the local solutions of the wave equation.

---

## §11 Summary

We summarize the main results of this paper:

1. **$R$-independence**: The phase wave equation and the spatial curvature are essentially orthogonal, and $R$ is a geometric quantity that can be read off a posteriori from the set of standing-wave solutions.

2. **Universality of the parameter description**: The parameters $(\theta_0, \lambda_{\text{base}}, A, N)$ of shape-invariant standing waves are described by the same skeleton from $S^1$ to $S^n$, with only the dimension of the positional phase varying.

3. **Conserved quantity**: The identity of a wave packet is carried not by the instantaneous wave height $u$ but by the nonlocal squared integral $\int |u|^2 \, dV$.

4. **Dimensional selection via Huygens' principle**: Shape-invariant standing waves can stably exist only on odd-dimensional spheres, and the minimal nontrivial dimension is $S^3$ (three-dimensional).

5. **Convergence of the dual argument**: The continuous side (Huygens' principle) and the discrete side (Frobenius' theorem) independently necessitate the same $3+1$ spacetime structure.

---

## References

- [P1] Noriaki Kihara, "Central Projection Between Hyperspheres," 2025.
- [P2] Noriaki Kihara, "Regularity of Central Projection Between Hyperspherical Shells," DOI: 10.5281/zenodo.19692192, 2026.
- [W8] Noriaki Kihara, "Six-Dimensional Hyperrectangular Parallelepiped and Particle Model," 2025.
- Hadamard, J., *Lectures on Cauchy's Problem in Linear Partial Differential Equations*, 1923.
- Frobenius, G., "Über lineare Substitutionen und bilineare Formen", *J. reine angew. Math.*, **84**, 1–63, 1877.
- Hurwitz, A., "Über die Composition der quadratischen Formen von beliebig vielen Variablen", *Nachr. Ges. Wiss. Göttingen*, 309–316, 1898.
