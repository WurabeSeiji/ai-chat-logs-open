# Shape-Invariant Traveling Waves on a 4-Dimensional Integer Lattice: Existence Conditions, Conservation Structure, and Self-Consistent Loop

**Author**: Noriaki Kihara  
**Affiliation**: WF System Co., Ltd.  
**ORCID**: 0009-0004-6753-4020  
**DOI**: https://doi.org/10.5281/zenodo.19731598  
**Version**: v1.0  
**Date**: April 2026

---

## Abstract

On an $n$-dimensional integer lattice $\mathbb{Z}^n$, we systematically characterize the conditions under which a wave that strictly preserves its amplitude distribution $|\psi|$ under translation — a shape-invariant traveling wave — can consistently exist, along three axes: spherical consistency, multiplicative closure, and mode decomposability. The existence conditions (D1)–(D4) for spherical standing waves are satisfied in perfect-square dimensions $n \in \{4, 9, 16, \ldots\}$, but when one simultaneously requires an associative multiplicative structure (Frobenius' theorem) and explicit mode counting (Jacobi's four-square formula), $n = 4$ is identified as the **unique dimension**. The shape-invariant traveling waves defined on this lattice simultaneously satisfy individuality, coexistence, crossing invariance, composition rules, and symmetry-group classification as labeled objects. Extension to the tangent bundle $\mathbb{Z}^4 \times \mathbb{Z}^4$ embeds a discrete flow, and the self-return loop via orthogonal interaction waves derives a speed bound from Nyquist stability. Furthermore, we show that central projection from a parent space is the unique geometric mechanism that introduces nonlinear effects into a linear system without invoking the concept of time, and we argue that these structures point to **de Sitter spacetime $\mathrm{dS}_4$** as a geometric necessity.

**Keywords**: shape-invariant traveling wave, integer lattice, Hurwitz quaternion order, Jacobi four-square formula, Nyquist stability, central projection, de Sitter spacetime

---

## §1 Introduction and Scope of Claims

This paper describes, as purely geometric, number-theoretic, and algebraic facts, how the dimensional dependence of waves on an $n$-dimensional integer lattice $\mathbb{Z}^n$ that strictly preserve their waveform under translation — shape-invariant traveling waves — is determined.

The claims of this paper are limited to the following four points:

1. The intersection of the existence conditions (D1)–(D4) for spherical standing waves and the structural conditions (F1)–(F5) for shape-invariant traveling waves identifies $n = 4$ as the unique consistent dimension.
2. Shape-invariant traveling waves on $\mathbb{Z}^4$ simultaneously satisfy the five conditions of individuality, coexistence, crossing invariance, composition rules, and symmetry-group classification as labeled objects.
3. Through tangent bundle extension and orthogonal interaction waves, a self-return loop is constructed, and the speed bound is geometrically derived from the Nyquist stability condition.
4. Central projection from a parent space is the unique geometric mechanism that introduces nonlinear effects into a linear system without invoking the concept of time.

The following elements are **outside the scope of this paper**:

- Dynamical stability (eigenvalue analysis for time evolution, convexity of energy functionals)
- Rigorous formulation as solutions of nonlinear wave equations
- Identification with physical particles, reproduction of quantum field theory

The objects called "shape-invariant traveling waves" in this paper are clearly distinguished from "solitons," which carry dynamical implications. "Shape-invariant" refers solely to the geometric property that $|\psi|$ is preserved under translation. "Traveling" refers only to covariance with respect to translation as a lattice automorphism, and does not include dynamical concepts of motion or propagation.

---

## §2 Number-Theoretic Privilege of the 4-Dimensional Integer Lattice

### §2.1 Lagrange's Four-Square Theorem and Completeness of the Squared-Radius Set

Consider the set of squared distances between lattice points in $\mathbb{Z}^n$: $\{|\mathbf{x}|^2 : \mathbf{x} \in \mathbb{Z}^n\}$.

For $n = 2$, values such as $3, 6, 7, \ldots$ are missing, and for $n = 3$, integers of the form $4^k(8m+7)$ are missing. By Lagrange's four-square theorem, only for $n \geq 4$ does

$$\{|\mathbf{x}|^2 : \mathbf{x} \in \mathbb{Z}^n\} = \mathbb{Z}_{\geq 0}$$

hold. $n = 4$ is the **minimum dimension** at which this completeness is achieved.

### §2.2 Jacobi's Four-Square Formula and Non-Emptiness of Spherical Lattice Points

For the number $r_n(N)$ of $\mathbb{Z}^n$ lattice points on the hypersphere of radius $\sqrt{N}$, Jacobi's four-square formula gives

$$r_4(N) = 8 \sum_{\substack{d \mid N \\ 4 \nmid d}} d > 0 \quad (N \geq 1)$$

$r_4(N)$ is positive for all $N \geq 1$ and, moreover, possesses a **closed explicit formula** in terms of divisor-sum functions.

$r_3(N)$ has zeros at $N = 4^k(8m+7)$, and the explicit formulas for $r_5, r_6, r_7$ involve cuspidal parts of automorphic forms and are far more complex. The property of being "positive for all $N$ and expressible in terms of divisor sums" holds only for $r_4$.

### §2.3 Integrality of the Equal-Component Direction

The length of the main diagonal of an $n$-dimensional unit hypercube is $\sqrt{n}$, and the length of a lattice segment along the equal-component direction $(\pm k, \pm k, \ldots, \pm k)$ is $\sqrt{n} \cdot k$. $\sqrt{n}$ is an integer only when $n$ is a perfect square, and $n = 4$ ($\sqrt{4} = 2$) is the **smallest non-trivial such dimension**.

### §2.4 Multiplicative Structure of the Hurwitz Quaternion Order

Within the quaternion algebra $\mathbb{H}$, taking

$$\mathcal{H} = \mathbb{Z}^4 \cup \left(\mathbb{Z}^4 + \left(\tfrac{1}{2}, \tfrac{1}{2}, \tfrac{1}{2}, \tfrac{1}{2}\right)\right)$$

yields a **left Euclidean domain** satisfying norm multiplicativity $|xy| = |x| \cdot |y|$.

By Frobenius' theorem, the only finite-dimensional associative division algebras over the reals are $\mathbb{R}, \mathbb{C}, \mathbb{H}$. The multiplicative integer lattices are limited to $\mathbb{Z}, \mathbb{Z}[i], \mathcal{H}$, with dimensions $1, 2, 4$. $n = 4$ is the **maximum dimension** admitting an associative multiplicative structure.

---

## §3 Existence Conditions for Spherical Standing Waves (D1)–(D4)

We decompose the conditions for standing waves consistent with the family of spheres $\{S(\sqrt{N})\}_{N \geq 0}$ on $\mathbb{Z}^n$ into the following four conditions.

**(D1) Integrality of nodal radii**: The radius sequence of the spherical family forming nodes is arranged at integer intervals along the equal-component direction. Requires $\sqrt{n} \in \mathbb{Z}$.

**(D2) Non-emptiness of nodal spheres**: For any nodal radius, lattice points exist on that sphere. Requires $r_n(N) > 0$ for all $N \geq 1$.

**(D3) Completeness of the nodal-radius-squared set**: The range of values taken by the squared nodal radii covers all of $\mathbb{Z}_{\geq 0}$. Requires $\{|\mathbf{x}|^2 : \mathbf{x} \in \mathbb{Z}^n\} = \mathbb{Z}_{\geq 0}$.

**(D4) Wavenumber duality**: The nodal lattice in real space and wavenumber space share an isomorphic discrete structure (self-duality of $\mathbb{Z}^n$).

**Dimensional dependence**:

| Condition | $n=2$ | $n=3$ | $n=4$ | $n=9$ | $n \geq 5$ general |
|-----------|-------|-------|-------|-------|---------------------|
| (D1) | × ($\sqrt{2}$) | × ($\sqrt{3}$) | ✓ | ✓ | Perfect squares only |
| (D2) | × | × | ✓ | ✓ | $n \geq 4$ |
| (D3) | × | × | ✓ | ✓ | $n \geq 4$ |
| (D4) | ✓ | ✓ | ✓ | ✓ | All $n$ |

**Proposition 3.1**: The set of dimensions simultaneously satisfying (D1)–(D4) is $\{k^2 : k \in \mathbb{Z}_{\geq 2}\} = \{4, 9, 16, 25, \ldots\}$, and the minimum element is $n = 4$.

---

## §4 Definition and Structural Conditions for Shape-Invariant Traveling Waves (F1)–(F5)

### §4.1 Definition of Shape-Invariant Traveling Waves

**Definition 4.1**: A function $\psi : \Lambda \to \mathbb{C}$ on a lattice $\Lambda \subset \mathbb{R}^n$ is a **shape-invariant traveling wave** if there exists a wavevector $\mathbf{k} \in \Lambda^*$ such that for any translation $\mathbf{a} \in \Lambda$,

$$\psi(\mathbf{x} + \mathbf{a}) = e^{i\langle \mathbf{k}, \mathbf{a} \rangle} \cdot \psi(\mathbf{x}), \quad \langle \mathbf{k}, \mathbf{a} \rangle \in \mathbb{Z}$$

holds.

This definition involves neither time, equations, nor energy, and is described solely in terms of a lattice, an inner product, and an exponential map.

### §4.2 Structural Conditions

The geometric conditions for shape-invariant traveling waves to "consistently exist":

**(F1) Phase consistency**: $\langle \mathbf{k}, \mathbf{a} \rangle \in \mathbb{Z}$ for all $\mathbf{a} \in \Lambda$ $\Leftrightarrow$ $\mathbf{k} \in \Lambda^*$

**(F2) Amplitude preservation**: $|\psi(\mathbf{x} + \mathbf{a})| = |\psi(\mathbf{x})|$ (a direct consequence of the definition)

**(F3) Composition closure**: $\psi_1 \cdot \psi_2$ is again a shape-invariant traveling wave $\Leftrightarrow$ $\mathbf{k}_1 + \mathbf{k}_2 \in \Lambda^*$

**(F4) Spherical consistency**: The magnitude of the wavevector $|\mathbf{k}|$ lies on the radius sequence of the spherical family $\Leftrightarrow$ (D1)–(D4)

**(F5) Multiplicative closure**: $\psi(\mathbf{x}) \cdot \psi(\mathbf{y})$ can be written as a wave on the lattice product $\mathbf{x} \cdot \mathbf{y}$ $\Leftrightarrow$ an associative multiplicative structure on $\Lambda$

### §4.3 Dimensional Sieve of the Conditions

- (F1)(F2)(F3): Hold for all $n$ when $\Lambda = \mathbb{Z}^n$, $\Lambda^* = \mathbb{Z}^n$
- (F4): $n \in \{4, 9, 16, 25, \ldots\}$ (perfect-square dimensions)
- (F5): $n \in \{1, 2, 4\}$ ($\mathbb{R}, \mathbb{C}, \mathbb{H}$; the octonions $n = 8$ are excluded due to non-associativity)

**Theorem 4.2**: $(F4) \cap (F5) = \{4\}$.

The integer lattice on which shape-invariant traveling waves simultaneously satisfy spherical consistency and multiplicative closure is **limited to** $\mathbb{Z}^4$.

---

## §5 Conserved Quantities and Invariant Structures (C1)–(C8)

From the definition of a shape-invariant traveling wave $\psi(\mathbf{x} + \mathbf{a}) = e^{i\langle \mathbf{k}, \mathbf{a} \rangle} \cdot \psi(\mathbf{x})$, the following quantities are **strictly conserved** under translation $\mathbf{a} \in \mathbb{Z}^4$.

**(C1) Amplitude distribution**: $|\psi(\mathbf{x} + \mathbf{a})| = |\psi(\mathbf{x})|$

**(C2) Nodal set configuration**: $\{\mathbf{x} : \psi(\mathbf{x}) = 0\}$ is an automorphism under translation by $\mathbf{a}$

**(C3) Phase difference structure**: $\langle \mathbf{k}, \mathbf{x}_i - \mathbf{x}_j \rangle \in \mathbb{Z}$ is preserved for any pair of lattice points

**(C4) Wavevector magnitude**: $|\mathbf{k}|$ is invariant under translation

**(C5) Mode multiplicity**: $r_4(|\mathbf{k}|^2)$ is invariant under translation

In four dimensions, additional invariant structures exist with respect to symmetry groups beyond translation:

**(C6) Conservation under quaternion multiplication**: Under $\mathbf{x} \mapsto \mathbf{x} \cdot u$ by a unit element $u$ ($|u| = 1$) of the Hurwitz order, the spherical average of $|\psi|$ is preserved. This is described as an action of $\mathrm{SO}(4) \cong (\mathrm{SU}(2) \times \mathrm{SU}(2))/\mathbb{Z}_2$.

**(C7) Conservation under $D_4$ triality**: The theta series $\theta_{D_4}(\tau)$ is invariant under cyclic permutation of the vector, spinor, and conjugate spinor representations.

**(C8) Conservation under modular transformations**: $\theta_{\mathbb{Z}^4}(\tau)$ transforms as a modular form of weight 2 under $\Gamma_0(4)$.

**Proposition 5.1**: (C1)–(C5) hold in all dimensions, but (C6)–(C8) hold only for $n = 4$. The shape-invariant traveling wave in four dimensions is the **unique structure** that simultaneously preserves invariants under translation, quaternion multiplication, triple duality, and modular transformations.

---

## §6 Coexistence of Multiple Waves and Crossing Invariance (G1)–(G5)

### §6.1 Superposition and Complete Separation

Consider multiple shape-invariant traveling waves $\psi_{\mathbf{k}_j}(\mathbf{x}) = e^{i\langle \mathbf{k}_j, \mathbf{x} \rangle}$ on $\mathbb{Z}^4$.

The operation of extracting each component from the superposition $\Psi(\mathbf{x}) = \sum_j c_j \cdot \psi_{\mathbf{k}_j}(\mathbf{x})$ can be performed exactly by the uniqueness of the discrete Fourier transform on $\mathbb{Z}^4$:

$$c_j = \lim_{N \to \infty} \frac{1}{(2N+1)^4} \sum_{\mathbf{x} \in [-N,N]^4} \Psi(\mathbf{x}) \cdot e^{-i\langle \mathbf{k}_j, \mathbf{x} \rangle}$$

Since the dual basis of $\mathbb{Z}^4$ forms a complete orthogonal system, any finite or countable number of labeled waves are **completely separated**.

### §6.2 Geometric Realization of Particle-Like Conditions

Regarding a shape-invariant traveling wave $\psi_{\mathbf{k}}(\mathbf{x})$ as a "geometric object identified by label $\mathbf{k}$," the following hold in four dimensions:

**(G1) Individuality**: Each $\mathbf{k} \in \mathbb{Z}^4$ provides a unique label, and objects are identified by a countable discrete set.

**(G2) Coexistence**: Any number of labeled objects consistently coexist in the same space (closure under superposition).

**(G3) Crossing invariance**: Even when multiple objects occupy the same region, each label is **completely recovered** by Fourier decomposition.

**(G4) Composition rule**: Composition is defined by label addition $\mathbf{k}_1 + \mathbf{k}_2 \in \mathbb{Z}^4$, and is closed as multiplication in the Hurwitz order.

**(G5) Classification system**: Labels are organized group-theoretically under $\mathrm{SO}(4)$, $D_4$ triality, and modular transformations.

**Theorem 6.1**: The unique dimension in which shape-invariant traveling waves can exist on an integer lattice while simultaneously satisfying (G1)–(G5) is $n = 4$.

**Remark**: (G1)–(G5) embed structural correspondences to the geometric conditions required of particles, but this paper claims only the fact that "the geometric structure satisfies the conditions" and does not make the identification "therefore they are particles."

---

## §7 Label Structure (L1)–(L5)

The wavevector $\mathbf{k} \in \mathbb{Z}^4$ of a shape-invariant traveling wave functions as an integer-valued label. In four dimensions, this label set simultaneously possesses the following five layers of structure.

**(L1) Group structure**: The label set $\mathbb{Z}^4$ is closed as a quaternion additive group and further possesses the multiplicative group structure of the Hurwitz order $\mathcal{H}$. Addition and multiplication between labels are group-theoretically closed.

**(L2) Geometric hierarchy**: $|\mathbf{k}|^2$ takes integer values, and the label space is naturally stratified into spherical shells by the Jacobi four-square formula $r_4(|\mathbf{k}|^2)$. The number of labels on each shell $S(\sqrt{N})$ is explicitly given by a divisor sum.

**(L3) Symmetry group action**: $\mathrm{SO}(4) \cong (\mathrm{SU}(2) \times \mathrm{SU}(2))/\mathbb{Z}_2$ acts on the label set, and labels are organized as representations of the 4-dimensional rotation group.

**(L4) Triality**: $D_4$ triality classifies labels into three types of representations (vector, spinor, and conjugate spinor), which are cyclically permutable.

**(L5) Modular transformations**: The theta series $\theta_{\mathbb{Z}^4}(\tau)$ transforms modularly under $\Gamma_0(4)$, establishing scale duality of the label distribution.

**Proposition 7.1**: The label structure simultaneously satisfying all of (L1)–(L5) is limited to $\mathbb{Z}^4$. (L1) holds for $n \in \{1, 2, 4\}$, (L2) admits a concise closed form only for $n = 4$, and (L4) holds only for $n = 4$.

---

## §8 Tangent Bundle Extension and Discrete Flow

### §8.1 Tangent Bundle Extension of the Label Space

We extend the label $\mathbf{k} \in \mathbb{Z}^4$ of a shape-invariant traveling wave to the pair of label and rate of change $(\mathbf{k}, \Delta\mathbf{k}) \in \mathbb{Z}^4 \times \mathbb{Z}^4$. This corresponds to the **discrete tangent bundle** $T\mathbb{Z}^4 = \mathbb{Z}^4 \times \mathbb{Z}^4$ of the label space $\mathbb{Z}^4$.

The extended wave is written as

$$\psi_{(\mathbf{k}, \Delta\mathbf{k})}(\mathbf{x}, t) = e^{i\langle \mathbf{k} + t \cdot \Delta\mathbf{k},\; \mathbf{x} \rangle}$$

Here $t \in \mathbb{Z}$ is a parameter of the discrete flow and is not physical time.

### §8.2 Discrete Flow

The self-update map $U : (\mathbf{k}, \Delta\mathbf{k}) \mapsto (\mathbf{k} + \Delta\mathbf{k}, \Delta\mathbf{k})$ generates a discrete flow on $\mathbb{Z}^4 \times \mathbb{Z}^4$. The wave itself embeds the rule $\Delta\mathbf{k}$ for how to advance its label, and the external concept of time appears only as the flow parameter $t$.

### §8.3 Structures Specific to Four Dimensions

The extended label space $\mathbb{Z}^4 \times \mathbb{Z}^4$ possesses the double structure $\mathcal{H} \times \mathcal{H}$ of the Hurwitz order. The self-update flow is organized as an $\mathrm{SO}(4)$ action, and the quaternion product operations on $(\mathbf{k}, \Delta\mathbf{k})$ remain closed within $\mathbb{Z}^4 \times \mathbb{Z}^4$.

---

## §9 Self-Return Loop and Nyquist Stability

### §9.1 Emission of Orthogonal Interaction Waves

From the central phase of a wave $\psi_{(\mathbf{k}, \Delta\mathbf{k})}$, an interaction wave $\chi$ is emitted in directions orthogonal to the wavevector $\mathbf{k}$ in four dimensions. By the spherical wave regularity of (D1)–(D4), $\chi$ is a 4-dimensional spherical wave that spreads in a shape-invariant manner on $\mathbb{Z}^4$.

### §9.2 Self-Recrossing Condition

For the trajectory of the wave $\tau \mapsto \mathbf{x}_{\mathrm{self}}(\tau)$ and the current position $\mathbf{x}_{\mathrm{self}}(t)$, the recrossing condition

$$|\mathbf{x}_{\mathrm{self}}(t) - \mathbf{x}_{\mathrm{self}}(\tau)| = c(t - \tau)$$

is satisfied for some $\tau < t$, meaning that a $\chi_\tau$ emitted in the past acts multiplicatively on the current $\psi$. If the wave is stationary ($\Delta\mathbf{k} = 0$), $\chi$ merely spreads outward and does not return to itself. Self-recrossing occurs only when the wave is in motion ($\Delta\mathbf{k} \neq 0$), causing a label transformation $(\mathbf{k}, \Delta\mathbf{k}) \to (\mathbf{k}', \Delta\mathbf{k}')$.

### §9.3 Nyquist Stability and Speed Bound

We analyze this self-return loop from a control-theoretic perspective. The phase structure of the open-loop transfer function:

- **Fundamental phase**: $-\pi/2$ (emission in the orthogonal direction)
- **Delay phase**: $-\omega \Delta t$ (recrossing delay)
- **Composite phase**: $-\pi/2 - \omega \Delta t$

The amplitude decay of 4-dimensional spherical waves is $1/r^{3/2}$, and the loop gain decays sufficiently fast with distance. The Nyquist trajectory does not encircle the $(-1, 0)$ point, and the loop is stable.

**Derivation of the speed bound**: From the effective delay $\Delta t \propto 1/(c - v_\parallel)$, as the wave's travel speed $v$ approaches the propagation speed $c$, we have $\Delta t \to 0$, and the divergence boundary frequency $\omega_c = \pi/(2\Delta t) \to \infty$. At $v = c$, all frequencies reach the divergence boundary.

**Proposition 9.1**: As a condition for maintaining Nyquist stability, $v \leq c$ is geometrically derived. The speed bound and displacement bound are unified as different approach paths to the same stability boundary $\angle G = -\pi, |G| = 1$.

---

## §10 Central Projection and the Geometric Origin of Nonlinear Effects

### §10.1 Limitation of Linear Systems

Within a single subjective space $\mathbb{Z}^4$, shape-invariant traveling waves **coexist without interaction** due to the closure property of linear superposition. Addition operations alone cannot transform the label set (by the uniqueness of Fourier decomposition).

### §10.2 Four Methods of Multiplicative Action

Methods for introducing nonlinear effects into a linear system:

**(i) Duality in phase space**: Formulate the system as a $U(1)$-valued function and read phase addition as amplitude multiplication.

**(ii) Transformation of the evaluation-point lattice**: Apply quaternion rotation $\mathbf{x} \mapsto \mathbf{x} \cdot u$ to the evaluation-point space, inducing a multiplicative action $\mathbf{k} \mapsto \mathbf{k} \cdot u^*$ on the label space. Specific to $n = 4$.

**(iii) Phase composition of two systems**: Phase composition through same-point evaluation of two independent systems.

**(iv) Central projection from a parent space**: Project a linear wave on the parent space $S^5(R)$ onto the child space via central projection $\pi : S^5 \to \mathbb{Z}^4$. Because the projection map is a fractional (nonlinear) map, the result appears on the child space as a function with **nonlinear phase**.

### §10.3 Privilege of Central Projection

Method (iv) is qualitatively different from the other three. Methods (i)–(iii) exploit the internal structure of the child space $\mathbb{Z}^4$ and produce changes within the range of the child space's symmetries. Method (iv) is an action from **outside the frame**:

- Operations in the parent space are mere addition (linear)
- On the child space, they appear as nonlinear changes of unknown origin
- The symmetries of the child space themselves can be deformed "from the outside"

**Proposition 10.1**: The unique method that introduces nonlinear effects into a linear system in a geometrically justified manner without introducing dynamics is central projection from a higher-dimensional parent space.

---

## §11 Convergence to de Sitter Spacetime

### §11.1 Soliton Stability on Closed Spheres

As shown in the companion paper [1], shape-invariant standing waves can exist stably on a closed $n$-dimensional sphere $S^n$ **only in odd dimensions**. On even-dimensional spheres, the hairy ball theorem (the nonexistence of a nowhere-vanishing vector field) causes the standing wave structure to break down.

### §11.2 Intersection of Two Arguments

| Argument | Target space | Privileged dimension | Basis |
|----------|-------------|---------------------|-------|
| This paper (discrete side) | Infinite flat lattice $\mathbb{Z}^n$ | $n = 4$ | (F4)∩(F5)={4} |
| Companion paper (continuous side) | Closed sphere $S^n$ | Odd dimensions | Huygens' principle + hairy ball theorem |

Both arguments hold independently, but they impose constraints on the construction of physical spacetime:

- **Spatial part**: Solitons are stable on a closed sphere → **odd dimension** → minimum is $S^3$
- **Temporal part**: Shape-invariant waves propagate along an infinitely open direction → **hyperbolic 1-dimension**
- **Overall**: (D1)–(D4) hold in a locally flat 4-dimensional setting → $3 + 1$ dimensions

### §11.3 Geometric Necessity of de Sitter Spacetime

De Sitter spacetime $\mathrm{dS}_4$ is defined as the hyperboloid in 5-dimensional Minkowski space $\mathbb{M}^{1,4}$:

$$-X_0^2 + X_1^2 + X_2^2 + X_3^2 + X_4^2 = R^2$$

The spatial part is a closed $S^3$ (odd-dimensional, soliton-stable), and the temporal part is hyperbolic (infinitely open, the arena for shape-invariant wave propagation).

**Proposition 11.1**: The only 4-dimensional spacetime structure that simultaneously satisfies the existence conditions for shape-invariant traveling waves (discrete side: $n = 4$ is unique) and the existence conditions for closed solitons (continuous side: odd dimensions only) is de Sitter spacetime $\mathrm{dS}_4$, described as the product structure of the space $S^3$ and hyperbolic time.

### §11.4 Roles of Parent and Child Spaces

Within the framework of the central projection series:

- **Child space** $S^4(R)$ (even-dimensional, interior $\mathbb{Z}^4$): The arena for shape-invariant traveling waves. Closed solitons are unstable on the surface.
- **Parent space** $S^5(R)$ (odd-dimensional): The source of closed solitons. Introduces nonlinear effects into the child space via central projection.

The structure in which the parent space is odd-dimensional and the spatial part of the child space is odd-dimensional is **necessarily** determined by the geometric difference between odd and even dimensions.

---

## §12 Summary

We organize the logical flow presented in this paper:

1. The existence conditions (D1)–(D4) for spherical standing waves on $\mathbb{Z}^n$ are satisfied in perfect-square dimensions $\{4, 9, 16, \ldots\}$ (§3).
2. Adding the structural conditions (F1)–(F5) for shape-invariant traveling waves, $(F4) \cap (F5) = \{4\}$ identifies $n = 4$ as the unique dimension (§4).
3. Shape-invariant traveling waves on $\mathbb{Z}^4$ possess conserved quantities (C1)–(C8), geometrically realize particle-like conditions (G1)–(G5), and are organized by label structures (L1)–(L5) (§5–§7).
4. Through tangent bundle extension $\mathbb{Z}^4 \times \mathbb{Z}^4$, a discrete flow is defined, and the self-return loop of orthogonal interaction waves derives the speed bound $v \leq c$ from Nyquist stability (§8–§9).
5. Central projection from a parent space is the unique geometric mechanism that introduces nonlinear effects into a linear system without invoking the concept of time (§10).
6. Two independent arguments from the discrete side ($n = 4$) and the continuous side (odd-dimensional spheres) point to de Sitter spacetime $\mathrm{dS}_4$ as a geometric necessity (§11).

All claims in this paper are stated as purely geometric, number-theoretic, and algebraic facts, and require no physical interpretation. Physical correspondences (particles, interactions, spacetime structure) remain open as **interpretations connectable** to these geometric facts, but do not affect the validity of the propositions themselves.

---

## References

[1] Noriaki Kihara, "Shape-Invariant Standing Waves on Closed Spheres: Dimensional Stability and the Geometric Necessity of 3+1 Spacetime," 2026.

[2] Noriaki Kihara, "A Three-Layer Model of the Universe via Central Projection," Zenodo, 2025. DOI: 10.5281/zenodo.15083729.

[3] Noriaki Kihara, "Regularity of Inter-Sphere Projection (Gnomonic Orthoprojection)," Zenodo, 2025. DOI: 10.5281/zenodo.15292757.

[4] Noriaki Kihara, "(R,Q) Mapping Definition Paper [M3]," Zenodo, 2025. DOI: 10.5281/zenodo.19692853.
