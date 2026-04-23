# Four Modes of Shape-Invariant Waves: Spin, Chirality, and Statistics Determined by Wave Vector Structure

**Author**: Noriaki Kihara  
**Affiliation**: WF System Co., Ltd.  
**ORCID**: 0009-0004-6753-4020  
**Version**: v1.0  
**Date**: April 2026  
**DOI**: 10.5281/zenodo.19709798

---

## Abstract

Building on the shape-invariant standing waves on the closed sphere $S^3$ introduced in the companion paper [1] and the shape-invariant traveling waves on $\mathbb{Z}^4$ formulated in [2], we classify shape-invariant waves based on the number of nonzero components $\kappa$ of the wave vector. Waves with $\kappa = 2$ (two nonzero components) constitute direction-constrained waves whose displacement direction is restricted to one spatial dimension, and the signed area of their two-axis plane admits only two states, thereby providing a geometric origin for the Pauli exclusion principle. Waves with $\kappa = 0$ (all spatial components zero) propagate spherically; when the temporal component $k_t = 0$, they are time-reversal symmetric and thus carry no intrinsic causal ordering. Waves with $\kappa = 1$ (only one spatial component nonzero, temporal component zero) are standing waves that possess scalar character due to the equivalence of three directions and fill the entire subjective space as a ground mode. Direction-constrained waves with $\kappa = 2$ possessing a (space, $t$) plane have chirality defined by the orientation of the $t$-axis, whereas waves possessing only (space, space) planes have degenerate chirality and may satisfy the Majorana condition. Central projection onto the child space induces geodesic curvature, and acceleration is expressed as a discrete update rule for the wave vector.

**Keywords**: shape-invariant wave, mode classification, direction-constrained wave, signed area, chirality, Pauli exclusion principle, Majorana condition, central projection, acceleration

---

## §1 Introduction and Scope of Claims

This paper demonstrates, within the framework of shape-invariant waves constructed in the companion papers [1] and [2], that the structure of the wave vector determines four qualitatively distinct modes, and that the geometric properties of each mode provide the origins of spin, chirality, and statistics.

The claims of this paper are limited to the following five points:

1. The number of nonzero components $\kappa$ of the wave vector classifies shape-invariant waves into three types: $\kappa \in \{0, 1, 2\}$. No stable modes exist for $\kappa \geq 3$.
2. Waves with $\kappa = 2$ admit only two states of signed area orientation in their two-axis plane, geometrically prohibiting the coexistence of three or more waves with the same wave vector.
3. For waves with $\kappa = 2$, the (space, $t$) plane defines chirality from the fixed direction of the $t$-axis, whereas chirality is degenerate in the (space, space) plane.
4. Spherical waves with $\kappa = 0$ and $k_t = 0$ carry no causal ordering, and the shortest geodesic path is retroactively constructed from the interaction outcome.
5. Central projection onto the child space induces geodesic curvature on wave trajectories and geometrically expresses acceleration.

The following elements are **outside the scope of this paper**:

- Identification with specific physical particles
- Numerical derivation of coupling constants
- Reproduction in the continuous limit of gauge groups
- Numerical derivation of mass values and mixing matrices

**Note on coordinate systems**: This paper primarily works on $\mathbb{Z}^4$ (four axes $xyzt$). However, to avoid the problem that spherical waves with $\kappa_s = 0$ become trivial (constant functions) within $\mathbb{Z}^4$, the extended axes introduced in [4] ($R$: radial direction of central projection, $Q$: discrete label axis) are invoked where necessary. See [4] for details on these extended axes.

---

## §2 3+1 Decomposition of the Wave Vector and Mode Classification

### §2.1 3+1 Decomposition

As shown in §11 of [2], the existence condition for shape-invariant traveling waves and the stability condition for standing waves on a closed sphere point to the structure of de Sitter spacetime $\mathrm{dS}_4$. Within this structure, the four axes of $\mathbb{Z}^4$ naturally decompose into 3+1:

- **Spatial components**: $\mathbf{k}_s = (k_1, k_2, k_3) \in \mathbb{Z}^3$ (corresponding to the three directions of $S^3$)
- **Temporal component**: $k_t \in \mathbb{Z}$ (corresponding to the hyperbolic time direction)

### §2.2 Number of Nonzero Components $\kappa$

For a wave vector $\mathbf{k} = (k_1, k_2, k_3, k_t)$, we define:

$$\kappa_s = |\{i \in \{1, 2, 3\} : k_i \neq 0\}| \quad (\text{number of nonzero spatial components})$$

$$\kappa_t = \begin{cases} 1 & (k_t \neq 0) \\ 0 & (k_t = 0) \end{cases} \quad (\text{temporal nonzero indicator})$$

$$\kappa = \kappa_s + \kappa_t \quad (\text{total number of nonzero components})$$

The possible values of $\kappa$ are $\{0, 1, 2, 3, 4\}$, but by §2.3, $\kappa \geq 3$ admits no stable modes.

### §2.3 Non-Existence of $\kappa \geq 3$

Waves with $\kappa \geq 3$ have simultaneous phase variations in three or more directions. The following two **independently sufficient** reasons show that they do not form stable modes (either one alone excludes $\kappa \geq 3$):

**(i) Conflict with the Nyquist stability condition** (sufficient on its own): By Proposition 9.1 of [2], the phase change along each axis is restricted to $|\phi_i| \leq \pi/2$. The total phase load of a self-returning loop maintaining nonzero phase changes simultaneously on $\kappa$ axes is at least $\kappa \cdot |\phi_{\min}|$, and for $\kappa = 3$, $3 \cdot |\phi_{\min}| > \pi$ (exceeding the stability boundary $\pi$). Therefore, for $\kappa \geq 3$, self-returning loops cannot form stable closed paths.

**(ii) Directional delocalization** (sufficient on its own): Waves with $\kappa_s \geq 2$ and $\kappa_t = 1$ propagate simultaneously in multiple spatial directions. They cannot concentrate displacement in a specific spatial region while satisfying the shape-invariance condition $|\psi(\mathbf{x} + \mathbf{a})| = |\psi(\mathbf{x})|$ in all directions, becoming indistinguishable from the background field.

**Proposition 2.1**: Shape-invariant traveling waves on $\mathbb{Z}^4$ with $\kappa \geq 3$ have no stable modes. Conditions (i) and (ii) are each independently sufficient to establish this. The stable modes of shape-invariant waves are limited to three types: $\kappa \in \{0, 1, 2\}$.

### §2.4 Mode Classification Table

| Mode | $\kappa$ | $\kappa_s$ | $\kappa_t$ | Geometric character |
|------|----------|------------|------------|---------------------|
| Direction-constrained wave | 2 | 1 | 1 | Linear propagation in the (one spatial axis, $t$) plane |
| Planar standing wave ${}^{*}$ | 2 | 2 | 0 | Standing pattern in the (space, space) plane |
| Spherical wave | 0 | 0 | 0 or 1 | No spatial direction; isotropic |
| Standing wave | 1 | 1 | 0 | Standing oscillation along one spatial axis; scalar |

${}^{*}$ The planar standing wave shares $\kappa = 2$ and therefore shares the signed-area property (exclusion principle) of §7 with the direction-constrained wave. However, having $\kappa_t = 0$, it is a standing structure without propagation and is qualitatively different from the direction-constrained wave. Details are given in §3.4.

---

## §3 Direction-Constrained Waves ($\kappa = 2$, $\kappa_s = 1$, $\kappa_t = 1$)

### §3.1 Linear Propagation

For $\kappa = 2$ with $\kappa_s = 1$, $\kappa_t = 1$, the wave vector can be written as $\mathbf{k} = (k_x, 0, 0, k_t)$ (the choice of spatial axis is without loss of generality by the $\mathrm{SO}(3)$ symmetry of $S^3$). From the shape-invariance condition:

$$\psi(\mathbf{x} + \mathbf{e}_1) = e^{ik_x} \cdot \psi(\mathbf{x}) \quad (\text{phase change in the $x$ direction})$$

$$\psi(\mathbf{x} + \mathbf{e}_2) = \psi(\mathbf{x}), \quad \psi(\mathbf{x} + \mathbf{e}_3) = \psi(\mathbf{x}) \quad (\text{invariant in the $y, z$ directions})$$

$$\psi(\mathbf{x} + \mathbf{e}_4) = e^{ik_t} \cdot \psi(\mathbf{x}) \quad (\text{phase change in the $t$ direction})$$

This wave has spatial displacement only in the $x$ direction and no displacement in the $y, z$ directions.

**Proposition 3.1**: A shape-invariant traveling wave with $\kappa = 2$, $\kappa_t = 1$ is a direction-constrained wave that propagates in only one of the three spatial directions. Its displacement is constrained to one spatial dimension.

### §3.2 Propagation Velocity and Displacement Constraint

The propagation velocity is determined by $v = k_x / k_t$. By Proposition 9.1 of [2] (Nyquist stability condition), $v \leq c$.

| Condition | Propagation velocity | Nature of displacement |
|-----------|---------------------|----------------------|
| $|k_x| = |k_t|$ | $v = c$ (maximum) | Light-speed propagation |
| $|k_x| < |k_t|$ | $v < c$ (subluminal) | Massive |
| $k_x = 0$ | $v = 0$ | Degenerates to $\kappa = 1$ (standing) |

The case $v < c$ occurs when $|k_x / k_t| < 1$, i.e., the spatial displacement per temporal component is limited. This displacement constraint is intrinsic as a component ratio of the wave vector itself, not an external parameter.

### §3.3 Generation Structure

Depending on which of the spatial axes $k_1, k_2, k_3$ is nonzero, three equivalent but distinguishable direction-constrained waves exist:

| Nonzero spatial axis | Wave vector |
|---------------------|-------------|
| $k_1$ | $(k_1, 0, 0, k_t)$ |
| $k_2$ | $(0, k_2, 0, k_t)$ |
| $k_3$ | $(0, 0, k_3, k_t)$ |

By the $\mathrm{SO}(3)$ symmetry of $S^3$, these are locally equivalent but are distinguished as wave vector labels. In Table A of [4], these three choices correspond to the generation structure. The following naming convention of "first generation, second generation, third generation" indicates correspondence with the interpretation layer (§9) of [4] and does not affect the validity of the propositions in this paper.

### §3.4 Degenerate Type ($\kappa_s = 2$, $\kappa_t = 0$)

For $\kappa = 2$ with $\kappa_s = 2$, $\kappa_t = 0$, the wave vector is $\mathbf{k} = (k_x, k_y, 0, 0)$. Since the temporal component is zero:

$$\psi(\mathbf{x}, t+1) = \psi(\mathbf{x}, t) \quad (\text{time-invariant})$$

This wave is a standing pattern in space without temporal evolution and is qualitatively different from the direction-constrained wave of §3.1. However, since $\kappa = 2$, it shares the signed-area properties shown in §7.

In Table A of [4], the degenerate type corresponds to neutrinos ($\nu_e, \nu_\mu, \nu_\tau$). Each generation's neutrino shares one spatial axis with the charged fermion of the same generation while placing the other nonzero axis on a different spatial axis rather than the time axis:

| Generation | Charged fermion | Shared axis | Neutrino | Proper axis |
|------------|----------------|-------------|----------|-------------|
| 1 | $(k_1, 0, 0, k_t)$ | $k_1$ | $(k_1, k_2, 0, 0)$ | $k_2$ |
| 2 | $(0, k_2, 0, k_t)$ | $k_2$ | $(0, k_2, k_3, 0)$ | $k_3$ |
| 3 | $(0, 0, k_3, k_t)$ | $k_3$ | $(k_1, 0, k_3, 0)$ | $k_1$ |

---

## §4 Spherical Waves ($\kappa_s = 0$)

### §4.1 Absence of Spatial Directionality

Waves with $\kappa_s = 0$ have all spatial components equal to zero: $\mathbf{k}_s = (0, 0, 0)$. From the shape-invariance condition:

$$\psi(\mathbf{x} + \mathbf{e}_s) = \psi(\mathbf{x}) \quad (\text{phase-invariant in all spatial directions},\; s \in \{1, 2, 3\})$$

A wave with no phase variation in any spatial direction has no preferred direction in space. When emitted from a point source, it propagates as an isotropic spherical wave expanding in all directions.

**Proposition 4.1**: Waves with $\kappa_s = 0$ are isotropic in space and propagate as spherical waves in all directions.

### §4.2 Representation on $\mathbb{Z}^4$

Within the framework of $\mathbb{Z}^4$, $\kappa = 0$ ($\kappa_s = 0$ and $\kappa_t = 0$) reduces to a trivial constant function. Non-trivial spherical waves are realized by having nonzero components on the extended axes introduced in [4] ($R$: radial direction of central projection, $Q$: discrete label axis).

Specifically, a spherical wave with $\kappa_t = 1$ has $\mathbf{k} = (0, 0, 0, k_t)$ and possesses phase variation only in the temporal direction. It is spatially isotropic but undergoes temporal evolution.

### §4.3 Time Symmetry and Retroactive Construction of Causality

A spherical wave with $k_t = 0$ is phase-invariant in the temporal direction as well:

$$\psi(\mathbf{x}, t + 1) = \psi(\mathbf{x}, t) \quad (\text{time-reversal symmetric})$$

**Proposition 4.2**: Waves with $k_t = 0$ are time-reversal symmetric and carry no intrinsic causal ordering of propagation.

As a consequence of this proposition, the interaction between two direction-constrained waves mediated by a spherical wave with $k_t = 0$ has the following properties:

1. The temporal order of emission and absorption is undefined.
2. Only the outcome of the interaction is physically determined.
3. From the outcome, the propagation path along the shortest geodesic is retroactively constructed.

This provides a geometric derivation of the principle of least action: the spherical wave does not "choose" the shortest path; rather, the shortest path is retroactively determined from the outcome of an interaction mediated by a causally unordered spherical wave.

### §4.4 Filling of the Subjective Space as a Ground Mode

Based on the standing wave condition $2\pi R = n \cdot \lambda_{\text{base}}$ from [1], when a spherical wave takes the ground mode ($n = 1$), a wave with wavelength $\lambda = 2\pi R$ extends over the entire $S^3$.

The ground-mode spherical wave with $\kappa_s = 0$:
- Uniformly fills the entire $S^3$
- Has no phase structure in any spatial direction (isotropic)
- Is temporally invariant (when $k_t = 0$) or oscillates uniformly (when $k_t \neq 0$)

This structure functions as the geometric representation of force "fields." A spherical wave emitted from a point source spreads over the entire $S^3$, extending the range of the force to the whole subjective space.

### §4.5 CP Symmetry

In Table A of [4], spherical waves with $k_t = 0$ (photon-type, gluon-type) are distinguished from those with $k_t \neq 0$ ($W^{\pm}$-type). By Proposition 4.2:

**Proposition 4.3**: CP symmetry can be violated only in interactions mediated by spherical waves with $k_t \neq 0$. Interactions mediated by spherical waves with $k_t = 0$ are automatically CP-symmetric.

In Table A of [4], the only gauge boson with $k_t \neq 0$ is $W^{\pm}$, which is consistent with the fact that CP violation is observed only in weak interactions. The naming conventions such as "photon-type" and "$W^{\pm}$-type" here indicate correspondence with the interpretation layer of [4]; Proposition 4.3 itself relies solely on the geometric condition of whether $k_t = 0$.

---

## §5 Standing Waves ($\kappa = 1$, $\kappa_s = 1$, $\kappa_t = 0$)

### §5.1 Standing Oscillation Along One Spatial Axis

For $\kappa = 1$ with $\kappa_s = 1$, $\kappa_t = 0$, the wave vector is $\mathbf{k} = (k_x, 0, 0, 0)$:

$$\psi(\mathbf{x} + \mathbf{e}_1) = e^{ik_x} \cdot \psi(\mathbf{x}) \quad (\text{phase change in the $x$ direction})$$

$$\psi(\mathbf{x} + \mathbf{e}_s) = \psi(\mathbf{x}) \quad (s \neq 1,\; \text{invariant in other directions})$$

$$\psi(\mathbf{x}, t + 1) = \psi(\mathbf{x}, t) \quad (\text{time-invariant})$$

This wave has an amplitude distribution along one spatial axis and is a standing structure that does not change in time.

### §5.2 Scalar Character

The spatial axes $k_1, k_2, k_3$ are equivalent by the $\mathrm{SO}(3)$ symmetry of $S^3$. While $\mathbf{k} = (k_x, 0, 0, 0)$ selects the $x$-axis, $\mathbf{k} = (0, k_y, 0, 0)$ and $\mathbf{k} = (0, 0, k_z, 0)$ exist equally.

The average over the three equivalent configurations is isotropic, which is the condition for a scalar field.

**Proposition 5.1**: Waves with $\kappa = 1$, $\kappa_t = 0$ define a scalar quantity by virtue of the equivalence of three directions under $\mathrm{SO}(3)$ symmetry.

### §5.3 Filling of the Subjective Space by the Ground Mode

By the standing wave condition from [1], the ground mode is a wave with wavelength $\lambda = 2\pi R$ that extends over the entire $S^3$. When a standing wave with $\kappa = 1$ takes the ground mode:

- It fills the entire $S^3$ with one wavelength
- It is time-invariant ($k_t = 0$)
- It has phase structure along only one spatial axis
- It is scalar by $\mathrm{SO}(3)$ averaging

This is the description of the simplest non-trivial standing wave filling a closed space and functions as the geometric representation of vacuum structure.

### §5.4 Uniqueness of $\kappa = 1$

$\kappa = 0$ is trivial (constant function) within $\mathbb{Z}^4$ and has no spatial phase structure. $\kappa = 2$ inherently possesses surface orientation (rotational structure) as shown in §6 and is not scalar.

**Proposition 5.2**: The minimal structure for a spatially non-trivial and scalar (non-rotational) shape-invariant wave is $\kappa = 1$, which is the unique scalar mode.

### §5.5 Comparison with the Spherical Wave Mode

The spherical wave with $\kappa = 0$ (§4.4) and the standing wave with $\kappa = 1$ both fill the entire $S^3$ as ground modes, but they are qualitatively different:

| | Spherical wave ($\kappa = 0$) | Standing wave ($\kappa = 1$) |
|--|-------------------------------|------------------------------|
| Spatial phase structure | None | Exists along one axis |
| Spin | 0 or higher (depends on extended axes) | 0 (scalar) |
| Propagation | Can spread spherically | Standing (does not propagate) |
| Role | Force mediation (field) | Vacuum structure (ground state) |

---

## §6 Geometric Origin of Chirality

### §6.1 Two-Axis Plane and Winding Direction

Waves with $\kappa = 2$ possess a plane defined by two nonzero axes. On this plane, the winding direction of the phase field (clockwise / counterclockwise) is defined. For the winding direction to be physically distinguishable, the orientation of the plane must be fixed externally.

### §6.2 (Space, $t$) Plane: Definition of Chirality

For $\mathbf{k} = (k_x, 0, 0, k_t)$, on the $(x, t)$ plane:

- The $t$-axis has its future direction fixed (arrow of time)
- The $x$-axis is a spatial coordinate (selectable by $\mathrm{SO}(3)$ rotation, but fixed once selected)

Since the direction of $t$ is fixed, the two winding directions of the phase on the $(x, t)$ plane are physically distinguishable:

$$\chi = \mathrm{sign}(k_x \cdot k_t) = \begin{cases} +1 & (\text{right-handed}) \\ -1 & (\text{left-handed}) \end{cases}$$

**Definition 6.1**: For a wave with $\kappa = 2$, $\kappa_s = 1$, $\kappa_t = 1$, the sign product $\chi = \mathrm{sign}(k_{s_0} \cdot k_t)$ of the nonzero spatial component $k_{s_0}$ (the unique spatial component with $k_{s_0} \neq 0$) and the temporal component $k_t$ is defined as the **chirality**.

**Proposition 6.1**: Waves with $\kappa = 2$ possessing a (space, $t$) plane have two states, left-handed and right-handed, by virtue of the fixed direction of the $t$-axis. This distinction cannot be removed by continuous spacetime transformations.

### §6.3 (Space, Space) Plane: Degeneracy of Chirality

For $\mathbf{k} = (k_x, k_y, 0, 0)$, on the $(x, y)$ plane:

- Both the $x$ and $y$ axes are spatial coordinates
- An $\mathrm{SO}(2)$ rotation $R(\pi)$ can reverse the orientation of the plane

Since the orientation of the $(x, y)$ plane can be freely changed by continuous rotation, clockwise and counterclockwise are physically indistinguishable.

**Proposition 6.2**: Waves with $\kappa = 2$ possessing only (space, space) planes have degenerate chirality; left-handed and right-handed are physically indistinguishable.

### §6.4 Majorana Condition

By the Feynman--Stuckelberg interpretation, an antiparticle is described as a particle traveling backward in time. Waves with $\kappa_t = 0$ are invariant under time reversal:

$$\psi(\mathbf{x}, -t) = \psi(\mathbf{x}, t)$$

The forward-in-time wave and the backward-in-time wave describe the same state.

**Proposition 6.3**: Waves with $\kappa = 2$, $\kappa_t = 0$ are time-reversal symmetric, and the particle state and antiparticle state are described by the same wave function (Majorana condition).

In Table A of [4], the only fermions with $\kappa_t = 0$ are neutrinos. Proposition 6.3 geometrically suggests the possibility that neutrinos are of Majorana type. The naming "neutrino" here indicates correspondence with the interpretation layer of [4]; Proposition 6.3 itself is derived solely from the geometric conditions $\kappa = 2$, $\kappa_t = 0$.

### §6.5 Consistency with Discrete Transformations

We define four discrete transformations as operations on the wave vector:

| Transformation | Action on wave vector | Chirality $\chi$ |
|----------------|----------------------|-------------------|
| Parity $P$ | $k_x \to -k_x$, $k_t$ unchanged | Reversed |
| Time reversal $T$ | $k_t \to -k_t$, $k_x$ unchanged | Reversed |
| Charge conjugation $C$ | $(k_x, k_t) \to (-k_x, -k_t)$ | **Unchanged** |
| $CPT$ | All components reversed then restored | **Unchanged** |

**Proposition 6.4**: The $CPT$ transformation is the reversal of all components of the wave vector and holds automatically as a lattice symmetry of $\mathbb{Z}^4$.

---

## §7 Geometric Origin of the Pauli Exclusion Principle

### §7.1 Signed Area

The two-axis plane defined by a wave with $\kappa = 2$ possesses an orientation. The orientation of a two-dimensional plane has **only two possibilities: positive (front) and negative (back)**.

Consider the winding number $n \in \mathbb{Z}$ of the phase field. By Proposition 9.1 of [2] (Nyquist stability condition), the phase change along each axis is restricted to $|\phi_i| \leq \pi/2$. The total phase change for one cycle on the two-axis plane is restricted to $|\phi_1| + |\phi_2| \leq \pi$, and a phase field with winding number $n$ requires a phase change of $2\pi |n|$ per cycle. From $|2\pi n| \leq \pi$, we get $|n| \leq 1/2$, and combined with the condition $n \in \mathbb{Z}$, $n$ is limited to $n \in \{-1, 0, +1\}$ ($n = 0$ is trivial). Therefore, fundamental shape-invariant waves have $n = \pm 1$, giving:

$$n = +1 \quad \Rightarrow \quad \text{signed area} = +A$$

$$n = -1 \quad \Rightarrow \quad \text{signed area} = -A$$

### §7.2 Maximum Coexistence Number of Identical States

Consider a situation where multiple waves with $\kappa = 2$ sharing the same wave vector $\mathbf{k}$ coexist in the same lattice region:

- First wave: signed area $+A$
- Second wave: signed area $-A$ (opposite orientation)
- Third wave: Only two orientations, $+A$ and $-A$, exist for the two-dimensional plane. Both are already occupied, and no independent third state exists.

**Proposition 7.1** (Exclusion principle): Shape-invariant waves with $\kappa = 2$ sharing the same wave vector can coexist in the same lattice region up to a maximum of two, because only two orientational states ($\pm$) exist for the two-axis plane.

### §7.3 Cases with $\kappa \neq 2$

$\kappa = 0$: No nonzero axes exist in space $\to$ no two-axis plane can be formed $\to$ no orientational constraint.

$\kappa = 1$: Only one axis $\to$ one axis short of forming a plane $\to$ no orientational constraint.

**Proposition 7.2**: Waves with $\kappa \neq 2$ have no signed area, and arbitrarily many waves with the same wave vector can coexist.

### §7.4 Geometric Correspondence of Spin-Statistics

Combining Propositions 7.1 and 7.2:

| $\kappa$ | Signed area | Maximum coexistence number | Statistics |
|----------|-------------|---------------------------|-----------|
| 2 | Present ($\pm$) | 2 | Fermi-type |
| 0, 1 | Absent | Unlimited | Bose-type |

The geometric correspondence of the spin-statistics theorem -- half-integer spin obeys Fermi statistics, integer spin obeys Bose statistics -- is obtained from the fact that "a two-dimensional plane has only two orientations."

### §7.5 Pair Condensation

A pair consisting of a wave with signed area $+A$ and a wave with $-A$ has zero net signed area:

$$(+A) + (-A) = 0$$

A pair with zero area is not subject to the constraint of Proposition 7.1.

**Proposition 7.3**: A pair of waves with $\kappa = 2$ (a combination of areas $\pm A$) has zero net signed area, and as a pair unit, is not subject to the coexistence constraint.

This corresponds to the geometric description of Cooper pair condensation in superconductors (the phenomenon where pairs of fermions exhibit bosonic behavior).

---

## §8 Child-Space Central Projection and Acceleration

### §8.1 Uniform Motion and Acceleration

We restate the discrete flow on the tangent bundle $\mathbb{Z}^4 \times \mathbb{Z}^4$ introduced in §8 of [2]:

$$\mathbf{k}_{t+1} = \mathbf{k}_t + \Delta\mathbf{k}_t \quad (\text{position update})$$

$$\Delta\mathbf{k}_{t+1} = \Delta\mathbf{k}_t + \Delta^2\mathbf{k}(\mathbf{k}_t) \quad (\text{velocity update})$$

The first equation describes uniform motion, and $\Delta^2\mathbf{k}$ in the second equation corresponds to acceleration. When $\Delta^2\mathbf{k} = 0$, the motion is uniform and rectilinear; when $\Delta^2\mathbf{k} \neq 0$, it is accelerated motion.

### §8.2 Global Acceleration: Curvature of the Subjective Space

As shown in [1], shape-invariant standing waves exist on $S^3$, and the curvature radius $R$ of $S^3$ is not an externally given constant but a quantity self-consistently determined from the ground-mode structure of the standing waves present in the system ([1] §3.2: the standing wave condition $2\pi R = n \cdot \lambda$ determines $R$ a posteriori). Geodesics (great circles) on $S^3$ follow the curvature. When these geodesics are centrally projected onto $\mathbb{Z}^4$, they appear as curves, and their curvature gives rise to the global acceleration:

$$\Delta^2\mathbf{k}_{\text{global}} = f(R, \mathbf{k}_t)$$

Here $f$ is a function depending on the curvature scalar $R$ of $S^3$ and the current position $\mathbf{k}_t$, determined from the Jacobian of the central projection.

### §8.3 Local Acceleration: Child-Space Central Projection

We utilize the central projection structure shown in §10 of [2]. In the central projection $\pi : S^5 \to \mathbb{Z}^4$ from the parent space $S^5(R)$ to the child space $\mathbb{Z}^4$:

- Wave trajectories in the parent space are geodesics (segments of great circles)
- The central projection is a fractional (nonlinear) map
- In the child space, geodesics appear as curves

The local acceleration is determined by the nonlinearity of the central projection:

$$\Delta^2\mathbf{k}_{\text{local}}(\mathbf{k}_t) = D^2\pi(\mathbf{k}_t) \cdot \Delta\mathbf{k}_t$$

Here $D^2\pi$ is the Hessian matrix of the central projection $\pi$, evaluated discretely on $\mathbb{Z}^4$.

**Proposition 8.1**: The second derivative $D^2\pi$ of the child-space central projection geometrically determines the acceleration of shape-invariant wave trajectories. This acceleration is not a force imposed externally but arises intrinsically from the curvature of the projection.

### §8.4 Geometric Unification

Global acceleration (curvature of $S^3$) and local acceleration (nonlinearity of the central projection) are descriptions of the same geometric mechanism at different scales:

- **Global**: The energy distribution of many waves determines the curvature radius $R$ of $S^3$ $\to$ $R$ affects all wave trajectories (integral effect)
- **Local**: The Hessian matrix $D^2\pi$ of the projection in the neighborhood of an individual wave produces local trajectory deflection (differential effect)

Both are connected as two aspects of a single mechanism because the central projection $\pi: S^5(R) \to \mathbb{Z}^4$ is a map that parametrically depends on the curvature radius $R$, and $R$ itself is self-consistently determined from the distribution of waves (§8.2). The determination of $R$ (global) and the evaluation of $D^2\pi$ at each point (local) are different aspects of the same $R$-parameterized map $\pi$.

**Proposition 8.2**: Acceleration is not an additional term in the wave equation but is intrinsically derived from the geometry of the central projection. Global curvature (determination of $R$) and local projection nonlinearity (point-by-point evaluation of $D^2\pi$) are manifestations at different scales of the single $R$-parameterized central projection map.

---

## §9 Summary

We summarize the main results of this paper:

1. **Mode classification**: Shape-invariant waves are classified into three types by the number of nonzero components $\kappa$ of the wave vector: $\kappa \in \{0, 1, 2\}$. No stable modes exist for $\kappa \geq 3$ (§2).

2. **Direction-constrained waves** ($\kappa = 2$): They propagate in only one spatial direction, and the signed area of the two-axis plane exists. The choice of spatial axis gives rise to the generation structure (§3).

3. **Spherical waves** ($\kappa = 0$): They have no spatial directionality and propagate isotropically. When $k_t = 0$, they are time-reversal symmetric and carry no intrinsic causal ordering. As a ground mode, they fill the entire subjective space and function as force fields (§4).

4. **Standing waves** ($\kappa = 1$): They are standing oscillations along one spatial axis and constitute the unique scalar mode by $\mathrm{SO}(3)$ symmetry. As a ground mode, they fill the entire subjective space and function as vacuum structure (§5).

5. **Chirality**: In the (space, $t$) plane, left-handed and right-handed states are distinguished by the direction of the $t$-axis. In the (space, space) plane, chirality is degenerate and the Majorana condition is satisfied. CP violation is limited to interactions with $k_t \neq 0$ (§6).

6. **Exclusion principle**: Because the orientation of a two-dimensional plane has only two states, at most two waves with $\kappa = 2$ can coexist, while waves with $\kappa \neq 2$ have no such restriction. Pairs ($\pm A$) have zero area and are bosonic (§7).

7. **Acceleration**: The second derivative of the child-space central projection geometrically determines acceleration, unifying global curvature and local nonlinearity (§8).

The correspondence with Table A of [4] is shown below:

| Mode | $\kappa$ | Table A correspondence | Spin | Statistics |
|------|----------|----------------------|------|-----------|
| Direction-constrained wave | 2 ($\kappa_s = 1, \kappa_t = 1$) | Charged fermions | 1/2 | Fermi |
| Planar standing wave | 2 ($\kappa_s = 2, \kappa_t = 0$) | Neutrinos | 1/2 | Fermi (Majorana) |
| Spherical wave | 0 | Gauge bosons ($\gamma, Z, W^{\pm}, g, G$) | 0, 1, 2 | Bose |
| Standing wave | 1 | Higgs $H^0$ | 0 | Bose |

All claims in this paper are stated as geometric facts, and the correspondence with Table A is presented as a "connectable interpretation." Identification with physical particles does not affect the validity of the propositions themselves.

---

## References

[1] Noriaki Kihara, "Shape-Invariant Standing Waves on Closed Spheres: Stability by Dimension and the Geometric Necessity of 3+1 Spacetime," 2026.

[2] Noriaki Kihara, "Shape-Invariant Traveling Waves on the 4-Dimensional Integer Lattice: Existence Conditions, Conservation Structures, and Self-Consistent Loops," 2026.

[3] Noriaki Kihara, "A Three-Layer Model of the Universe via Central Projection," Zenodo, 2025. DOI: 10.5281/zenodo.15083729.

[4] Noriaki Kihara, "Set Structure and Combinatorial Properties of the 6-Dimensional Hyperrectangle," Zenodo, 2025. DOI: 10.5281/zenodo.19688521.
