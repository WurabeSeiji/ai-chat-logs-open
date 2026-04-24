# Interactions of Shape-Invariant Waves: Axis-Directional Displacement Transfer, Wave Packet Deformation, and Retroactive Construction of Causality

**Author**: Noriaki Kihara  
**Affiliation**: WF System Co., Ltd.  
**ORCID**: 0009-0004-6753-4020  
**Version**: v2.1  
**Date**: April 2026  
**DOI**: 10.5281/zenodo.19709800

---

## Abstract

Building on the three-mode classification of shape-invariant waves introduced in the companion paper [3] (direction-constrained waves $\kappa = 2$, spherical waves $\kappa = 0$, standing waves $\kappa = 1$), we provide a unified description of interactions between different modes. The fundamental mechanism of interaction is defined as **axis-directional displacement transfer** -- the change of a direction-constrained wave's state vector $(\mathbf{k}, \Delta\mathbf{k})$ along the direction of the non-zero axis carried by a spherical wave. Interactions are classified into four fundamental vertex types (absorption, emission, pair annihilation, pair creation), and we prove that under conservation laws, absorption $\Leftrightarrow$ emission and pair annihilation $\Leftrightarrow$ pair creation are reversible. Within this unified framework, we systematically describe all fundamental vertices for the four forces (electromagnetism, strong force, weak force, gravity). Coupling with a standing wave ($\kappa = 1$) modifies the propagation speed ratio $|k_x / k_t|$ of a direction-constrained wave, providing a geometric origin for mass. Under the conserved quantity $\int |\psi|^2 \, dV$, energy transfer deforms the amplitude and width of a wave packet inversely, and this deformation gives a dynamical description of "wave packet collapse." Interactions mediated by spherical waves with $k_t = 0$ carry no intrinsic causal ordering, and the shortest geodesic path and the localization of the wave packet are retroactively constructed from the interaction outcome.

**Keywords**: shape-invariant wave, fundamental vertex, axis-directional displacement transfer, wave packet deformation, origin of mass, retroactive construction of causality, conservation law

---

## S1 Introduction and Scope of Claims

This paper provides a unified description of interactions between waves of different modes within the framework of shape-invariant waves constructed in the companion papers [1]--[3], and as consequences derives the origin of mass, the classification of forces, wave packet deformation, and causal structure geometrically.

The claims of this paper are limited to the following six points:

1. All interactions are described uniformly as displacement transfer along the non-zero axes of spherical waves. The difference among types of force reduces solely to the difference in the axis along which the transfer occurs.
2. Interactions are classified into four fundamental vertex types (absorption, emission, pair annihilation, pair creation), and under conservation laws, absorption $\Leftrightarrow$ emission and pair annihilation $\Leftrightarrow$ pair creation are reversible.
3. Coupling with a standing wave ($\kappa = 1$) modifies the speed ratio of a direction-constrained wave, providing a geometric origin for mass.
4. Under the conserved quantity $\int |\psi|^2 \, dV$, energy transfer deterministically deforms the shape of a wave packet, and this deformation corresponds to "wave packet collapse" in observation.
5. Interactions mediated by spherical waves with $k_t = 0$ carry no intrinsic causal ordering, and the path and localization are retroactively constructed from the outcome.
6. The sign of the interaction (attraction / repulsion) between direction-constrained waves is determined geometrically by the relative sign of their $t$-axis components.

The following elements are **outside the scope of this paper**:

Quantitative derivations:

- Numerical derivation of coupling constants (specific values of $\alpha$, $G_F$, $\alpha_s$, $G$)
- Quantitative computation of scattering cross sections and decay widths
- Renormalization and scale dependence of running coupling constants
- Rigorous reconstruction of gauge groups in the continuum limit

Details of specific mechanisms:

- Mass acquisition mechanism for $W^{\pm}$ and $Z^0$ bosons (details of indirect coupling not mediated by spatial axes)
- Mass acquisition mechanism for 1st- and 2nd-generation fermions (details of indirect coupling when the non-zero axis of the standing wave is not directly shared)

**Note on physical nomenclature**: This paper uses physical names (photon, graviton, electron, flavor, color charge, etc.) assuming correspondence with Table A of [4], but this is for interpretive convenience; the geometric content of each proposition does not depend on the physical identification. All propositions are derived solely from the component structure of wave vectors ($\kappa$, $k_t$, type of axis).

---

## S2 Unified Definition of Interaction

### S2.1 Two-Layer Structure

In the tangent bundle $\mathbb{Z}^4 \times \mathbb{Z}^4$ introduced in [2] S8, the state of a direction-constrained wave ($\kappa = 2$) is:

$$(\mathbf{k},\; \Delta\mathbf{k}) \in \mathbb{Z}^4 \times \mathbb{Z}^4$$

- $\mathbf{k}$: **Static label** (type of wave, corresponding to Table A in [4])
- $\Delta\mathbf{k}$: **Dynamic state** (corresponding to energy-momentum)

### S2.2 Axis-Directional Displacement Transfer

**Definition 2.1** (Axis-directional displacement transfer): When a spherical wave ($\kappa = 0$) interacts with a direction-constrained wave ($\kappa = 2$), the state of the direction-constrained wave changes as follows:

$$(\mathbf{k},\; \Delta\mathbf{k}) \;\longrightarrow\; (\mathbf{k} + \delta\mathbf{k},\; \Delta\mathbf{k} + \delta(\Delta\mathbf{k}))$$

where the displacements $\delta\mathbf{k}$ and $\delta(\Delta\mathbf{k})$ satisfy the following conditions:

**(I1) Axiality**: The non-zero components of $\delta\mathbf{k}$ and $\delta(\Delta\mathbf{k})$ exist only along the non-zero axis directions in Table A of [4] for the spherical wave. Here, "non-zero axis directions" includes not only continuous displacements along the $t$ and $R$ axes but also discrete label transitions on the $Q$ axis ($Q_i \to Q_j$). That is, for gluons, (I1) requires that the displacement transfer acts only on the $Q$ label change, preserving all other axis components.

**(I2) Conservation law**: The total sum of tangent bundle states of all waves involved in the interaction is conserved:

$$\sum_i (\mathbf{k}_i + \Delta\mathbf{k}_i) = \text{const}$$

**(I3) Wave packet conservation**: The conserved quantity $\int |\psi|^2 \, dV$ of each wave that exists before and after the interaction is individually conserved. For waves newly created by pair creation (Type D), the initial value of the conserved quantity is set at the time of creation. For waves that are annihilated by pair annihilation (Type C), their conserved quantity is transferred to other waves at the time of annihilation.

### S2.3 Static Label Change and Dynamic State Change

Axis-directional displacement transfer encompasses two qualitatively different types of change:

| Type of change | Target | Physical meaning | Example |
|---------------|--------|-----------------|---------|
| $\delta\mathbf{k} \neq 0$ | Static label $\mathbf{k}$ | Conversion of wave type | Flavor conversion, color charge exchange |
| $\delta(\Delta\mathbf{k}) \neq 0$ | Dynamic state $\Delta\mathbf{k}$ | Change in energy-momentum | Scattering, acceleration |

The type of spherical wave determines which (or both) is dominant.

### S2.4 Unification of the Four Forces

From the non-zero axis structure of spherical waves in Table A of [4], the four forces are described as applications of the same operation to different axes:

| Force | Spherical wave | Non-zero axis | Dominant change |
|-------|---------------|--------------|-----------------|
| Electromagnetism | $\gamma\; (0,0,0,0,R,0)$ | $R$ | $\delta(\Delta\mathbf{k})$ (momentum change) |
| Strong force | $g\; (0,0,0,0,0,Q_i)$ ${}^{*}$ | $Q$ | $\delta\mathbf{k}$ (label change) |
| Weak force | $W^{\pm}\; (0,0,0,\pm t,0,0)$ | $t$ | $\delta\mathbf{k}$ (label change) |
| Gravity | $G\; (0,0,0,\pm t,\pm R,0)$ | $t, R$ | $\delta(\Delta\mathbf{k})$ (energy change) |

${}^{*}$ Whereas $\gamma$, $W^{\pm}$, and $G$ are spherical waves with static non-zero components on the $t$ or $R$ axes, $g$ mediates color charge exchange through a spherical wave with discrete label value $Q_i$ on the $Q$ axis transitioning to a different label value $Q_j$. That is, the "non-zero axis" of $g$ acts not as a fixed value but as a change in the $Q$ label (detailed in S5).

**Proposition 2.1**: All four forces are described as axis-directional displacement transfer per Definition 2.1, and the difference among types of force reduces solely to the difference in the axis along which the transfer occurs.

### S2.5 $R$-Axis Discretization Condition (Prerequisite)

This paper treats direction-constrained waves as integer lattice waves on $\mathbb{Z}^4$ (S2.1). For this framework to hold, the curvature radius $R$ must satisfy a number-theoretic constraint.

Theorem 5.2 of prior work [W5] showed that for a 4-dimensional subjective space to be densely regularized with integer values under 5-dimensional central projection, $R^2$ (in units of the lattice constant) must be a positive integer ($R^2 \in \mathbb{Z}_{>0}$). This is a number-theoretic constraint based on Lagrange's four-square theorem, and if this condition is violated, it becomes impossible to close all phases at integer values.

In [4] S9.8, the discretization $R = \pm n^2 \cdot L_0$ ($n \in \mathbb{N}$, $L_0$ is the lattice constant) is derived from this condition. That is, the $R$ axis takes not continuous real values but only **discrete values of perfect squares**.

**(I4) Consequence of $R$-axis discreteness**: Under the above prerequisite ($R^2 \in \mathbb{Z}_{>0}$), the values that the displacement $\delta(\Delta k)_R$ along the $R$-axis direction can take are automatically constrained by the requirement that the $R$ values before and after the displacement both satisfy the perfect-square condition. That is, $\delta(\Delta k)_R$ is restricted to the form of the difference of two perfect squares: $n_1^2 - n_2^2 = (n_1 - n_2)(n_1 + n_2)$.

The specific consequences of this discretization condition for individual propositions in this paper — in particular, how it constrains the allowed energy levels in S4 (electromagnetic interaction involving $R$-axis displacement) and S7 (gravitational interaction involving $(t, R)$-axis displacement) — remain unresolved and are outside the scope of this paper. This section is placed as a preventive reference point for tracking the consequences of the number-theoretic constraint in subsequent research. For details of the discretization condition, see [W5] S5.5 and [4] S9.8.

---

## S3 Classification and Reversibility of Fundamental Vertices

### S3.1 Four Fundamental Vertex Types

The primitive patterns of interaction (fundamental vertices) are classified into the following four types. $f$ denotes a direction-constrained wave ($\kappa = 2$) and $B$ denotes a spherical wave ($\kappa = 0$).

**Type A (Absorption)**: $f + B \longrightarrow f'$

A direction-constrained wave $f$ absorbs a spherical wave $B$. $B$ ceases to exist as an independent wave, and the $(\mathbf{k}_B + \Delta\mathbf{k}_B)$ carried by $B$ contributes to the state change of $f$:

$$(\mathbf{k}_f + \Delta\mathbf{k}_f) + (\mathbf{k}_B + \Delta\mathbf{k}_B) = (\mathbf{k}_{f'} + \Delta\mathbf{k}_{f'})$$

**Type B (Emission)**: $f \longrightarrow f' + B$

A direction-constrained wave $f$ emits a spherical wave $B$. $B$ is newly created as an independent wave, and $(\mathbf{k}_B + \Delta\mathbf{k}_B)$ of $B$ is supplied from the state change of $f$:

$$(\mathbf{k}_f + \Delta\mathbf{k}_f) = (\mathbf{k}_{f'} + \Delta\mathbf{k}_{f'}) + (\mathbf{k}_B + \Delta\mathbf{k}_B)$$

**Type C (Pair Annihilation)**: $f + \bar{f} \longrightarrow B(s)$

A direction-constrained wave $f$ and its anti-wave $\bar{f}$ annihilate, generating one or more spherical waves $B$. This is possible because the sum of static labels becomes zero:

$$\mathbf{k}_f + \mathbf{k}_{\bar{f}} = \mathbf{0} \;\Longrightarrow\; \kappa = 2 \to \kappa = 0$$

**Type D (Pair Creation)**: $B(s) \longrightarrow f + \bar{f}$

From one or more spherical waves, a pair $(f, \bar{f})$ of direction-constrained waves is created. A pair of static labels $(\mathbf{k}, -\mathbf{k})$ is created from zero:

$$\kappa = 0 \to \kappa = 2 \;\Longleftarrow\; \mathbf{0} \to \mathbf{k}_f + \mathbf{k}_{\bar{f}} = \mathbf{0}$$

### S3.2 Reversibility Theorem

**Proposition 3.1** (Reversibility): For any process $\mathcal{P}$ satisfying conditions (I1)--(I3), the time-reversed process $\bar{\mathcal{P}}$ also satisfies (I1)--(I3).

*Proof*: (I1) is a structural condition on axis direction and does not depend on the time direction. (I2) is conservation of a sum, which holds when the initial and final states are exchanged. (I3) is conservation of $\int |\psi|^2 \, dV$ for each wave and does not depend on the direction of the process. $\square$

**Corollary 3.1**: Type A (absorption) and Type B (emission) are reverse processes of each other, and Type C (pair annihilation) and Type D (pair creation) are reverse processes of each other.

**Remark 3.1**: When spherical waves with $k_t \neq 0$ ($W^{\pm}$, $G$) are involved, both the forward and reverse processes are permitted, but due to the temporal asymmetry introduced by $k_t \neq 0$, a causal temporal ordering exists (detailed in S9.5).

### S3.3 Composite Processes

More complex interactions are constructed by composing fundamental vertices.

**Exchange**: Composition of Type B + Type A. $f_1$ emits $B$, and $f_2$ absorbs $B$:

$$f_1 + f_2 \;\xrightarrow{B\text{ emission}}\; f_1' + B + f_2 \;\xrightarrow{B\text{ absorption}}\; f_1' + f_2'$$

**Scattering**: Composition of Type A + Type B. $f$ absorbs $B_1$ and emits $B_2$:

$$f + B_1 \;\xrightarrow{\text{absorption}}\; f'' \;\xrightarrow{\text{emission}}\; f' + B_2$$

**Proposition 3.2**: All interactions are described as finite compositions of the four fundamental vertex types (Types A--D).

---

## S4 Electromagnetic Interaction

The photon $\gamma = (0, 0, 0, 0, +, 0)$ is a spherical wave with a non-zero component on the $R$ axis. The $R$ axis corresponds to the radial direction of the central projection and acts as a change in the dynamic state $\Delta\mathbf{k}$ on $\mathbb{Z}^4$. By the discretization condition on the $R$ axis (S2.5, (I4)), displacement transfer along the $R$ direction via the photon takes only discrete values. As a candidate consequence of (I4), only discrete allowed values may exist for the energy-momentum carried by the photon.

### S4.1 Absorption: $f + \gamma \to f'$

When a direction-constrained wave (electron type) $(k_x, 0, 0, k_t)$ absorbs a photon:

$$(\mathbf{k},\; \Delta\mathbf{k}) + (\mathbf{k}_\gamma,\; \Delta\mathbf{k}_\gamma) \;\longrightarrow\; (\mathbf{k},\; \Delta\mathbf{k} + \delta(\Delta\mathbf{k})_\gamma)$$

The photon ceases to exist as an independent wave, and its energy-momentum is transferred to the dynamic state of the direction-constrained wave. The static label $\mathbf{k}$ does not change (the electron remains an electron).

### S4.2 Emission: $f \to f' + \gamma$

By Proposition 3.1 (reversibility), emission is defined as the reverse process of absorption:

$$(\mathbf{k},\; \Delta\mathbf{k}) \;\longrightarrow\; (\mathbf{k},\; \Delta\mathbf{k} - \delta(\Delta\mathbf{k})_\gamma) + (\mathbf{k}_\gamma,\; \delta(\Delta\mathbf{k})_\gamma)$$

A portion of the dynamic state of the direction-constrained wave separates, and a new photon is created as an independent wave. The static label $\mathbf{k}$ is conserved.

**Proposition 4.1**: Interactions with the photon (both absorption and emission) preserve the static label of the direction-constrained wave and change only the dynamic state $\Delta\mathbf{k}$.

### S4.3 Exchange and Attraction/Repulsion

The electromagnetic interaction between two direction-constrained waves $\psi_1$, $\psi_2$ is described as photon exchange (S3.3):

$$\psi_1 \;\xrightarrow{\gamma\text{ emission}}\; \psi_1' + \gamma \;\xrightarrow{\gamma\text{ absorption}}\; \psi_1' + \psi_2'$$

By conservation law (I2):

$$\Delta\mathbf{k}_1 + \Delta\mathbf{k}_2 = \Delta\mathbf{k}_1' + \Delta\mathbf{k}_2'$$

Consider the action of $\delta(\Delta\mathbf{k})_\gamma$ on the $t$-axis components in photon exchange. The photon has a non-zero component on the $R$ axis ((I1)), but the signs of $k_{t,1}$ and $k_{t,2}$ of the two direction-constrained waves involved in the interaction constrain the direction of momentum change under conservation law (I2). In the process where $\psi_1$ emits a photon and $\psi_2$ absorbs it, $\delta(\Delta\mathbf{k})_\gamma$ is transferred from $\psi_1$ to $\psi_2$. When $k_{t,1}$ and $k_{t,2}$ have the same sign, the relative momentum change under (I2) conservation is directed to separate the two waves (repulsion); when the signs are opposite, it is directed to bring them together (attraction).

**Concrete example**: When $k_{t,1}, k_{t,2}$ are both positive, $\psi_1$ emits a photon with $R$-axis displacement $\delta(\Delta k)_R$ and $\psi_2$ absorbs it; by (I2), $\Delta k_{R,1}$ decreases by $\delta(\Delta k)_R$ while $\Delta k_{R,2}$ increases by the same amount. Two waves with the same sign of $k_t$ have co-directional temporal evolution in the radial structure of the central projection ([5] S5), and this momentum redistribution acts spatially to separate the two waves (repulsion). For opposite signs, the contra-directional temporal evolution causes the redistribution to act in the direction of attraction.

The $t$-axis component $k_t$ corresponds to the charge structure in Table A of [4]:

| $k_t$ of $\psi_1$ | $k_t$ of $\psi_2$ | Sign product | Interaction |
|-------------------|-------------------|-------------|-------------|
| $+$ | $+$ | $+1$ | Repulsion |
| $+$ | $-$ | $-1$ | Attraction |
| $-$ | $-$ | $+1$ | Repulsion |

**Proposition 4.2**: The sign of the electromagnetic interaction between two direction-constrained waves is determined geometrically by the relative sign $\mathrm{sign}(k_{t,1} \cdot k_{t,2})$ of their $t$-axis components. The specific correspondence "same sign = repulsion, opposite sign = attraction" depends on the interpretive layer of the correspondence between $k_t$ and charge sign in Table A of [4]; the geometric content of the proposition is that "the relative sign uniquely determines the directionality of the interaction."

### S4.4 Pair Annihilation: $f + \bar{f} \to \gamma(s)$

A direction-constrained wave and its antiparticle are in a relationship where all components have opposite signs ($C$ transformation in [3] S6.5):

$$\psi: (k_x, 0, 0, k_t) \quad \longleftrightarrow \quad \bar{\psi}: (-k_x, 0, 0, -k_t)$$

When both interact in the same lattice region, the sum of static labels becomes zero:

$$(k_x, 0, 0, k_t) + (-k_x, 0, 0, -k_t) = (0, 0, 0, 0)$$

A pair of $\kappa = 2$ waves transitions to a $\kappa = 0$ state. By conservation law (I2), the energy of the vanished static labels is emitted as spherical waves (photons).

**Proposition 4.3**: Pair annihilation of a direction-constrained wave and an anti-direction-constrained wave is a process in which a pair of $\kappa = 2$ waves is converted to $\kappa = 0$ spherical waves, made possible by the sum of static labels becoming zero.

### S4.5 Pair Creation: $\gamma(s) \to f + \bar{f}$

By Proposition 3.1 (reversibility), pair creation is defined as the reverse process of pair annihilation. From the energy of a spherical wave, a pair of direction-constrained waves with static labels $(\mathbf{k}, -\mathbf{k})$ is created:

$$\kappa = 0 \;\longrightarrow\; \kappa = 2:\quad (0, 0, 0, 0) \;\longrightarrow\; (k_x, 0, 0, k_t) + (-k_x, 0, 0, -k_t)$$

By conservation law (I2), the dynamic state $\Delta\mathbf{k}_\gamma$ of the spherical wave is distributed to the created pair of direction-constrained waves.

**Proposition 4.4**: Pair creation is the reverse process of pair annihilation; a pair of $\kappa = 2$ waves is created from the energy of a spherical wave. The static labels of the created direction-constrained waves are restricted to antisymmetric $(\mathbf{k}, -\mathbf{k})$.

---

## S5 Strong Interaction and Confinement

### S5.1 $Q$-Axis Transition and the 3-Bit Structure of Color Charge

The gluon $g = (0, 0, 0, 0, 0, Q_i \to Q_j)$ is defined as a transition on the $Q$ axis. By condition (I1), the displacement transfer acts only in the $Q$-axis direction:

$$\delta\mathbf{k}_g: \quad Q \;\longrightarrow\; Q'$$

The spatial, temporal, and $R$ components of the direction-constrained wave do not change. Only the discrete label $Q$ changes.

In Table A of [4], $Q$ is encoded as 3 bits $(c_1, c_2, c_3)$:

- $c_2 c_3$: Color charge ($00$ = no color, $01, 10, 11$ = three colors)
- $c_1$: Weak isospin

The $Q$ transition of gluons is a change of the $c_2 c_3$ bits, corresponding to color charge exchange. The three non-zero values of $c_2 c_3$ ($01, 10, 11$) yield $3 \times 3 - 1 = 8$ independent gluons.

### S5.2 Absorption: $q + g \to q'$

When a quark-type direction-constrained wave with color charge $(c_2 c_3) = (01)$ (red) absorbs gluon $g_{01 \to 10}$:

$$(k_x, 0, 0, k_t, 0, Q_{01}) + g_{01 \to 10} \;\longrightarrow\; (k_x, 0, 0, k_t, 0, Q_{10})$$

The gluon ceases to exist as an independent wave, and the color charge of the quark changes from red ($01$) to green ($10$). The spatial, temporal, and $R$ components are conserved.

### S5.3 Emission: $q \to q' + g$

By Proposition 3.1 (reversibility), emission is defined as the reverse process of absorption:

$$(k_x, 0, 0, k_t, 0, Q_{10}) \;\longrightarrow\; (k_x, 0, 0, k_t, 0, Q_{01}) + g_{10 \to 01}$$

The color charge of the quark changes, and a gluon carrying the difference is newly created.

**Proposition 5.1**: Interactions with gluons (both absorption and emission) change only the $Q$ label (color charge) of the direction-constrained wave, preserving the spatial, temporal, $R$ components, and dynamic state.

### S5.4 Pair Annihilation and Pair Creation: $q + \bar{q} \leftrightarrow g(s)$

In a quark-antiquark pair, the spatial axes ($x, y, z$) have inverted signs, while the $Q$ value (color charge) and the sign of the $t$ axis remain unchanged ([4] S9.9). Therefore, the antiquark carries the same color charge $(c_2 c_3)$ as the corresponding quark, and the bitwise composition (XOR) of color charges reduces to zero:

$$(c_2 c_3) \oplus (c_2 c_3) = (00)$$

Pair annihilation occurs, and the $\kappa = 2$ wave pair is converted to $\kappa = 0$ gluons. As the reverse process (pair creation), a quark-antiquark pair with the same color $(c_2 c_3)$ is generated from gluon energy:

$$q(c_2 c_3) + \bar{q}(c_2 c_3) \;\longleftrightarrow\; g(s)$$

**Proposition 5.2**: Quark-antiquark pairs carrying the same color charge can undergo pair annihilation and pair creation via gluons, as the XOR composition of their color charges reduces to zero.

### S5.5 Confinement Condition

Within the framework of subspace central projection introduced in [3] S8, to form an independent subjective space, a collection of waves must be **color-neutral**.

**Definition 5.1** (Color neutrality condition): A set $\{\psi_i\}$ of direction-constrained waves is color-neutral if the composition of $c_2 c_3$ bits reduces to $(0, 0)$. Here, bit composition is defined as bitwise exclusive OR (XOR, addition over $\mathbb{F}_2$):

$$(a_1, a_2) \oplus (b_1, b_2) = (a_1 \oplus b_1,\; a_2 \oplus b_2) \quad (\oplus: \text{mod}\; 2 \text{ addition})$$

Realizations of color neutrality:

| Configuration | Color combination | Name |
|--------------|------------------|------|
| 3 waves | $(01) \oplus (10) \oplus (11) = (00)$ | Baryon type |
| 2 waves | $(c_2 c_3) \oplus (c_2 c_3) = (00)$ | Meson type |
| 1 wave | $(c_2 c_3) = (00)$ | Lepton type |

**Remark 5.1**: In the meson type, the quark and antiquark carry the same color $(c_2 c_3)$ (S5.4), so color neutrality is achieved by the XOR composition of two identical colors. In standard QCD, color neutrality is described as the annihilation of color and anti-color ($\mathrm{SU}(3)$: $\mathbf{3} \otimes \bar{\mathbf{3}} \supset \mathbf{1}$), whereas in this framework, XOR algebra over $\mathbb{F}_2^2$ provides the corresponding operation. The algebraic structures differ, but the result -- color neutrality -- is identical.

**Proposition 5.3**: A set of direction-constrained waves that does not satisfy the color neutrality condition cannot form an independent subjective space and is confined within a color-neutral set.

### S5.6 Asymptotic Freedom

Inside a color-neutral set (inside the shared subjective space), the curvature of the subspace central projection is small, and direction-constrained waves move almost freely. As one approaches the boundary of the set, the projection curvature increases, making it impossible to cross the boundary.

**Proposition 5.4**: The interior of a color-neutral set has small projection curvature (asymptotic freedom), while the curvature diverges at the boundary (confinement).

---

## S6 Weak Interaction and Flavor Conversion

### S6.1 $t$-Axis Displacement

$W^{\pm} = (0, 0, 0, \pm t, 0, 0)$ is a spherical wave with a non-zero component on the $t$ axis. By condition (I1):

$$\delta\mathbf{k}_W: \quad k_t \;\longrightarrow\; k_t' = k_t + \delta k_t$$

The change in the $t$-axis component corresponds to the conversion between up-type ($k_t > 0$) and down-type ($k_t < 0$) in Table A of [4]. Simultaneously, the $c_1$ bit (weak isospin) of the $Q$ axis also changes. $c_1$ and $\mathrm{sign}(k_t)$ are independent degrees of freedom, and $W^{\pm}$ mediates both. The $c_2 c_3$ bits (color charge) of the $Q$ axis are independent degrees of freedom orthogonal to the $t$ axis and $c_1$, and are conserved under interactions with $W^{\pm}$.

The $t$-axis displacement by $W^{\pm}$ acts as a sign inversion when $k_t$ is non-zero in both the initial and final states (quark system, S6.2--6.3). When the transition involves $k_t = 0$ (lepton system), axis reconfiguration accompanying $\kappa = 2$ preservation occurs (S6.7).

### S6.2 Absorption: $f + W \to f'$

Transformation of a direction-constrained wave by $W^+$ absorption:

$$\underbrace{(k_x, 0, 0, -k_t, 0, Q_d)}_{\text{down-type}} + \underbrace{(0, 0, 0, +, 0, 0)}_{W^+} \;\longrightarrow\; \underbrace{(k_x, 0, 0, +k_t, 0, Q_u)}_{\text{up-type}}$$

$W^+$ ceases to exist as an independent wave, the sign of the $t$ axis of the direction-constrained wave is inverted ($-k_t \to +k_t$), and the $c_1$ bit transitions $0 \to 1$ ($Q_d \to Q_u$). Color charge $c_2 c_3$ is conserved.

### S6.3 Emission: $f \to f' + W$

By Proposition 3.1 (reversibility), emission is defined as the reverse process of absorption:

$$\underbrace{(k_x, 0, 0, +k_t, 0, Q_u)}_{\text{up-type}} \;\longrightarrow\; \underbrace{(k_x, 0, 0, -k_t, 0, Q_d)}_{\text{down-type}} + \underbrace{(0, 0, 0, +, 0, 0)}_{W^+}$$

The flavor of the direction-constrained wave changes from up-type to down-type ($t$-axis sign inversion, $c_1: 1 \to 0$), and $W^+$ is newly created as an independent wave. By conservation law (I2), the change in the $t$-axis component is emitted as $W^+$.

**Independence of $Q$-axis color charge and the $t$ axis**: In Table A of [4], the $Q$ axis is encoded as 3 bits $(c_1, c_2, c_3)$, where $c_1$ represents weak isospin and $c_2 c_3$ represents color charge (S5.1). $c_1$ and the sign of the $t$ axis $\mathrm{sign}(k_t)$ are independent degrees of freedom, and $W^{\pm}$ changes both simultaneously (S6.5). On the other hand, $c_2 c_3$ (color charge) is an independent degree of freedom orthogonal to the $t$ axis and $c_1$, and (I1) requires displacement transfer only along the non-zero axis of the spherical wave (the $t$ axis), so color charge $c_2 c_3$ is conserved under interactions with $W^{\pm}$.

**Proposition 6.1**: Interactions with $W^{\pm}$ (both absorption and emission) change the $t$-axis component (flavor) and the $c_1$ bit (weak isospin) of the $Q$ axis of the direction-constrained wave. The $c_2 c_3$ bits (color charge) of the $Q$ axis are independent degrees of freedom orthogonal to the $t$ axis and $c_1$, and are conserved under interactions with $W^{\pm}$.

### S6.4 $W$ Decay and Inverse Decay

The decay of $W^{\pm}$ itself is a special case of pair creation (Type D):

$$W^+ \;\longrightarrow\; f_u + \bar{f}_d$$

$W^+$ is annihilated, and an up-type direction-constrained wave $f_u$ and a down-type anti-direction-constrained wave $\bar{f}_d$ are created. By conservation law (I2), the $t$-axis component and dynamic state of $W^+$ are distributed to the created wave pair.

Reverse process (inverse decay):

$$f_u + \bar{f}_d \;\longrightarrow\; W^+$$

A direction-constrained wave and an anti-direction-constrained wave are annihilated, and $W^+$ is created. This corresponds to the $t$-axis version of pair annihilation (Type C).

**Proposition 6.2**: $W$ decay is the $t$-axis version of Type D, and inverse decay is the $t$-axis version of Type C. Both are causal processes with $k_t \neq 0$, and decay and inverse decay are temporally distinguishable.

### S6.5 Chirality Selection Rule

For the chirality $\chi = \mathrm{sign}(k_{s_0} \cdot k_t)$ ($k_{s_0}$: non-zero spatial component) defined in Definition 6.1 of [3], we derive the coupling condition of $W^{\pm}$. In Table A of [4], the $c_1$ bit of the $Q$ axis encodes weak isospin and is an independent degree of freedom from the $t$ axis. $W^{\pm}$ mediates the change of the $t$ axis and the transition of $c_1$ respectively:

- $W^+$: $c_1 = 0 \to 1$ (down-type $\to$ up-type)
- $W^-$: $c_1 = 1 \to 0$ (up-type $\to$ down-type)

Using the isospin sign $(-1)^{c_1}$ of $c_1$, chirality is reformulated as:

$$\chi = \mathrm{sign}((-1)^{c_1} \cdot k_t)$$

**Case of $W^+$ absorption**: $W^+ + \psi_d \to \psi_u$. The initial state is down-type ($c_1 = 0$, $k_t < 0$):

$$\chi(\psi_d) = \mathrm{sign}((-1)^0 \cdot (-|k_t|)) = \mathrm{sign}((+1) \cdot (-)) = -1$$

**Case of $W^-$ absorption**: $W^- + \psi_u \to \psi_d$. The initial state is up-type ($c_1 = 1$, $k_t > 0$):

$$\chi(\psi_u) = \mathrm{sign}((-1)^1 \cdot (+|k_t|)) = \mathrm{sign}((-1) \cdot (+)) = -1$$

In both cases, $\chi = -1$ (left-handed) is derived without additional assumptions. Because $c_1$ and $k_t$ are independent degrees of freedom that distinguish $W^+$ and $W^-$, there is no need to assume sign reversal of $k_{s_0}$.

**Proposition 6.3**: The coupling condition of $W^{\pm}$ is determined by the sign of the product of the isospin sign $(-1)^{c_1}$ of the initial state and $k_t$. For both $W^+$ ($c_1: 0 \to 1$) and $W^-$ ($c_1: 1 \to 0$), $\chi = -1$ (left-handed) holds, and right-handed ($\chi = +1$) direction-constrained waves do not couple to $W^{\pm}$.

### S6.6 CP Violation

As shown in Proposition 6.4 of [3], CP symmetry can only be violated in interactions with $k_t \neq 0$. $W^{\pm}$ has $k_t = \pm 1 \neq 0$ and is the only one among the four types of spherical waves with $k_t \neq 0$.

**Proposition 6.4**: Only $W^{\pm}$, which has a non-zero component on the $t$ axis, mediates CP symmetry violation. This is a direct consequence of $k_t \neq 0$ introducing temporal asymmetry.

### S6.7 Axis Reconfiguration in the Lepton System

The $W^{\pm}$ absorption and emission described in S6.2--6.3 treated the sign inversion of $k_t$ (where $k_t \neq 0$ in both initial and final states). This is the typical process for the quark system ($c_2 c_3 \neq 00$). In the lepton system ($c_2 c_3 = 00$), the axis configurations of charged leptons and neutrinos differ qualitatively, so the action of $W^{\pm}$ involves axis reconfiguration.

In Table A of [4], fermions ($\kappa = 2$) with $c_2 c_3 = 00$ (no color charge) have the following two types of axis configurations:

| Type | Non-zero axes | $k_t$ | $c_1$ | $Q$ | Example |
|------|--------------|-------|-------|-----|---------|
| Charged lepton | 1 spatial axis + $t$ | $\neq 0$ | 1 | 4 | $e^-: (k_x, 0, 0, k_t, 0, 4)$ |
| Neutrino | 2 spatial axes | $= 0$ | 0 | 0 | $\nu: (k_{s_1}, 0, k_{s_2}, 0, 0, 0)$ |

Charged lepton $\to$ neutrino conversion by $W^+$ absorption:

$$\underbrace{(k_x, 0, 0, k_t, 0, 4)}_{\substack{\text{charged lepton} \\ \kappa = 2:\; (s, t),\; c_1 = 1}} + \underbrace{(0, 0, 0, +, 0, 0)}_{W^+} \;\longrightarrow\; \underbrace{(k_x, 0, k_{s_2}, 0, 0, 0)}_{\substack{\text{neutrino} \\ \kappa = 2:\; (s, s),\; c_1 = 0}}$$

This process cannot be fully described by (I1) alone. (I1) requires displacement along the non-zero axis of $W^+$ (the $t$ axis) and permits the transition $k_t \to 0$. However, the transition to $k_t = 0$ reduces the number of non-zero axes from $\kappa = 2$ to $\kappa = 1$, so the stability condition for $\kappa = 2$ (Theorem 5.1 of [3]) demands compensatory activation of a spatial axis.

That is, the action of $W^{\pm}$ on the lepton system is determined as the simultaneous satisfaction of two conditions:

1. **(I1)**: Displacement along the $t$ axis ($k_t \to 0$ or $0 \to k_t$)
2. **$\kappa = 2$ preservation**: The fermion condition maintains the number of non-zero axes at $\kappa = 2$

Reverse process (neutrino $\to$ charged lepton conversion):

$$\underbrace{(k_x, 0, k_{s_2}, 0, 0, 0)}_{\substack{\text{neutrino} \\ \kappa = 2:\; (s, s),\; c_1 = 0}} + \underbrace{(0, 0, 0, -, 0, 0)}_{W^-} \;\longrightarrow\; \underbrace{(k_x, 0, 0, k_t, 0, 4)}_{\substack{\text{charged lepton} \\ \kappa = 2:\; (s, t),\; c_1 = 1}}$$

$W^-$ supplies a $t$-axis component ($0 \to k_t$), and $\kappa = 2$ preservation requires that one of the spatial axes ($k_{s_2}$) becomes deactivated.

**Proposition 6.5**: The action of $W^{\pm}$ on the lepton system is described as the simultaneous satisfaction of $t$-axis displacement (I1) and $\kappa = 2$ preservation. The conversion between charged leptons (1 spatial axis + $t$) and neutrinos (2 spatial axes) preserves $\kappa = 2$ through the interchange of $t$-axis and spatial-axis components.

**Remark 6.1**: The choice of which spatial axis is activated (or deactivated) in the axis reconfiguration is not uniquely determined by $\kappa = 2$ preservation alone. For example, from $e^-: (k_x, 0, 0, k_t, 0, 4)$, a neutrino generated by $W^+$ absorption could be either $(k_x, k_y, 0, 0, 0, 0)$ or $(k_x, 0, k_z, 0, 0, 0)$, both satisfying $\kappa = 2$. This degree of freedom may constitute the geometric origin of neutrino mixing (a structure corresponding to the PMNS matrix), but its elucidation is outside the scope of this paper.

---

## S7 Gravitational Interaction

### S7.1 Displacement Transfer by the Graviton

The graviton $G = (0, 0, 0, \pm t, \pm R, 0)$ has non-zero components on both the $t$ and $R$ axes. By condition (I1), the displacement transfer acts in both the $t$ and $R$ directions, with the $R$ direction subject to the discretization condition (I4) (S2.5):

$$\delta(\Delta\mathbf{k})_G = (\delta(\Delta k)_t,\; \delta(\Delta k)_R)$$

This means simultaneous changes in energy ($t$ component) and spatial configuration ($R$ component). As a candidate consequence of (I4), only discrete allowed values may exist for the energy-momentum carried by the graviton along the $R$ direction.

### S7.2 Absorption: $f + G \to f'$

When a direction-constrained wave absorbs a graviton:

$$(\mathbf{k},\; \Delta\mathbf{k}) + (\mathbf{k}_G,\; \Delta\mathbf{k}_G) \;\longrightarrow\; (\mathbf{k},\; \Delta\mathbf{k} + \delta(\Delta\mathbf{k})_G)$$

The graviton ceases to exist as an independent wave, and its energy-momentum is transferred to the dynamic state of the direction-constrained wave. The simultaneous change in both the $t$ and $R$ directions distinguishes this from the photon ($R$ only).

### S7.3 Emission: $f \to f' + G$

By Proposition 3.1 (reversibility), emission is defined as the reverse process of absorption:

$$(\mathbf{k},\; \Delta\mathbf{k}) \;\longrightarrow\; (\mathbf{k},\; \Delta\mathbf{k} - \delta(\Delta\mathbf{k})_G) + (\mathbf{k}_G,\; \delta(\Delta\mathbf{k})_G)$$

A portion of the dynamic state of the direction-constrained wave separates, and a new graviton is created as an independent wave.

**Proposition 7.1**: Interactions with the graviton (both absorption and emission) change only the dynamic state $\Delta\mathbf{k}$ and do not act on the static label $\mathbf{k}$.

### S7.4 Universal Coupling

We restate the description of acceleration introduced in [3] S8. The second-order difference of the subspace central projection determines the acceleration:

$$\Delta^2\mathbf{k}(\mathbf{k}_t) = D^2\pi(\mathbf{k}_t) \cdot \Delta\mathbf{k}_t$$

Graviton exchange is the discrete realization of this $\Delta^2\mathbf{k}$. Since the graviton acts on the dynamic state $\Delta\mathbf{k}$:

- It is independent of the static label $\mathbf{k}$ (regardless of wave type)
- It acts on all waves with $\Delta\mathbf{k} \neq 0$ (all waves carrying energy)

**Proposition 7.2**: All shape-invariant waves with $\Delta\mathbf{k} \neq 0$ couple to the graviton (universal coupling).

**Remark 7.1**: Pair annihilation ($f + \bar{f} \to G$) and pair creation ($G \to f + \bar{f}$) via gravity are permitted in principle under conservation law (I2). However, since the graviton has $k_t \neq 0$, these are causal processes, and since coupling is mediated only through the dynamic state, the coupling strength is remarkably weak compared to other forces. Quantitative evaluation of these processes is outside the scope of this paper.

---

## S8 Origin of Mass: Coupling with Standing Waves

### S8.1 Geometric Structure of Coupling

The standing wave mode ($\kappa = 1$, [3] S5) fills the entire subjective space $S^3$ as a fundamental mode along one spatial axis. Direction-constrained waves ($\kappa = 2$) propagate through this standing wave background.

**Remark 8.1**: Higgs coupling is outside the scope of the fundamental vertex classification of S3 (Types A--D, which involve spherical waves with $\kappa = 0$). The standing wave ($\kappa = 1$) is not a spherical wave ($\kappa = 0$), and its coupling with direction-constrained waves is described independently in this section as a modification of the speed ratio, not through fundamental vertices.

In Table A of [4], the non-zero axis of the Higgs boson is assigned to the $z$ axis ($(0, 0, +, 0, 0, 0)$). The $z$ axis corresponds to the 3rd generation in the generation labels (generation 1: $x$, generation 2: $y$, generation 3: $z$, [4] S9.9). The conclusions of S8 of this paper (direction of the generation mass hierarchy) depend on this correspondence between generation numbers and spatial axes, and presuppose that the axis assignment of $H^0$ in Table A of [4] and the generation label assignment are mutually consistent. With this assignment, the standing wave is expressed as $\mathbf{k}_H = (0, 0, k_H, 0)$.

A 3rd-generation direction-constrained wave $\mathbf{k} = (0, 0, k_z, k_t)$ and the standing wave $\mathbf{k}_H = (0, 0, k_H, 0)$ directly share the spatial $z$ axis. On this shared axis, the phase structures of the waves overlap, producing the strongest coupling. 1st-generation ($\mathbf{k} = (k_x, 0, 0, k_t)$) and 2nd-generation ($\mathbf{k} = (0, k_y, 0, k_t)$) direction-constrained waves do not share the $z$ axis with the standing wave, so they lack direct coupling and only indirect coupling through the global topological structure of $S^3$ is possible.

### S8.2 Modification of the Speed Ratio

The propagation speed of a direction-constrained wave without coupling is $v = k_{s_0} / k_t$ ($k_{s_0}$: non-zero spatial component). Coupling with the standing wave affects the phase structure along the shared spatial axis direction, modifying the effective spatial component:

$$k_{s_0} \;\longrightarrow\; k_{s_0}^{\text{eff}} = k_{s_0} - \delta k_H$$

where $\delta k_H > 0$ is proportional to the coupling strength with the standing wave. The sign of $\delta k_H$ (the direction of reducing the effective spatial component) is derived from the properties of the standing wave. The standing wave has $\kappa = 1$, $k_t = 0$, and possesses no temporal phase propagation. A wave without propagation phase cannot advance in phase with a propagating wave, and coupling can only act in the cancellation direction of the spatial phase. Furthermore, the standing wave fills the entire $S^3$ as a fundamental mode oscillation, so direction-constrained waves cannot propagate while avoiding the standing wave background. Therefore, for direction-constrained waves sharing the $z$ axis, coupling in the cancellation direction is the unique and inevitable consequence. The propagation speed after coupling is:

$$v^{\text{eff}} = \frac{k_{s_0}^{\text{eff}}}{k_t} = \frac{k_{s_0} - \delta k_H}{k_t} < \frac{k_{s_0}}{k_t} = v$$

The speed decreases -- this is the **geometric origin of mass**.

**Proposition 8.1**: Coupling with a standing wave ($\kappa = 1$) reduces the effective spatial component of a direction-constrained wave, limiting the propagation speed to $v < c$. This speed limitation is observed as mass.

### S8.3 Waves That Do Not Couple

Spherical waves ($\kappa_s = 0$) have no spatial components. The standing wave has phase structure along one spatial axis, but spherical waves have no spatial axis to share.

**Proposition 8.2**: Spherical waves with $\kappa_s = 0$ have no coupling axis with the standing wave, yielding $\delta k_H = 0$. The speed remains $v = c$ without modification.

This is consistent with the gauge bosons with $\kappa_s = 0$ ($\gamma$, $g$) having zero mass in Table A of [4]. $W^{\pm}$ and $Z^0$ have $\kappa_s = 0$ and therefore lack direct coupling through spatial axes as described by Proposition 8.2, but they may acquire mass through an indirect standing wave coupling mechanism via the $t$ or $R$ axis. Proposition 8.2 addresses only coupling through direct spatial axis sharing; the details of indirect coupling mechanisms via the $t$ or $R$ axis are outside the scope of this paper.

### S8.4 Mass Hierarchy

The coupling strength $\delta k_H$ depends on the overlap integral along the spatial axis between the direction-constrained wave and the standing wave. With the standing wave's non-zero axis assigned to $z$ (the 3rd-generation axis) in S8.1, the ordering of coupling strengths among generations is structurally determined:

| Generation | Direction-constrained wave | Axis sharing with standing wave | Coupling strength $\delta k_H$ |
|-----------|---------------------------|-------------------------------|-------------------------------|
| 3rd | $(0, 0, k_z, k_t)$ | Direct sharing of $z$ axis | Maximum |
| 1st / 2nd | $(k_x, 0, 0, k_t)$ / $(0, k_y, 0, k_t)$ | No shared axis (indirect coupling) | Small (mechanism unresolved) |

The 3rd generation directly shares the spatial axis with the standing wave, yielding the maximum $\delta k_H$ and acquiring the largest mass. Empirically, the 3rd generation is the heaviest ($m_t \gg m_c \gg m_u$), which is consistent with this structural prediction.

The 1st and 2nd generations do not share the $z$ axis and therefore have only indirect coupling through the global topological structure of $S^3$. However, since the $x$ and $y$ axes are locally equivalent under $\mathrm{SO}(3)$ symmetry, the difference $m_2 \neq m_1$ cannot be derived within the framework of this paper. The origin of the mass difference between the 1st and 2nd generations depends on an additional mechanism that breaks this $\mathrm{SO}(3)$ symmetry (the detailed structure of indirect coupling), but its elucidation is outside the scope of this paper (see S1).

What this paper structurally claims is only that "$m_3 \gg m_{1,2}$" (the 3rd generation is overwhelmingly heavier than the 1st and 2nd generations).

---

## S9 Wave Packet Deformation and Retroactive Construction of Causality

### S9.1 Conservation-Law Description of Wave Packet Deformation

In the interactions of S4--S7, when a direction-constrained wave receives energy, under the conserved quantity $\int |\psi|^2 \, dV = E_0$ ((C1) from [1] S5, [2] S5), the amplitude $A$ and characteristic width $w$ of the wave packet satisfy:

$$A^2 \cdot w^d = E_0 = \text{const} \quad (d = 3: \text{spatial dimension of } S^3)$$

A wave packet that receives energy:

$$\Delta E > 0 \;\Rightarrow\; A \uparrow,\; w \downarrow \quad (\text{amplitude increase, width decrease})$$

$$\Delta E < 0 \;\Rightarrow\; A \downarrow,\; w \uparrow \quad (\text{amplitude decrease, width increase})$$

**Proposition 9.1**: Under the conserved quantity $\int |\psi|^2 \, dV$, energy transfer deforms the amplitude and width of a wave packet inversely. Energy increase is accompanied by localization (width decrease) of the wave packet.

### S9.2 Wave Packet Deformation of Spherical Waves

Spherical waves ($\kappa = 0$) also couple to the graviton (Proposition 7.2: coupling occurs if $\Delta\mathbf{k} \neq 0$). When a spherical wave receives energy:

- The amplitude of the spherical wave increases and its spatial spread decreases
- In the limit, the spherical wave becomes point-like localized

This corresponds to the geometric description of gravitational lensing. When a spherical wave (photon) passes through the gravitational field of a massive body, the wave packet is deformed through interaction with gravitons, bending the propagation path.

### S9.3 Dynamical Description of "Observation"

In conventional quantum mechanics, "wave packet collapse through observation" is introduced axiomatically. In this framework:

1. A direction-constrained wave ($\kappa = 2$) interacts with a spherical wave ($\kappa = 0$)
2. Energy transfer increases the amplitude of the wave packet
3. By conservation law (I3) and Proposition 9.1, the width decreases (localization)
4. A sufficiently localized wave packet is described as "detected as a particle"

"Observation" is not a special operation but merely one example of wave packet deformation through interaction. No external "observer" is required, and wave packet deformation is a deterministic process following conservation laws.

**Proposition 9.2**: Wave packet collapse is a conservation-law-governed wave packet deformation through interaction, and there is no need to postulate a special process called "observation."

### S9.4 Retroactive Construction of Causality

As shown in Proposition 4.2 of [3], spherical waves with $k_t = 0$ are time-reversal symmetric and carry no intrinsic causal ordering. The "deterministic" of Proposition 9.2 and the "retroactive" below are not contradictory. They belong to different layers:

- **Deterministic**: Conservation laws (I2)(I3) hold strictly, and the result of wave packet deformation is uniquely determined
- **Retroactive**: Due to $k_t = 0$, causal ordering (the temporal order of "cause $\to$ effect") is not intrinsic, so which event is the "cause" and which is the "effect" is indeterminate

That is, the structure "conservation laws are strictly followed, but causal ordering is absent" arises naturally from the time-reversal symmetry of $k_t = 0$. Combining this property with wave packet deformation:

**(R1) Retroactive determination of path**: A spherical wave ($k_t = 0$) propagates isotropically in all directions (Proposition 4.1 of [3]). Viewed from the result of the interaction, the spherical wave **appears to have** propagated along the shortest geodesic. Since causal ordering does not exist, the result retroactively determines the path.

**(R2) Retroactive determination of localization**: Wave packet deformation is a deterministic process following conservation laws (Proposition 9.2). However, since $k_t = 0$ eliminates causal ordering, "the interaction occurred and therefore localization happened" and "localization happened and therefore the interaction appears to have occurred" are indistinguishable.

**(R3) Principle of least action**: As a synthesis of (R1) and (R2), the variational principle of natural law (the principle of least action) is not an axiom but is derived as a consequence of the causal symmetry of interactions mediated by spherical waves with $k_t = 0$.

**Proposition 9.3**: In interactions mediated by spherical waves with $k_t = 0$, (a) the shortest propagation path, (b) wave packet localization, and (c) the variational principle are all retroactively constructed from the absence of causal ordering.

### S9.5 The Case of $k_t \neq 0$

Interactions mediated by $W^{\pm}$ ($k_t \neq 0$) and $G$ ($k_t \neq 0$) are not causally symmetric. (R1)--(R3) do not apply:

- The path is causally determined (emission first, absorption later)
- Time-reversal symmetry is broken
- CP non-conservation becomes possible (Proposition 6.4, $W^{\pm}$ only)

| $k_t$ of spherical wave | Causal ordering | Path selection | CP | Corresponding force |
|------------------------|----------------|---------------|-----|-------------------|
| $k_t = 0$ | None | Retroactive (geodesic) | Conserved | Electromagnetism, strong force |
| $k_t \neq 0$ | Present | Causal | May be violated | Weak force, gravity |

Combined with Remark 3.1 of S3.2: when $k_t = 0$, absorption and emission are temporally indistinguishable (the reversibility of S3 is completely symmetric), whereas when $k_t \neq 0$, even though they are reversible, they are temporally distinguishable.

---

## S10 Summary

The main results of this paper are summarized:

1. **Unified interaction mechanism**: All forces are described as displacement transfer (I1)--(I3) along non-zero axes of spherical waves, and the type of force reduces solely to the difference in axes (S2).

2. **Fundamental vertices and reversibility**: Interactions are classified into four fundamental vertex types (absorption, emission, pair annihilation, pair creation), and under conservation laws, absorption $\Leftrightarrow$ emission and pair annihilation $\Leftrightarrow$ pair creation are reversible (S3).

3. **Electromagnetic interaction**: The photon acts on the $R$ axis and changes only the dynamic state $\Delta\mathbf{k}$. All processes -- absorption, emission, exchange, pair annihilation, pair creation -- are systematically described. The relative sign of $t$-axis components determines attraction/repulsion (S4).

4. **Strong interaction**: Gluons mediate $Q$-axis transitions and exchange color charge. Absorption, emission, pair annihilation, pair creation are systematically described. The color neutrality condition provides the origin of confinement, and asymptotic freedom holds inside the set (S5).

5. **Weak interaction**: $W^{\pm}$ mediates $t$-axis displacement and $c_1$ (weak isospin) transition simultaneously. In the quark system, this acts as a sign inversion of $k_t$ (and $c_1: 0 \leftrightarrow 1$); in the lepton system, the simultaneous satisfaction of $t$-axis displacement and $\kappa = 2$ preservation produces axis reconfiguration between charged leptons and neutrinos. In addition to absorption and emission, $W$ decay and inverse decay are unified as Type D/C. $k_t \neq 0$ provides the origin of CP violation and the left-handed chirality selection rule (S6).

6. **Gravitational interaction**: The graviton acts on both $(t, R)$ axes and couples to all waves with $\Delta\mathbf{k} \neq 0$ (universal coupling). Absorption and emission are systematically described (S7).

7. **Origin of mass**: Coupling with the standing wave ($\kappa = 1$, $z$ axis) reduces the speed ratio $|k_{s_0} / k_t|$ of direction-constrained waves, producing mass as $v < c$. The 3rd generation, which directly shares the $z$ axis, has the strongest coupling strength, providing the origin of the inter-generation mass hierarchy $m_3 \gg m_{1,2}$. Spherical waves with $\kappa_s = 0$ lack a coupling axis and have zero mass (S8).

8. **Wave packet deformation and causal structure**: Wave packet collapse is a conservation-law-governed deformation and "observation" is unnecessary. Interactions mediated by spherical waves with $k_t = 0$ are causally symmetric, and the path, localization, and variational principle are retroactively constructed. The case $k_t \neq 0$ has causal temporal ordering (S9).

The correspondence with Table A of [4] is shown below:

| Interaction | Spherical wave | Non-zero axis | Effect on direction-constrained wave | $k_t$ | Causality | Fundamental vertices |
|------------|---------------|--------------|-------------------------------------|-------|-----------|---------------------|
| Higgs coupling | Standing wave $H^0$ | $z$ axis | Decrease in $v$ (mass) | 0 | -- | -- |
| Electromagnetic | $\gamma$ | $R$ | $\Delta\mathbf{k}$ change (scattering) | 0 | Retroactive | A, B, C, D |
| Strong | $g$ | $Q$ | $Q$ change (color exchange) | 0 | Retroactive | A, B, C, D |
| Weak | $W^{\pm}$ | $t$ | $k_t$ change (flavor conversion) | $\neq 0$ | Causal | A, B, C, D |
| Gravity | $G$ | $t, R$ | $\Delta\mathbf{k}$ change (acceleration) | $\neq 0$ | Causal | A, B |

All claims in this paper are stated as geometric facts, and the physical correspondences are presented as "connectable interpretations."

---

## References

[1] Noriaki Kihara, "Shape-Invariant Standing Waves on Closed Spheres -- Stability per Dimension and the Geometric Necessity of 3+1 Spacetime," 2026.

[2] Noriaki Kihara, "Shape-Invariant Traveling Waves on 4-Dimensional Integer Lattices -- Existence Conditions, Conservation Structure, and Self-Consistent Loop," 2026.

[3] Noriaki Kihara, "Four Modes of Shape-Invariant Waves -- Spin, Chirality, and Statistics Determined by Wave Vector Structure," 2026.

[4] Noriaki Kihara, "Set-Theoretic Structure and Combinatorial Properties of the 6-Dimensional Hypercube," Zenodo, 2025. DOI: 10.5281/zenodo.19688521.

[5] Noriaki Kihara, "Three-Layer Model of the Universe via Central Projection," Zenodo, 2025. DOI: 10.5281/zenodo.15083729.
