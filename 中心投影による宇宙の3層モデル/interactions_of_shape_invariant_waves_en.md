# Interactions of Shape-Invariant Waves: Axis-Directional Displacement Transfer, Wave Packet Deformation, and Retroactive Construction of Causality

**Author**: Noriaki Kihara  
**Affiliation**: WF System Co., Ltd.  
**ORCID**: 0009-0004-6753-4020  
**Version**: v1.0  
**Date**: April 2026  
**DOI**: 10.5281/zenodo.19709800

---

## Abstract

Building on the three-mode classification of shape-invariant waves introduced in the companion paper [3] (direction-constrained waves $\kappa = 2$, spherical waves $\kappa = 0$, standing waves $\kappa = 1$), we provide a unified description of interactions between different modes. The fundamental mechanism of interaction is defined as **axis-directional displacement transfer** -- the change of a direction-constrained wave's state vector $(\mathbf{k}, \Delta\mathbf{k})$ along the direction of the non-zero axis carried by a spherical wave -- and the four forces (electromagnetism, weak force, strong force, gravity) are unified as applications of this single operation to different axes. Coupling with a standing wave ($\kappa = 1$) modifies the propagation speed ratio $|k_x / k_t|$ of a direction-constrained wave, providing a geometric origin for mass. Under the conserved quantity $\int |\psi|^2 \, dV$, energy transfer deforms the amplitude and width of a wave packet inversely, and this deformation gives a dynamical description of "wave packet collapse." Interactions mediated by spherical waves with $k_t = 0$ carry no intrinsic causal ordering, and the shortest geodesic path and the localization of the wave packet are retroactively constructed from the interaction outcome.

**Keywords**: shape-invariant wave, axis-directional displacement transfer, wave packet deformation, origin of mass, retroactive construction of causality, conservation law

---

## S1 Introduction and Scope of Claims

This paper provides a unified description of interactions between waves of different modes within the framework of shape-invariant waves constructed in the companion papers [1]--[3], and as consequences derives the origin of mass, the classification of forces, wave packet deformation, and causal structure geometrically.

The claims of this paper are limited to the following five points:

1. All interactions are described uniformly as displacement transfer along the non-zero axes of spherical waves. The difference among types of force reduces solely to the difference in the axis along which the transfer occurs.
2. Coupling with a standing wave ($\kappa = 1$) modifies the speed ratio of a direction-constrained wave, providing a geometric origin for mass.
3. Under the conserved quantity $\int |\psi|^2 \, dV$, energy transfer deterministically deforms the shape of a wave packet, and this deformation corresponds to "wave packet collapse" in observation.
4. Interactions mediated by spherical waves with $k_t = 0$ carry no intrinsic causal ordering, and the path and localization are retroactively constructed from the outcome.
5. The sign of the interaction (attraction / repulsion) between direction-constrained waves is determined geometrically by the relative sign of their $t$-axis components.

The following elements are **outside the scope of this paper**:

- Numerical derivation of coupling constants (specific values of $\alpha$, $G_F$, $\alpha_s$, $G$)
- Quantitative computation of scattering cross sections and decay widths
- Renormalization and scale dependence of running coupling constants
- Rigorous reconstruction of gauge groups in the continuum limit

**Note on physical nomenclature**: This paper uses physical names (photon, graviton, electron, flavor, color charge, etc.) assuming correspondence with Table A of [4], but this is for interpretive convenience; the geometric content of each proposition does not depend on the physical identification. All propositions are derived solely from the component structure of wave vectors ($\kappa$, $k_t$, type of axis).

---

## S2 Unified Definition of Interaction

### S2.1 Two-Layer Structure

In the tangent bundle $\mathbb{Z}^4 \times \mathbb{Z}^4$ introduced in [2] S8, the state of a direction-constrained wave ($\kappa = 2$) is:

$$(\mathbf{k},\; \Delta\mathbf{k}) \in \mathbb{Z}^4 \times \mathbb{Z}^4$$

- $\mathbf{k}$: **static label** (type of wave, corresponding to Table A of [4])
- $\Delta\mathbf{k}$: **dynamic state** (corresponding to energy and momentum)

### S2.2 Axis-Directional Displacement Transfer

**Definition 2.1** (Axis-directional displacement transfer): When a spherical wave ($\kappa = 0$) interacts with a direction-constrained wave ($\kappa = 2$), the state of the direction-constrained wave changes as follows:

$$(\mathbf{k},\; \Delta\mathbf{k}) \;\longrightarrow\; (\mathbf{k} + \delta\mathbf{k},\; \Delta\mathbf{k} + \delta(\Delta\mathbf{k}))$$

where the displacements $\delta\mathbf{k}$ and $\delta(\Delta\mathbf{k})$ satisfy the following conditions:

**(I1) Axis-directionality**: The non-zero components of $\delta\mathbf{k}$ and $\delta(\Delta\mathbf{k})$ exist only in the direction of the non-zero axis of the spherical wave in Table A of [4].

**(I2) Conservation law**: The total sum of the tangent bundle states of all waves involved in the interaction is conserved:

$$\sum_i (\mathbf{k}_i + \Delta\mathbf{k}_i) = \text{const}$$

**(I3) Wave packet conservation**: The conserved quantity $\int |\psi|^2 \, dV$ of each wave is individually conserved.

### S2.3 Static Label Change and Dynamic State Change

Axis-directional displacement transfer encompasses two qualitatively distinct types of change:

| Type of change | Target | Physical meaning | Example |
|-----------|------|----------|-----|
| $\delta\mathbf{k} \neq 0$ | Static label $\mathbf{k}$ | Conversion of wave type | Flavor transition, color charge exchange |
| $\delta(\Delta\mathbf{k}) \neq 0$ | Dynamic state $\Delta\mathbf{k}$ | Change in energy/momentum | Scattering, acceleration |

The type of spherical wave determines which (or both) of these is dominant.

### S2.4 Unification of the Four Forces

From the non-zero axis structure of spherical waves in Table A of [4], the four forces are described as applications of the same operation to different axes:

| Force | Spherical wave | Non-zero axis | Dominant change |
|----|-------|---------|----------|
| Electromagnetism | $\gamma\; (0,0,0,0,R,0)$ | $R$ | $\delta(\Delta\mathbf{k})$ (momentum change) |
| Weak force | $W^{\pm}\; (0,0,0,\pm t,0,0)$ | $t$ | $\delta\mathbf{k}$ (label change) |
| Strong force | $g\; (0,0,0,0,0,Q_i)$ ${}^{*}$ | $Q$ | $\delta\mathbf{k}$ (label change) |
| Gravity | $G\; (0,0,0,\pm t,\pm R,0)$ | $t, R$ | $\delta(\Delta\mathbf{k})$ (energy change) |

${}^{*}$ While $\gamma$, $W^{\pm}$, and $G$ are spherical waves with static non-zero components on the $t$ or $R$ axis, $g$ is a spherical wave carrying a discrete label value $Q_i$ on the $Q$ axis that mediates color charge exchange by transitioning to a wave with a different label value $Q_j$. That is, the "non-zero axis" of $g$ acts not as a fixed value but as a change in the $Q$ label (detailed in S6).

**Proposition 2.1**: All four forces are described as axis-directional displacement transfer per Definition 2.1, and the difference among types of force reduces solely to the difference in the axis along which the transfer occurs.

---

## S3 Origin of Mass: Coupling with Standing Waves

### S3.1 Geometric Structure of Coupling

The standing wave mode ($\kappa = 1$, [3] S5) fills the entire subjective space $S^3$ as a fundamental mode along one spatial axis. A direction-constrained wave ($\kappa = 2$) propagates through this standing wave background.

A direction-constrained wave $\mathbf{k} = (k_x, 0, 0, k_t)$ and a standing wave $\mathbf{k}_H = (k_{H}, 0, 0, 0)$ share the spatial $x$-axis. On this shared axis, the phase structures of the waves overlap, giving rise to coupling.

### S3.2 Modification of the Speed Ratio

The propagation speed of a direction-constrained wave without coupling is $v = k_x / k_t$. Coupling with a standing wave affects the phase structure in the $x$-axis direction, modifying the effective spatial component:

$$k_x \;\longrightarrow\; k_x^{\text{eff}} = k_x - \delta k_H$$

where $\delta k_H > 0$ is proportional to the coupling strength with the standing wave. The sign is negative (subtraction) because the background phase structure of the standing wave overlaps in antiphase with the spatial phase of the direction-constrained wave, coupling in the direction that cancels the effective spatial phase variation. The propagation speed after coupling is:

$$v^{\text{eff}} = \frac{k_x^{\text{eff}}}{k_t} = \frac{k_x - \delta k_H}{k_t} < \frac{k_x}{k_t} = v$$

The speed decreases -- this is the geometric origin of **mass**.

**Proposition 3.1**: Coupling with a standing wave ($\kappa = 1$) reduces the effective spatial component of a direction-constrained wave, restricting its propagation speed to $v < c$. This speed restriction is observed as mass.

### S3.3 Waves That Do Not Couple

Spherical waves ($\kappa_s = 0$) have no spatial components. Although a standing wave has a phase structure along one spatial axis, there is no shared spatial axis with a spherical wave.

**Proposition 3.2**: Spherical waves with $\kappa_s = 0$ have no coupling axis with a standing wave, yielding $\delta k_H = 0$. Their speed remains $v = c$ without modification.

This is consistent with the fact that gauge bosons with $\kappa_s = 0$ ($\gamma$, $g$) are massless in Table A of [4]. Although $W^{\pm}$ and $Z^0$ have $\kappa_s = 0$ and thus lack the direct coupling via spatial axis described in Proposition 3.2, they can acquire mass through an indirect standing wave coupling mechanism via the $t$ or $R$ axis. Proposition 3.2 addresses only coupling through direct sharing of a spatial axis; the details of indirect coupling mechanisms via the $t$ or $R$ axis are outside the scope of this paper.

### S3.4 Mass Hierarchy

The coupling strength $\delta k_H$ depends on the overlap integral between the direction-constrained wave and the standing wave on the spatial axis. Due to the generation structure (selection of spatial axes $k_1, k_2, k_3$) shown in [3] S3.3, the overlap integrals can differ among the three generations.

On $S^3$, the three spatial axes are locally equivalent under $\mathrm{SO}(3)$ symmetry, but they can become non-equivalent in their global coupling with the standing wave mode. This suggests a geometric origin for mass differences between generations, but the derivation of specific mass ratios is outside the scope of this paper.

---

## S4 Electromagnetic Interaction and Scattering

### S4.1 Displacement Transfer by Photons

The photon $\gamma = (0, 0, 0, 0, +, 0)$ is a spherical wave with a non-zero component on the $R$ axis. The $R$ axis corresponds to the radial direction of central projection, and on $\mathbb{Z}^4$ it acts as a change in the dynamic state $\Delta\mathbf{k}$.

When a direction-constrained wave (electron type) $(k_x, 0, 0, k_t)$ absorbs a photon:

$$\Delta\mathbf{k} \;\longrightarrow\; \Delta\mathbf{k} + \delta(\Delta\mathbf{k})_{\gamma}$$

The static label $\mathbf{k}$ does not change (the electron remains an electron). Only the dynamic state changes:

- Change of momentum direction (scattering angle)
- Change of energy (absorption of photon energy)

**Proposition 4.1**: Interaction with a photon preserves the static label of a direction-constrained wave and changes only the dynamic state $\Delta\mathbf{k}$.

### S4.2 Electromagnetic Interaction Between Direction-Constrained Waves

The electromagnetic interaction between two direction-constrained waves $\psi_1$, $\psi_2$ is described as photon exchange:

$$\psi_1 \;\xrightarrow{\text{photon emission}}\; \psi_1' + \gamma \;\xrightarrow{\text{photon absorption}}\; \psi_1' + \psi_2'$$

By conservation law (I2):

$$\Delta\mathbf{k}_1 + \Delta\mathbf{k}_2 = \Delta\mathbf{k}_1' + \Delta\mathbf{k}_2'$$

### S4.3 Sign of Interaction: Attraction and Repulsion

The relative sign of the $t$-axis components of two direction-constrained waves determines the directionality of the interaction.

The $t$-axis component $k_t$ corresponds to the charge structure in Table A of [4]. Photon exchange between two waves with same-sign $k_t$ acts repulsively, while opposite-sign $k_t$ acts attractively.

| $k_t$ of $\psi_1$ | $k_t$ of $\psi_2$ | Sign product | Interaction |
|-------------------|-------------------|--------|---------|
| $+$ | $+$ | $+1$ | Repulsion |
| $+$ | $-$ | $-1$ | Attraction |
| $-$ | $-$ | $+1$ | Repulsion |

**Proposition 4.2**: The sign of the electromagnetic interaction between two direction-constrained waves is determined geometrically by the relative sign of their $t$-axis components, $\mathrm{sign}(k_{t,1} \cdot k_{t,2})$. The specific correspondence "same sign = repulsion, opposite sign = attraction" is an interpretation layer dependent on the correspondence between $k_t$ and charge sign in Table A of [4]; the geometric content of the proposition is that the relative sign uniquely determines the directionality of the interaction.

### S4.4 Pair Annihilation

A direction-constrained wave and its antiparticle have all components with opposite signs ($C$ transformation in [3] S6.5):

$$\psi: (k_x, 0, 0, k_t) \quad \longleftrightarrow \quad \bar{\psi}: (-k_x, 0, 0, -k_t)$$

When both interact in the same lattice region, the sum of the static labels vanishes:

$$(k_x, 0, 0, k_t) + (-k_x, 0, 0, -k_t) = (0, 0, 0, 0)$$

A pair of $\kappa = 2$ waves transitions to a $\kappa = 0$ state. By conservation law (I2), the energy of the vanished static labels is emitted as spherical waves (photons).

**Proposition 4.3**: Pair annihilation of a direction-constrained wave and an anti-direction-constrained wave is the process in which a pair of $\kappa = 2$ waves is converted into $\kappa = 0$ spherical waves, made possible by the vanishing of the sum of static labels.

---

## S5 Gravitational Interaction and Wave Packet Deformation

### S5.1 Displacement Transfer by Gravitons

The graviton $G = (0, 0, 0, \pm t, \pm R, 0)$ has non-zero components on both the $t$ axis and the $R$ axis. By condition (I1), displacement transfer acts in both the $t$ and $R$ directions:

$$\delta(\Delta\mathbf{k})_G = (\delta(\Delta k)_t,\; \delta(\Delta k)_R)$$

This means simultaneous change of energy ($t$ component) and spatial configuration ($R$ component).

### S5.2 Universal Coupling

We recall the description of acceleration introduced in [3] S8. The second-order difference of the child-space central projection determines the acceleration:

$$\Delta^2\mathbf{k}(\mathbf{k}_t) = D^2\pi(\mathbf{k}_t) \cdot \Delta\mathbf{k}_t$$

Graviton exchange is the discrete realization of this $\Delta^2\mathbf{k}$. Since gravitons act on the dynamic state $\Delta\mathbf{k}$:

- They do not depend on the static label $\mathbf{k}$ (regardless of the type of wave)
- They act on all waves with $\Delta\mathbf{k} \neq 0$ (all waves carrying energy)

**Proposition 5.1**: Interaction with gravitons acts only on the dynamic state $\Delta\mathbf{k}$ and does not act on the static label $\mathbf{k}$. Therefore all shape-invariant waves with $\Delta\mathbf{k} \neq 0$ couple to gravitons (universal coupling).

### S5.3 Wave Packet Deformation

Energy transfer through interaction with gravitons deforms the shape of a wave packet. Under the conserved quantity $\int |\psi|^2 \, dV = E_0$ ([1] S5, [2] S5, condition (C1)), the amplitude $A$ and characteristic width $w$ of the wave packet satisfy:

$$A^2 \cdot w^d = E_0 = \text{const} \quad (d = 3: \text{spatial dimension of } S^3)$$

For a wave packet that receives energy:

$$\Delta E > 0 \;\Rightarrow\; A \uparrow,\; w \downarrow \quad (\text{amplitude increase, width decrease})$$

$$\Delta E < 0 \;\Rightarrow\; A \downarrow,\; w \uparrow \quad (\text{amplitude decrease, width increase})$$

**Proposition 5.2**: Under the conserved quantity $\int |\psi|^2 \, dV$, energy transfer deforms the amplitude and width of a wave packet inversely. An increase in energy is accompanied by localization (decrease in width) of the wave packet.

### S5.4 Wave Packet Deformation of Spherical Waves

Spherical waves ($\kappa = 0$) also couple to gravitons (Proposition 5.1: they couple if $\Delta\mathbf{k} \neq 0$). When a spherical wave receives energy:

- The amplitude of the spherical wave increases and its spatial extent decreases
- In the extreme limit, the spherical wave becomes point-like and localized

This corresponds to a geometric description of the gravitational lensing effect. When a spherical wave (photon) passes through the gravitational field of a massive body, the wave packet is deformed by interaction with gravitons, causing the propagation path to bend.

---

## S6 Strong Interaction and Confinement

### S6.1 $Q$-Axis Transition

The gluon $g = (0, 0, 0, 0, 0, Q_i \to Q_j)$ is defined as a transition on the $Q$ axis. By condition (I1), displacement transfer acts only in the $Q$-axis direction:

$$\delta\mathbf{k}_g: \quad Q \;\longrightarrow\; Q'$$

The spatial components, time component, and $R$ component of the direction-constrained wave do not change. Only the discrete label $Q$ changes.

### S6.2 Three-Bit Structure of Color Charge

In Table A of [4], $Q$ is encoded as three bits $(c_1, c_2, c_3)$:

- $c_2 c_3$: color charge ($00$ = colorless, $01, 10, 11$ = three colors)
- $c_1$: weak isospin

The $Q$ transition of a gluon is a change in the $c_2 c_3$ bits, corresponding to color charge exchange. Transitions among the three non-zero values ($01, 10, 11$) of $c_2 c_3$ give $3 \times 3 - 1 = 8$ independent gluons.

### S6.3 Confinement Condition

In the framework of child-space central projection introduced in [3] S8, forming an independent subjective space requires that the collection of waves be **color-neutral**.

**Definition 6.1** (Color neutrality condition): A collection of direction-constrained waves $\{\psi_i\}$ is color-neutral if the composition of the $c_2 c_3$ bits reduces to $(0, 0)$. Here the composition of bits is defined as bitwise exclusive OR (XOR, addition over $\mathbb{F}_2$):

$$(a_1, a_2) \oplus (b_1, b_2) = (a_1 \oplus b_1,\; a_2 \oplus b_2) \quad (\oplus: \text{mod}\; 2 \text{ addition})$$

Realizations of color neutrality:

| Configuration | Color combination | Name |
|------|----------|------|
| Three waves | $(01) \oplus (10) \oplus (11) = (00)$ | Baryon type |
| Two waves | $(c_2 c_3) \oplus (c_2 c_3) = (00)$ | Meson type |
| One wave | $(c_2 c_3) = (00)$ | Lepton type |

**Proposition 6.1**: A collection of direction-constrained waves that does not satisfy the color neutrality condition cannot form an independent subjective space and is confined within a color-neutral collection.

### S6.4 Asymptotic Freedom

Inside a color-neutral collection (inside the shared subjective space), the curvature of the child-space central projection is small, and direction-constrained waves move nearly freely. As the boundary of the collection is approached, the projection curvature increases, making it impossible to cross the boundary.

**Proposition 6.2**: The interior of a color-neutral collection has small projection curvature (asymptotic freedom), while at the boundary the projection curvature diverges (confinement).

---

## S7 Weak Interaction and Flavor Transition

### S7.1 $t$-Axis Displacement

$W^{\pm} = (0, 0, 0, \pm t, 0, 0)$ is a spherical wave with a non-zero component on the $t$ axis. By condition (I1):

$$\delta\mathbf{k}_W: \quad k_t \;\longrightarrow\; k_t' = k_t + \delta k_t$$

The change in the $t$-axis component corresponds to the up-type / down-type conversion in Table A of [4]. Simultaneously, the $c_1$ bit of the $Q$ axis also changes.

### S7.2 Concrete Example of Flavor Transition

Transformation of a direction-constrained wave by $W^+$ absorption:

$$\underbrace{(k_x, 0, 0, -k_t, 0, Q_d)}_{\text{down-type}} + \underbrace{(0, 0, 0, +, 0, 0)}_{W^+} \;\longrightarrow\; \underbrace{(k_x, 0, 0, +k_t, 0, Q_u)}_{\text{up-type}}$$

The sign of the $t$ axis is inverted ($-k_t \to +k_t$), and the isospin bit $c_1$ of $Q$ changes.

### S7.3 $k_t \neq 0$ and CP Violation

As shown in Proposition 6.4 of [3], CP symmetry can be violated only in interactions with $k_t \neq 0$. $W^{\pm}$ has $k_t = \pm 1 \neq 0$, and is the only one among the four types of spherical waves to have $k_t \neq 0$.

**Proposition 7.1**: Only $W^{\pm}$, which has a non-zero component on the $t$ axis, mediates CP symmetry violation. This is a direct consequence of $k_t \neq 0$ introducing asymmetry in the time direction.

### S7.4 Chirality Selection Rule

For the chirality $\chi = \mathrm{sign}(k_{s_0} \cdot k_t)$ ($k_{s_0}$: non-zero spatial component) defined in Definition 6.1 of [3], the $t$-axis structure of $W^{\pm}$ imposes a coupling condition with a specific chirality.

$W^{\pm} = (0, 0, 0, \pm t, 0, 0)$ has a phase structure only on the $t$ axis. Coupling with a direction-constrained wave $(k_{s_0}, 0, 0, k_t)$ is established only when the $t$-axis phase of $W^{\pm}$ and the $t$-axis phase of the direction-constrained wave have opposite signs ($\mathrm{sign}(k_t) \neq \mathrm{sign}(k_{t,W})$). By Definition 6.1 of [3], $\chi = \mathrm{sign}(k_{s_0} \cdot k_t)$, so this coupling condition corresponds to $\chi = -1$ (left-handed).

**Proposition 7.2**: Coupling with $W^{\pm}$ is restricted to direction-constrained waves with chirality $\chi = -1$ (left-handed). Right-handed ($\chi = +1$) waves do not couple with $W^{\pm}$.

---

## S8 Wave Packet Deformation and Retroactive Construction of Causality

### S8.1 Conservation-Law Description of Wave Packet Deformation

We generalize the mechanism of wave packet deformation introduced in S5.3. In any interaction where a direction-constrained wave receives energy, under the conserved quantity $\int |\psi|^2 \, dV = E_0$:

$$A_{\text{after}}^2 \cdot w_{\text{after}}^d = A_{\text{before}}^2 \cdot w_{\text{before}}^d = E_0$$

Concentration of energy density (increase in $A$) necessarily accompanies spatial localization (decrease in $w$).

### S8.2 Dynamical Description of "Observation"

In conventional quantum mechanics, "wave packet collapse upon observation" is introduced axiomatically. In the present framework:

1. A direction-constrained wave ($\kappa = 2$) interacts with a spherical wave ($\kappa = 0$)
2. Energy transfer increases the amplitude of the wave packet
3. By conservation law (I3), the width decreases (localization)
4. A sufficiently localized wave packet is described as "detected as a particle"

"Observation" is not a special operation but merely an instance of wave packet deformation due to interaction. No external "observer" is required, and the deformation of the wave packet is a deterministic process governed by conservation laws.

**Proposition 8.1**: Wave packet collapse is a conservation-law-governed wave packet deformation due to interaction, and there is no need to postulate a special process called "observation."

### S8.3 Retroactive Construction of Causality

As shown in Proposition 4.2 of [3], spherical waves with $k_t = 0$ are time-reversal symmetric and carry no intrinsic causal ordering. The "deterministic" in Proposition 8.1 and the "retroactive" below do not contradict each other. They belong to different layers:

- **Deterministic**: Conservation laws (I2)(I3) hold strictly, and the outcome of wave packet deformation is uniquely determined
- **Retroactive**: Because $k_t = 0$ entails no causal ordering (no temporal order of "cause then effect"), which event is the "cause" and which is the "effect" remains indeterminate

That is, the structure of "strictly obeying conservation laws, yet possessing no causal ordering" arises naturally from the time-reversal symmetry of $k_t = 0$. Combining this property with wave packet deformation:

**(R1) Retroactive determination of path**: A spherical wave ($k_t = 0$) propagates isotropically in all directions (Proposition 4.1 of [3]). Viewed from the outcome of an interaction, the spherical wave **appears to have** propagated along the shortest geodesic. Since no causal ordering exists, the outcome retroactively determines the path.

**(R2) Retroactive determination of localization**: Wave packet deformation is a deterministic process governed by conservation laws (Proposition 8.1). However, because $k_t = 0$ entails no causal ordering, "the interaction occurred, therefore localization happened" and "localization occurred, therefore the interaction appears to have happened" are indistinguishable.

**(R3) Principle of least action**: As the synthesis of (R1) and (R2), the variational principle (principle of least action) of natural law is derived not as an axiom but as a consequence of the causal symmetry of interactions mediated by spherical waves with $k_t = 0$.

**Proposition 8.2**: In interactions mediated by spherical waves with $k_t = 0$, (a) the shortest path of propagation, (b) the localization of the wave packet, and (c) the variational principle are all retroactively constructed from the absence of causal ordering.

### S8.4 The Case of $k_t \neq 0$

Interactions mediated by $W^{\pm}$ ($k_t \neq 0$) are not causally symmetric. (R1)--(R3) do not apply:

- The path is causally determined (emission first, absorption later)
- Time-reversal symmetry is broken
- CP non-conservation becomes possible (Proposition 7.1)

| $k_t$ of spherical wave | Causal ordering | Path selection | CP | Corresponding force |
|---------------|---------|---------|-----|----------|
| $k_t = 0$ | None | Retroactive (geodesic) | Conserved | Electromagnetism, strong force |
| $k_t \neq 0$ | Present | Causal | Can be violated | Weak force |

---

## S9 Summary

We summarize the main results of this paper:

1. **Unified interaction mechanism**: All forces are described as displacement transfer along the non-zero axes of spherical waves via (I1)--(I3), and the type of force reduces solely to the difference in the axis (S2).

2. **Origin of mass**: Coupling with a standing wave ($\kappa = 1$) decreases the speed ratio $|k_x / k_t|$ of a direction-constrained wave, producing mass as $v < c$. Spherical waves with $\kappa_s = 0$ have no coupling axis and remain massless (S3).

3. **Electromagnetic scattering and pair annihilation**: The photon changes only the dynamic state $\Delta\mathbf{k}$. The relative sign of the $t$-axis components determines attraction / repulsion. Pair annihilation of a direction-constrained wave and its anti-wave is the conversion $\kappa = 2 \to \kappa = 0$ (S4).

4. **Gravity and wave packet deformation**: The graviton acts on both the $(t, R)$ axes and couples to all waves with $\Delta\mathbf{k} \neq 0$. Under the conserved quantity, energy transfer deforms the amplitude and width of a wave packet inversely (S5).

5. **Strong force and confinement**: Mediates $Q$-axis transitions (color charge exchange), and the color neutrality condition provides the condition for forming an independent subjective space (S6).

6. **Weak force and CP violation**: Mediates $t$-axis displacement (flavor transition), and $k_t \neq 0$ is the origin of CP violation and the left-handed chirality selection rule (S7).

7. **Wave packet deformation and causal structure**: Wave packet collapse is a conservation-law-governed deformation, and "observation" is unnecessary. Interactions mediated by spherical waves with $k_t = 0$ are causally symmetric, and the path, localization, and variational principle are retroactively constructed (S8).

The correspondence with Table A of [4] is shown below:

| Interaction | Spherical wave | Non-zero axis | Effect on direction-constrained wave | $k_t$ | Causality |
|---------|--------|---------|-----------------|-------|------|
| Higgs coupling | Standing wave $H^0$ | Spatial 1-axis | Decrease of $v$ (mass) | 0 | -- |
| Electromagnetic | $\gamma$ | $R$ | $\Delta\mathbf{k}$ change (scattering) | 0 | Retroactive |
| Strong | $g$ | $Q$ | $Q$ change (color exchange) | 0 | Retroactive |
| Weak | $W^{\pm}$ | $t$ | $k_t$ change (flavor transition) | $\neq 0$ | Causal |
| Gravity | $G$ | $t, R$ | $\Delta\mathbf{k}$ change (acceleration) | $\neq 0$ | Causal |

All claims of this paper are stated as geometric facts, and the physical correspondences are presented as "connectable interpretations."

---

## References

[1] Noriaki Kihara, "Shape-Invariant Standing Waves on Closed Spheres: Stability by Dimension and the Geometric Necessity of 3+1 Spacetime," 2026.

[2] Noriaki Kihara, "Shape-Invariant Traveling Waves on a 4-Dimensional Integer Lattice: Existence Conditions, Conservation Structure, and Self-Consistent Loop," 2026.

[3] Noriaki Kihara, "Four Modes of Shape-Invariant Waves: Spin, Chirality, and Statistics Determined by Wave Vector Structure," 2026.

[4] Noriaki Kihara, "Set Structure of a 6-Dimensional Hyperrectangle and Its Combinatorial Properties," Zenodo, 2025. DOI: 10.5281/zenodo.19688521.

[5] Noriaki Kihara, "Three-Layer Model of the Universe via Central Projection," Zenodo, 2025. DOI: 10.5281/zenodo.15083729.
