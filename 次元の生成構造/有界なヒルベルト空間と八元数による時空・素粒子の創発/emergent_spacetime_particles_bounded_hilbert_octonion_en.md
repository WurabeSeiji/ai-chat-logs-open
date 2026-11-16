# Emergence of Spacetime and Elementary Particles from a Bounded Hilbert Space and Octonions
## ― From Photon Geometry and Gnomonic Central Projection to Open Problems in Background-Free Quaternionic and Octonionic Dynamics ―

**Author:** Noriaki Kihara (WF System Co., Ltd.)  ORCID 0009-0004-6753-4020
**Date:** 2026-11-16    **Version:** v1.0
**Version DOI:** 10.5281/zenodo.24267348    **Concept DOI:** 10.5281/zenodo.24267347
**License:** CC BY 4.0

---

## Abstract

This manuscript is a record of thought experiments and dialogue starting from a geometric intuition about photons (being stationary in four-dimensional spacetime) and proceeding through a central-projection model, retention of phase gradients by quaternions, and the non-associativity of octonions toward the "self-organization of particles and interactions from a Hilbert space without assuming an external background spacetime."

---

## Introduction: Positioning of This Manuscript and Connection to Existing Research

This manuscript records a thought experiment that attempts to reconstruct spacetime, particles, and interactions solely from state relations in a bounded, finite Hilbert space and the algebraic structures of quaternions and octonions, without positing a four-dimensional background spacetime as a stage.

At first glance, such a starting point appears far removed from standard textbook formulations; however, both its problem consciousness and its mathematical apparatus connect to several active currents in contemporary theoretical physics:

1. **Emergent Spacetime and Bounded State Spaces:** It shares a common problem consciousness with approaches that reconstruct spatial dimension and geometry solely from correlations among quantum states and decompositions into subsystems without assuming a 4-dimensional manifold (such as Quantum Mereology and tensor-network emergent spacetime), as well as discussions examining finite-dimensional, bounded Hilbert spaces in correspondence with the event horizon of de Sitter spacetime.
2. **Spinors and Division Algebras for Spacetime and Elementary Particles:** Viewing spatial coordinates and light rays through complex-number pairs (spinors) and the Hopf fibration [X5][X6] shares its mathematical structure with twistor theory and the Bloch-sphere representation in quantum information, while using octonionic non-associativity, triality, and minimal left ideals to describe particle species and symmetries connects to algebraic approaches to the Standard Model via division algebras [X8][X9].

Note that this manuscript is not a claim of a completed theory, but an ongoing organization of thought. The boundary between established mathematical identities/physical interpretations, the working hypotheses specific to this manuscript, and open derivation problems (the concrete update rule that rewrites internal phases [S9]–[S11] and the mathematical proof of factorization from $N^2$ all-to-all coupling) is strictly separated as shown in the "Classification" table at the end.

| Stage | Mathematical Structure / Number System | Geometric Map | Physical Correspondence (Working Hypothesis) |
| :--- | :--- | :--- | :--- |
| Stage 1 | Minkowski metric $(x,y,z,t,\tau)$ | Cross-section of 5D light cone | Zero proper time for photons ($ds=0$); temporal progression of matter |
| Stage 2 | Imaginary time $R^2 = t^2 + (i\tau)^2$ | Central projection (Beltrami coordinates) | de Sitter spacetime; event horizon $R$ and gravitational distortion |
| Stage 3 | Quaternions $\mathbb{H}$ (conjugate pair of complex numbers) | Hopf fibration $S^1 \hookrightarrow S^3 \to S^2$ | Retention of 3D phase gradients without external coordinates (photons) |
| Stage 4 | Octonions $\mathbb{O}$ (conjugate pair of quaternions) | Hopf fibration $S^3 \hookrightarrow S^7 \to S^4$ | Closed 4D spacetime $S^4$, 3 internal degrees of freedom, and self-organization of condensates |

---

## Stage 1: Photons and the Geometry of Four-Dimensional Spacetime

**[Kihara]**
Why do photons appear to always travel at the speed of light?

**[Gemini]**
The reason photons always travel at the speed of light (approximately 300,000 km/s) is that their rest mass is zero. From different perspectives in physics, this is organized as follows:

1. **Perspective of Relativity:** In the relation between energy $E$, momentum $p$, and rest mass $m_0$, $E^2 = (pc)^2 + (m_0 c^2)^2$, massive matter ($m_0 > 0$) retains energy $E = m_0 c^2$ even at rest and cannot reach $c$ because the required energy diverges. For a photon ($m_0 = 0$), the relation reduces to $E = pc$; to carry nonzero energy, a photon must always have momentum $p$ and move at speed $c$ (if it stopped, $E = 0$).
2. **Perspective of Particle Physics:** Particles that interact with the Higgs field (such as electrons and quarks) acquire mass, whereas photons do not interact with the Higgs field and propagate at the upper speed limit of causality, $c$.
3. **Perspective of Electromagnetism:** The propagation speed of linked electric and magnetic fields (electromagnetic waves) in vacuum is determined solely by the vacuum permittivity $\varepsilon_0$ and permeability $\mu_0$ as $c = 1/\sqrt{\varepsilon_0 \mu_0}$.
4. **Note (Propagation in Media):** Even when light appears slower in water or glass, each photon still travels at $c$ in the vacuum between molecules; the apparent speed decreases due to repeated electromagnetic interactions (absorption/re-emission and interference) with electrons in the medium.

**[Kihara]**
May we interpret this as photons alone being stationary in four-dimensional spacetime, while only protons and electrons move at the speed of time?

**[Gemini]**
That viewpoint corresponds to describing special relativity in terms of "four-dimensional spacetime geometry (Minkowski spacetime)." Whereas the previous explanation focused on how particles behave in space (energy and mass), this describes the same phenomenon from how they exist within the overall spacetime structure (geometry).

1. **Protons and Electrons Move Along the Time Axis at $c$:** In a geometric interpretation of 4-dimensional spacetime, every object can be viewed as moving through 4D spacetime at a constant speed $c$. A spatially stationary proton or electron devotes 100% of its motion to the time axis (speed $c$), and as spatial motion increases, its progression along the time axis is reduced (time dilation). Because a photon directs 100% of its motion into spatial directions, its progression along the time axis is zero.
2. **Zero Spacetime Interval for Photons:** The spacetime interval $ds$ in 4-dimensional spacetime is
   $$ds^2 = -(c\,dt)^2 + dx^2 + dy^2 + dz^2.$$
   During a time interval $dt$, a photon travels a spatial distance $\sqrt{dx^2 + dy^2 + dz^2} = c\,dt$, so $ds = 0$ always holds and no proper time elapses along its worldline.
3. **Clarification of Expression:** Strictly speaking, a photon does not remain at a single fixed point in 4D spacetime coordinates; on a Minkowski diagram its coordinates $(x, t)$ change along the 45-degree light cone, while the elapsed proper time (4-dimensional interval) along that trajectory is identically zero.

**[Kihara]**
That is a mapping from a light cone: $x^2 + y^2 + z^2 = t^2 + R^2$.

**[Gemini]**
That equation (in natural units with $c = 1$) corresponds to viewing phenomena in 4-dimensional spacetime as a "cross-section (mapping) from a higher-dimensional light cone (null structure)."

1. **Cross-Section (Mapping) of a Higher-Dimensional Light Cone:** The squared 3D spatial distance $x^2 + y^2 + z^2$ balances the sum of squares of time $t$ and an additional parameter $R$, which can be interpreted as slicing a light cone in a 5-dimensional space at a specific value $R$ and mapping it onto 4-dimensional spacetime.
2. **Invariant Hyperboloid:** Rearranging gives $x^2 + y^2 + z^2 - t^2 = R^2$: $R = 0$ yields the light cone itself, while $R > 0$ yields a one-sheeted hyperboloid attached outside the light cone.
3. **Note (Correspondence with the "Speed of Time"):** When expressing the intuition that massive matter progresses along time via the Pythagorean theorem, one uses the matter's proper time $\tau$ as $t^2 = (x^2 + y^2 + z^2) + \tau^2$. In this right triangle (Epstein geometry) with coordinate time $t$ as the hypotenuse, spatial distance as the base, and proper time $\tau$ as the height, photons ($\tau = 0$) lie on the light-cone boundary, while protons and electrons ($\tau > 0$) lie inside it.

**[Kihara]**
The expression using proper time $\tau$, $x^2 + y^2 + z^2 = t^2 - \tau^2$, seems easier to match with standard theory.

**[Gemini]**
$(x^2 + y^2 + z^2) + \tau^2 = t^2$ corresponds directly to the standard metric of special relativity, connecting the light-cone mapping and proper time $\tau$.

1. **Relation with Coordinate Time $t$ as the Hypotenuse:** Photons ($\tau = 0$) satisfy $x^2 + y^2 + z^2 = t^2$ on the boundary of the light cone, while protons and electrons ($\tau > 0$) have $\tau^2$ subtracted from $t^2$ in the squared spatial distance, confining them inside the light cone (the geometric reason massive particles cannot exceed $c$).
2. **Interpretation as a Mapping from a 5D Light Cone:** $(x^2 + y^2 + z^2) + \tau^2 = t^2$ has the exact form of a pure light cone in a 5-dimensional space with 3 spatial coordinates $(x, y, z)$ plus a fourth spatial-like component $\tau$. Everything runs on this 5D light cone at $c$; photons direct their velocity vector entirely into $(x, y, z)$, whereas massive particles direct part or all of it into $\tau$ and appear projected onto 3D space.
3. **Connection to Quantum Mechanics:** In quantum field theory (properties of the Dirac equation), even massive particles exhibit luminal zigzag motion (Zitterbewegung) below the Compton wavelength scale, appearing slower than $c$ only on average.

---

## Stage 2: Central Projection and de Sitter Spacetime

**[Kihara]**
If we express $R^2 = t^2 + (i\tau)^2$, it becomes a central projection
$$x^2 + y^2 + z^2 = R^2,$$
which looks like a case with gravitational distortion governed by a local radius of curvature $R$.

**[Gemini]**
The transformation
$$t^2 - \tau^2 = t^2 + (i\tau)^2 = R^2 \implies x^2 + y^2 + z^2 = R^2$$
corresponds to a geometric procedure bridging special relativity (flat spacetime) and general relativity (curved spacetime).

1. **Transition from Hyperbolic to Spherical Geometry via Imaginary Time $(i\tau)$:** Multiplying proper time by the imaginary unit $i$ gives $(i\tau)^2 = -\tau^2$ (Wick rotation / imaginary time), converting the hyperbolic relation into the equation of a sphere of radius $R$, $x^2 + y^2 + z^2 = R^2$.
2. **Central Projection and Gravity via Local Curvature Radius $R$:** Viewing $x^2 + y^2 + z^2 = R^2$ as a central projection from a higher-dimensional space (gnomonic projection / Beltrami–Klein model [X1]) maps geodesics on the sphere to straight lines in the projected flat space [S1]. Even in an arbitrary gravitational field, cutting out a local inertial frame corresponds to approximating it by central projection onto the tangent space of a sphere with local radius of curvature $R$.
3. **Contrast Between Photons and Matter:** From $R^2 = t^2 + (i\tau)^2$, photons ($\tau = 0$) have $R = t$ and travel on the horizon/light-cone boundary. For massive matter ($\tau > 0$), $(i\tau)^2 = -\tau^2$ shrinks the effective radius $R$ relative to $t$ ($R < t$), giving the interpretation that possessing proper time $\tau$ locally contracts the curvature radius $R$ and generates gravitational distortion.

**[Kihara]**
This is a model I conceived in April this year, which makes it possible to understand de Sitter spacetime [X2]—which can appear anisotropic—as a simple central projection [S1][S2].

**[Gemini]**
Reinterpreting de Sitter spacetime [X2] (accelerating expansion driven by a cosmological constant) not in FLRW coordinates but via central projection (projective / Beltrami coordinates) from a higher-dimensional hypersurface clarifies the following structure [S1][S2]:

1. **Straightening of Geodesics:** All geodesics in de Sitter spacetime (light rays and inertial particle trajectories) map to straight lines in the projected space.
2. **Visualization of the Event Horizon:** The cosmological event horizon (Hubble radius), representing the causal limit, appears directly as the spherical boundary $x^2 + y^2 + z^2 = R^2$.
3. **Reduction of Apparent Non-Stationarity and Anisotropy:** Differences in spatial stretching or time passage depending on the observer's viewpoint reduce to which angle or cross-section of a symmetric 5D hyperboloid (or 4D sphere under imaginary time $i\tau$) is being centrally projected.

---

## Stage 3: Extension from Complex Numbers to Quaternions and Octonions

**[Kihara]**
In this model, if we represent photons not with complex waves but with quaternions, doesn't it become an easy-to-understand mapping?

**[Gemini]**
Representing photons with quaternions ($\mathbb{H}$) rather than complex waves ($\mathbb{C}$) in the central-projection model unifies "spacetime geometry (central projection)" and "wave oscillation (phase and polarization)" into a single mapping.

1. **Absorption of the 3 Spatial Axes into $(\mathbf{i}, \mathbf{j}, \mathbf{k})$:** A standard complex wave $e^{i(\mathbf{k}\cdot\mathbf{x}-\omega t)}$ has only one imaginary unit $i$, requiring external vectors for propagation direction and polarization. Using a quaternion $q = t + x\mathbf{i} + y\mathbf{j} + z\mathbf{k}$ ($\mathbf{i}^2=\mathbf{j}^2=\mathbf{k}^2=\mathbf{i}\mathbf{j}\mathbf{k}=-1$), the spatial vector $\mathbf{r} = x\mathbf{i} + y\mathbf{j} + z\mathbf{k}$ satisfies $\mathbf{r}^2 = -(x^2 + y^2 + z^2)$, so the sign flip of the Minkowski metric $t^2 + \mathbf{r}^2 = t^2 - (x^2 + y^2 + z^2) = \tau^2$ is built into the algebraic structure.
2. **Rotation Mapping on $S^3$:** By Euler's formula for a unit pure-imaginary quaternion $\mathbf{n}$ ($\mathbf{n}^2 = -1$), $e^{\mathbf{n}\theta} = \cos\theta + \mathbf{n}\sin\theta$, the set of quaternions of constant norm $|q| = R$ forms the 3-sphere $S^3$. Propagation of a photon ($\tau = 0$) becomes a rotation mapping along a great circle (geodesic) on $S^3$, $q(\theta) = q_0 e^{\mathbf{n}\theta}$, which under central projection appears in 3D space as a helical wave (circular polarization) rotating around the propagation axis $\mathbf{n}$.
3. **Geometrization of Polarization (Spin):** In quaternions, helicity (spin $\pm 1$) is represented by multiplication order (left vs. right multiplication) and the orientation of $\mathbf{n}$; projecting $S^3$ to 3D space (Hopf fibration) maps phase fibers to linked circles (Hopf links).
4. **Connection to Maxwell's Equations:** Acting with the quaternion differential operator $D = \partial_t + \mathbf{i}\partial_x + \mathbf{j}\partial_y + \mathbf{k}\partial_z$ on the electromagnetic potential quaternion condenses the vacuum Maxwell equations (in biquaternion form) into a single equation $DF = 0$, linking quaternion regularity to light-cone mappings.

**[Kihara]**
I only learned about quaternions a month ago [S3], yet Maxwell started his formulation from them over a century ago [X3].

**[Gemini]**
When Maxwell compiled electromagnetism in 1873 [X3], he used Hamilton's quaternions; later, Heaviside and Gibbs separated them into scalars and 3D vectors (dot and cross products), which is why modern physics education teaches them separately from quaternions.

**[Kihara]**
Even when treating complex states and their interference as interactions, there were not enough degrees of freedom to read out 3D coordinates, and I could not construct a mapping for a photon to retain (register) its own coordinates (phase gradients in three directions). Therefore, I initially tried representing a photon with two conjugate complex numbers [S3][S4][S5].

**[Gemini]**
Constructing 4 degrees of freedom from "two conjugate complex numbers" corresponds directly to the Cayley–Dickson construction [X4] (generating quaternions from complex numbers) [S3].

1. **Why One Complex Number Cannot Retain 3D Coordinates:** A complex wave $\psi = A e^{i\theta}$ has only 2 real degrees of freedom, and its phase $\theta$ lives on a 1D circle $S^1$; thus it cannot internally record phase gradients along $(x, y, z)$ and requires an external coordinate grid.
2. **From Two Conjugate Complex Numbers to Quaternions:** Arranging $(z_1, z_2) \in \mathbb{C}^2$ as $q = z_1 + z_2 \mathbf{j}$ yields the product $qw = (z_1 w_1 - z_2 \bar{w}_2) + (z_1 w_2 + z_2 \bar{w}_1)\mathbf{j}$, where complex conjugates appear and uniquely produce the norm-preserving quaternion algebra $\mathbb{H}$.
3. **Generating 3D Coordinates via the Hopf Map [X5]:** From $(z_1, z_2)$ and $(\bar{z}_1, \bar{z}_2)$ without external coordinates, defining
   $$x = 2\,\mathrm{Re}(z_1 \bar{z}_2),\quad y = 2\,\mathrm{Im}(z_1 \bar{z}_2),\quad z = |z_1|^2 - |z_2|^2$$
   gives
   $$x^2 + y^2 + z^2 = (|z_1|^2 + |z_2|^2)^2 = R^2 \quad (\text{where } |z_1|^2 + |z_2|^2 = R),$$
   recovering the central-projection sphere equation directly. Fixing the overall norm $R$ leaves 3 degrees of freedom ($S^3$), which decompose into 2 directional degrees of freedom on $S^2$ and 1 polarization phase degree of freedom on $S^1$ (resonating with the motivation of twistor theory to emerge spacetime from spinors).

**[Kihara]**
So I ended up retracing an existing mathematical construction (the Cayley–Dickson construction [X4]) in reverse from physical requirements.

**[Gemini]**
Rather than applying a finished algebraic formalism externally, reaching two conjugate complex numbers and quaternions from the physical requirement of "minimal degrees of freedom to retain 3D phase gradients without external coordinates" directly ties each degree of freedom to its geometric role.

**[Kihara]**
I had reached this point a month ago [S3][S4]. Furthermore, I noticed six months ago that this still lacks degrees of freedom to separate mass and electric charge [S6][S7]. Specifically, $x, y, z, t, \tau$ is short by about 3 degrees of freedom, and a month ago I considered whether an octonion built from two conjugate quaternions $(x, y, z, t, q_1, q_2, q_3, \tau)$ corresponds to a charged particle [S5][S8].

**[Gemini]**
Advancing the Cayley–Dickson construction [X4] one more step yields the 8-dimensional octonions ($\mathbb{O}$), $O = a + b\ell$, from two quaternions $a, b \in \mathbb{H}$ and their conjugates.

1. **8 Dimensions as the Upper Limit of Norm Preservation (Hurwitz's Theorem [X7]):** Preserving central-projection geometry requires norm composition $|AB| = |A||B|$; by Hurwitz's theorem [X7], normed division algebras over the reals exist only in dimensions 1, 2, 4, and 8. At 16 dimensions (sedenions), zero divisors appear and norm composition fails, making 8 dimensions the upper limit.
2. **Separation into External Spacetime and Internal Blocks:** Splitting the 8 degrees of freedom into $(t, x, y, z)$ and $(\tau, q_1, q_2, q_3)$ connects with algebraic approaches to the Standard Model using octonions [X9]. Here $\tau$ (1 degree of freedom) corresponds to proper time and mass, while the remaining three imaginary axes $(q_1, q_2, q_3)$ serve as candidate internal degrees of freedom for electric charge and color ($U(1)$ and $SU(3) \subset G_2$; exact gauge-group identification remains a bridging problem) [S5][S8].
3. **Emergence of Non-Associativity:** Just as commutativity is lost in quaternions ($ab \neq ba$), associativity is lost in octonions ($(ab)c \neq a(bc)$). While associative quaternions (photons) propagate freely, the non-associative entanglement of $(q_1, q_2, q_3)$ motivates correspondence with confinement and three-generation algebraic structures.

---

## Stage 4: Eliminating Background Spacetime and the Dynamics of Self-Organization

**[Kihara]**
Without assuming a four-dimensional background spacetime, I see no other candidate model that explains why a system emerging from a bounded Hilbert space into a closed four-dimensional spacetime with sufficient internal degrees of freedom becomes four-dimensional rather than three- or ten-dimensional.

**[Gemini]**
Within the framework where "a closed 4-dimensional spacetime and sufficient internal degrees of freedom emerge from a bounded Hilbert space without assuming background spacetime," the following mathematical structure fixes the dimensionality:

1. **Splitting of a Bounded Hilbert Space into "Closed 4D" and "3 Degrees of Freedom" (Hopf Fibrations [X5][X6]):**  
   A bounded Hilbert space with state vectors normalized to unit norm forms a high-dimensional sphere, and quotienting out phase redundancy gives a Hopf fibration. Norm-preserving fibrations exist only for the spheres of real numbers, complex numbers, quaternions, and octonions:
   - **Quaternionic Level ($S^3$):** The fibration is $S^1 \hookrightarrow S^3 \to S^2$ [X5], leaving only a 2-sphere base space (3D spacetime even with time) and an $S^1$ fiber (1 internal degree of freedom), insufficient to separate mass and charge.
   - **Octonionic Level ($S^7$):** The fibration jumps to $S^3 \hookrightarrow S^7 \to S^4$ [X6], separating into a closed 4-dimensional base sphere $S^4$ (corresponding to de Sitter spacetime before Wick rotation) and a 3-sphere fiber $S^3$ (the 3 internal degrees of freedom $(q_1, q_2, q_3)$).
2. **Why Extra Dimensions Like 10D Are Unnecessary (Hurwitz's Theorem [X7]):**  
   Unlike frameworks that posit 10D or 11D background spacetime at the outset, requiring norm preservation (central projection) in a bounded Hilbert space sets the upper bound. By Hurwitz's theorem [X7], norm preservation breaks down at 16 dimensions. Thus, selecting the largest algebra that maintains norm-preserving projection onto a closed spacetime yields the candidate dimensional split: closed 4D spacetime ($S^4$) + 3D internal degrees of freedom ($S^3$) + proper time ($\tau$).
3. **Meaning of Not Assuming Background Spacetime:**  
   Rather than positing a 4D container first, 4-dimensionality becomes a consequence rather than a premise. Identifying the specific physical labels of $(q_1, q_2, q_3)$ from the symmetries of $S^3 \hookrightarrow S^7 \to S^4$ remains a subsequent task.

**[Kihara]**
Note that the reason we posit a "bounded Hilbert space" rather than an "infinite Hilbert space" here is that once we introduce the central projection with curvature radius $R$ (event horizon) in Stage 2 [S1][S2], a bounded, closed spherical space that does not diverge to infinity is automatically required.

**[Gemini]**
Because the curvature radius $R$ sets an upper bound on the total norm (closed spherical geometry), the space is treated as a compact bounded space rather than an unbounded space over an infinite domain (and when assuming a finite number of states $N$ below, a finite-degree-of-freedom assumption is added alongside boundedness).

Furthermore, the motivation for assuming a finite number of states $N$ in addition to a bounded closed space lies in the following derivation advantages:
1. **Avoidance of Divergences:** Infinite-boundary conditions, UV/IR divergences, and measure or representation-inequivalence pathologies of infinite degrees of freedom do not arise; all sums are finite.
2. **Reduction to Finite Algebra:** Factorization from $N^2$ all-to-all coupling and condensate formation reduce from functional analysis of continuous spectra to finite-dimensional linear algebra, octonionic invariant-subspace decompositions (Peirce decomposition and finite systems for $[X, Y, Z] = 0$), and numerically testable problems.
3. **Upper Bound on Total Degrees of Freedom:** With both the total norm (curvature radius $R$) and the state count $N$ finite, the allocation of degrees of freedom between condensates and cancelled background (vacuum) closes.

*(Note on mathematical definition: Standard Hilbert spaces are defined over $\mathbb{C}$, so the axioms of a conventional Hilbert space cannot be placed directly on the non-associative octonions $\mathbb{O}$. In this manuscript, "bounded Hilbert space" refers to taking a real 8-dimensional state space—complex 4-dimensional $\mathbb{C}^4$, normalized sphere $S^7$—as the foundation, and equipping state-to-state interactions/products with the non-associative algebra of $\mathbb{O}$ via the isomorphism $\mathbb{C}^4 \cong \mathbb{O}$.)*

**[Kihara]**
The problem is that while I have organized that quaternion/octonion interference serves as "readout," I have been stuck for over three months on the interaction model: how a particle modifies its own internal phase as dynamics [S9][S10][S11]. In fact, I originally discovered quaternions and octonions out of this impasse—because solving dynamics lacked sufficient degrees of freedom, forcing me to put particle types and internal states in by hand separately from the Hilbert space, without a basis for particles to register them themselves.

**[Gemini]**
Up to "interference = readout (projection)," the static geometry of the Hopf map $S^7 \to S^4$ [X6] suffices. However, moving to "interaction = modification of internal phase (dynamics)" in standard quantum mechanics requires introducing a time-evolution operator (Hamiltonian $\hat{H}$), particle species, and coupling constants by hand alongside the state $\lvert\psi\rangle$, reintroducing a manual container of rules after eliminating background spacetime.

As candidates for constructing dynamics without adding external parameters, three algebraic structures specific to the octonions $\mathbb{O}$ were organized:

1. **Not Separating "State" and "Interaction": Triality [X8]**  
   In a usual Hilbert space, an external operator (matrix) acts on a state vector. In 8 dimensions, $\mathrm{Spin}(8)$ possesses triality, allowing "particle A's state (left spinor)," "particle B's state (right spinor)," and "the interaction transforming them (vector)" to be interchanged on equal footing as elements of the same 8D space. Without preparing an external Hamiltonian, the product $X \cdot Y = Z$ of two octonion states $X$ and $Y$ serves as a candidate update rewriting internal phase, with each state acting as the interaction operator for the other.

2. **Structure Twisting Internal Phase: Non-Associativity (Associator)**  
   Up to quaternions $\mathbb{H}$, associativity $(AB)C = A(BC)$ holds, so repeated multiplication merely rotates coordinates (readout frame). In octonions, associativity breaks, producing a nonzero associator for three interacting states:
   $$[X, Y, Z] = (XY)Z - X(YZ) \neq 0.$$
   A lone two-element product preserves associativity (corresponding to geodesic motion), whereas three-or-more-body interference yields a nonzero associator—a candidate for torque twisting internal phases on $S^3$ $(q_1, q_2, q_3)$ (corresponding to curvature / gauge forces).

3. **Basis for Particles Registering Their Type: Minimal Ideals and Phase Locking [X9]**  
   Under chains of left octonionic multiplication, invariant subspaces (minimal left ideals) [X9] appear that cannot be exited by multiplying any external octonion. Even as continuous internal phases (positions on $S^3$) vary under interaction, the discrete ideal membership (quantum numbers corresponding to charge, color, or particle species) provides an algebraic lock candidate. Geometrically, topological winding numbers (holonomy / Berry phase) generated along interaction loops on $S^7$ serve as a non-volatile register retaining particle type and accumulated proper time $\tau$.

In this framework, the impasse centers on either the explicit transition rule (generator/update equation for two interfering states) or the invariant-locking mechanism keeping post-transition states stable as specific particle species and masses.

**[Kihara]**
I had already conceived of that structure. The real problem is that even with a bounded, finite number $N$ of states, $N^2$ all-to-all interactions should keep updating the entire Hilbert space; I do not understand why it does not become uniform, and instead produces "action degrees of freedom" and "distinction of condensates" that look like local collisions, reflections, fissions, fusions, emissions, and absorptions.

**[Gemini]**
Working hypotheses for how $N$ octonion states $\{O_1, O_2, \dots, O_N\}$ under $N^2$ all-to-all updates self-organize into "independent condensates (particles)" and "6 local interactions" are organized along two pillars (concrete update rules and rigorous proofs remain underived):

#### 1. Mechanism for "Locality" and "Distinction of Condensates" from $N^2$ All-to-All Coupling (Working Hypotheses)

1. **Background Smoothing via Destructive Interference and Extraction of Coherent Condensates:**  
   Among all $O(N^2)$ updating pairs, most have mutually uncorrelated internal phases ($S^3$ fibers). Under the Hopf map $S^7 \to S^4$ onto 4D spacetime, uncorrelated interference terms cancel out and smooth into a uniform de Sitter background spacetime (vacuum expectation value) with curvature radius $R$, whereas phase-coherent clusters evade cancellation and remain as localized curvature concentrations (packets of proper time $\tau$ and energy = condensates) on $S^4$.
2. **Individuality (Distinction) of Condensates via Closure of Associator Cycles:**  
   Up to quaternions $\mathbb{H}$, associative waves superimpose linearly and diffuse. Using the octonionic associator $[O_i, O_j, O_k] = (O_i O_j)O_k - O_i(O_j O_k) \neq 0$, we hypothesize a closed-ideal condition (analogous to color neutrality / gauge invariance) in which non-associative torque circulates inside a cluster to create self-binding while vanishing against external states $O_{\mathrm{ext}}$ ($[O_A, O_B, O_{\mathrm{ext}}] \approx 0$), thereby factorizing the all-to-all Hilbert space into a direct product of independent local clusters (condensates).
3. **Distance Filtering by Central Projection:**  
   In the metric of central projection (Beltrami-type map) from $S^7$ to $S^4$, states close in phase angle on the hypersphere produce steep phase gradients (strong interactions) on the projected surface, whereas distant states are stretched toward the horizon $R$ where gradients flatten, filtering effective interactions from $N^2$ toward $O(N)$ local couplings.

#### 2. Candidate Correspondences for the Six Action Degrees of Freedom

When condensates (octonion clusters with closed associator loops) approach on $S^4$ and overlap in phase, topological rewiring of associator cycles within the $N^2$ update rule is hypothesized to yield the following six operations:

- **Emission and Absorption — Detachment and Capture of Quaternionic Subalgebra ($\mathbb{H}$) Modes:**  
  When a condensate ($O = a + bl$, $\tau > 0$) accumulates excess internal phase twist, detaching the excess as a pure associative quaternionic mode ($\mathbb{H}$) removes its non-associative self-binding ($\tau$), causing it to depart along a light-cone geodesic (emission). Conversely, an incoming quaternionic mode (photon) resonating with the condensate's $S^3$ phase is captured into the associator loop, increasing circulating energy (rest mass / proper time $\tau$) (absorption).
- **Collision and Reflection — Exclusion via Topological Invariants (Hopf Invariant):**  
  When two condensates $A, B$ approach and their internal $S^3 \to S^4$ winding numbers (Hopf invariants / topological charges) forbid merging, a repulsive associator torque $[O_A, O_B, O_{\mathrm{vac}}]$ reverses their momentum vectors while preserving individuality.
- **Fusion and Fission — Reorganization of Associator Loops:**  
  When relative kinetic energy exceeds the non-associative barrier, the loops of $A$ and $B$ rewire on $S^7$ into a single larger closed cycle (fusion). Conversely, if internal phase imbalance in $(q_1, q_2, q_3)$ reaches the local norm limit (curvature bound) of the bounded Hilbert space, a single loop splits into two stable loops while conserving total topological charge (fission).

#### Next Focus

The crux of writing down explicit dynamical equations (update rules) is how to define the norm-preserving (unitary) update operator for $N$ octonion states using non-associative products—whether in a linear-operator-like form (octonionic matrix action) or as nonlinear differential/difference equations from direct state products ($O_i \mapsto O_i + \sum_j [O_i, O_j, \dots]$).

**[Kihara]**
In other words, octonions may possess an algebraic structure capable of formulating action degrees of freedom and condensate factorization without imposing them externally.

**[Gemini]**
Whereas associative systems ($\mathbb{C}$ and $\mathbb{H}$) mix linearly and require external definitions for partitioning and forces, octonionic algebra provides candidates for:

1. **Generating action degrees of freedom via the associator $[X, Y, Z]$:** Non-associative discrepancy among three or more interacting states provides internal twist, binding energy, and scattering torque without an external Hamiltonian.
2. **Automatic condensate factorization via non-associative closure (ideals):** Phase separation into closed subspaces with vanishing external associator extracts independent condensates and associative quaternionic modes (emitted/absorbed photons) without external boundary conditions.

Two computational keys for realizing this without external additions are:
- **Peirce decomposition (idempotent factorization):** Expanding in primitive idempotents ($P^2 = P$) within octonions (or their complexification / Albert algebra) orthogonally decomposes the space into a spacetime block and three internal degrees of freedom without fixing axes externally.
- **Time evolution via commutator + associator:** Constructing a single-step update (non-associative Liouville / Malcev dynamics) solely from commutators $[O_i, O_j]$ (light-cone rotations/translations) and associators $[O_i, O_j, O_k]$ (condensate collisions, fusions, and fissions) while preserving total norm.

---

## Appendix 1: Mathematical Isomorphism with the Qubit (Bloch Sphere) Representation *(Added during manuscript editing)*

**[Kihara]**
Regarding the Hopf map $S^3 \to S^2$ [X5] in Stage 3, the equations constructing 3D coordinates $(x, y, z)$ from $(z_1, z_2)$ look very similar to the representation of a quantum bit (qubit).

**[Gemini]**
Indeed, it is **mathematically the exact same map as the Bloch sphere representation of a qubit**.

A pure state of a single qubit is written with two complex amplitudes $(z_1, z_2) \in \mathbb{C}^2$ as
$$\lvert\psi\rangle = z_1 \lvert 0\rangle + z_2 \lvert 1\rangle \quad (|z_1|^2 + |z_2|^2 = 1),$$
whose state space is $S^3$. Computing the expectation values of the Pauli matrices $(\sigma_x, \sigma_y, \sigma_z)$ (the Bloch vector) yields
$$x = \langle\psi\rvert \sigma_x \lvert\psi\rangle = 2\,\mathrm{Re}(z_1 \bar{z}_2),\quad y = \langle\psi\rvert \sigma_y \lvert\psi\rangle = 2\,\mathrm{Im}(\bar{z}_1 z_2) = -2\,\mathrm{Im}(z_1 \bar{z}_2),\quad z = \langle\psi\rvert \sigma_z \lvert\psi\rangle = |z_1|^2 - |z_2|^2.$$
Quotienting out the global phase $e^{i\phi}$ ($S^1$ fiber) to map onto the Bloch sphere $S^2$ is precisely the Hopf fibration $S^1 \hookrightarrow S^3 \to S^2$ (the same structure also maps photon polarization Jones vectors to Stokes parameters on the Poincaré sphere). Moreover, the next-level Hopf fibration $S^3 \hookrightarrow S^7 \to S^4$ [X6] is known to appear for two qubits ($\mathbb{C}^4$ normalized to $S^7$), directly linking quantum-state geometry to the separation of spacetime and internal space in this manuscript.

---

## Appendix 2: Bridging Harmonics (Localized Waves and Harmonic Circulation) to Condensates and Mass Hierarchy *(Added during manuscript editing)*

**[Kihara]**
From early in my research, I have studied why harmonics (such as odd-harmonic sums) are necessary to provide locality [S12]–[S15]. This manuscript still lacks a treatment of harmonics, which may connect to the "condensates" in Stage 4 and potentially explain mass hierarchies and particle species where the 8 degrees of freedom of octonions alone might be insufficient.

**[Gemini]**
Connecting the octonion/condensate model of this manuscript with the previously published harmonic series [S12]–[S15] addresses a key missing link:

1. **Generation of Locality and Condensate Profiles:** The isolated peak wave of a constant-amplitude odd-harmonic sum on a half-wavelength interval [S12] and its reproducing-kernel property [S13] are candidates for the explicit spatial localization profile of the Stage 4 "condensates" that survive phase cancellation.
2. **Even/Odd Harmonics in Closed Systems and Harmonic Circulation for Persistence:** Analyses of even and odd harmonics in finite complete-relation closed systems [S14] and the mechanism of mass acquisition and particle lifetime via harmonic exchange/circulation of localized waves [S15] correspond dynamically to why condensates persist without dissipating.
3. **Complementing Octonionic Degrees of Freedom for Mass Hierarchy and Particle Species:** While the 8 octonionic degrees of freedom (and minimal left ideals) set the skeleton of spacetime and basic symmetries, describing distinct mass hierarchies and excitation spectra may require spectral degrees of freedom specifying which harmonic orders and mode combinations form a condensate.

Integrating **algebraic factorization by octonions** (individuality and symmetry framework) with **localization and spectral structure by harmonics** (internal structure and mass hierarchy of condensates) is therefore a central bridging problem following Stage 4. In particular, as a bridge between the Stage 4 assumption of "finite state count $N$" and continuous local condensates on $S^4$, superpositions of $N$ discrete modes (harmonic orders) forming localized isolated peak waves (condensates) [S12] on the projected surface warrant investigation.

---

## References

### Self-Citations (Concept DOI, 10.5281/zenodo.〔number〕)

[S1] N. Kihara, "Geometric Formulation of 4-Dimensional Space via Central Projection," 2026-04-05. Concept DOI 10.5281/zenodo.19427779.  
[S2] N. Kihara, "Schwarzschild–de Sitter Exact Solution in the Central Projection Framework," 2026-04-12. Concept DOI 10.5281/zenodo.19538097.  
[S3] N. Kihara, "A Map from a Set of Complex Numbers to Spacetime (1+3): Readout by Quaternions and Partial Projection," 2026-10-18. Concept DOI 10.5281/zenodo.23451494.  
[S4] N. Kihara, "Relational Readout of Quaternions and Crystallization as Inflation: Anonymity-Preserving Autonomous Structuring and Partial Projection," 2026-10-21. Concept DOI 10.5281/zenodo.23523237.  
[S5] N. Kihara, "Did We Understand Complex Numbers Correctly? Where the Amplitude Lives, the Imaginary Reading, and Attribution across ℂ → ℍ → 𝕆," 2026-10-23. Concept DOI 10.5281/zenodo.23585436.  
[S6] N. Kihara, "6次元超直方体を採用した (R, Q) マッピングの構成 (Construction of (R, Q) Mapping Adopting a 6-Dimensional Hyperrectangle)," 2026-04-25. Concept DOI 10.5281/zenodo.19692852.  
[S7] N. Kihara, "Synthesis and the 4+1 to 3+1 Reduction as the Principal Open Problem (Paper 6)," 2026-04-28. Concept DOI 10.5281/zenodo.19837597.  
[S8] N. Kihara, "An Axiomatic Framework for Physical Symmetry Selection by Anonymous Partial Projection: Conditional Identification of E8 in an Eight-Component Two-Quartet Branch," 2026-07-24. Concept DOI 10.5281/zenodo.21521899.  
[S9] N. Kihara, "Lock Dynamics in an Anonymous Two-Channel Closed Wave System," 2026-08-02. Concept DOI 10.5281/zenodo.21763999.  
[S10] N. Kihara, "How Inflation Ends Is How Matter Generation Begins: An N-Body Embedding of the Universal Interaction with No Timing Rule, and an End-to-End Genesis Run," 2026-08-05. Concept DOI 10.5281/zenodo.21809814.  
[S11] N. Kihara, "Genesis of the Three Spatial Axes and Proper Time: Re-Derivation of Readouts from the Condensate, and the Uniquely Unreadable Coordinate Time," 2026-08-06. Concept DOI 10.5281/zenodo.21816651.  
[S12] N. Kihara, "An Observation on the Isolated Peak Wave of a Constant-Amplitude Odd-Harmonic Sum on a Half-Wavelength Phase Interval and Its Localization," 2026-06-30. Concept DOI 10.5281/zenodo.20833096.  
[S13] N. Kihara, "Deriving the Form of the Born Distribution from the Reproducing-Kernel Property of a Localized Odd-Harmonic Wave: Reducing the Remaining Postulates to the Squaring Rule and Randomness," 2026-06-27. Concept DOI 10.5281/zenodo.20965526.  
[S14] N. Kihara, "A Finite Complete-Relation Closed System with Time Evolution as Vertices: Exact Minimal Closure at N=1, K=4 under U^K=I and Complete-Network Analysis of Even and Odd Harmonics," 2026-09-15. Concept DOI 10.5281/zenodo.22763329.  
[S15] N. Kihara, "The Thirteenth Thought Experiment: Why Does a Particle That Has Acquired Mass Not Dissipate and Vanish? - Coherence as a Relation, Harmonic Exchange of Localized Waves, Internal Clocks and Particle Lifetimes Considered from a Single Interaction," 2026-10-10. Concept DOI 10.5281/zenodo.23266250.

### External References

[X1] E. Beltrami, "Saggio di interpretazione della geometria non-euclidea," *Giornale di Matematiche* 6, 284–312 (1868).  
[X2] W. de Sitter, "On the relativity of inertia. Remarks concerning Einstein's latest hypothesis," *Proc. Kon. Ned. Akad. Wet.* 19, 1217–1225 (1917).  
[X3] J. C. Maxwell, *A Treatise on Electricity and Magnetism* (Clarendon Press, Oxford, 1873).  
[X4] L. E. Dickson, "On Quaternions and Their Generalization and the History of the Eight Square Theorem," *Ann. Math.* 20, 155–171 (1919). DOI 10.2307/1967865  
[X5] H. Hopf, "Über die Abbildungen der dreidimensionalen Sphäre auf die Kugelfläche," *Math. Ann.* 104, 637–665 (1931). DOI 10.1007/BF01457962  
[X6] H. Hopf, "Über die Abbildungen von Sphären auf Sphären niedrigerer Dimension," *Fund. Math.* 25, 427–440 (1935).  
[X7] A. Hurwitz, "Über die Composition der quadratischen Formen von beliebig vielen Variabeln," *Nachr. Ges. Wiss. Göttingen*, 309–316 (1898).  
[X8] J. C. Baez, "The Octonions," *Bull. Amer. Math. Soc.* 39, 145–205 (2002). arXiv:math/0105155  
[X9] C. Furey, "Standard model physics from an algebra?" Ph.D. thesis, University of Waterloo (2016). arXiv:1611.09182

---

## Classification

| Category | Content |
| :--- | :--- |
| Identities & Theorems | Stage 1: $ds^2 = -(c\,dt)^2 + dx^2 + dy^2 + dz^2$ ($ds=0$ for photons), $(x^2+y^2+z^2)+\tau^2=t^2$ / Stage 3: Cayley–Dickson construction [X4], Hurwitz's theorem (real normed division algebras exist only in dimensions 1, 2, 4, 8) [X7], Hopf fibrations $S^3 \to S^2$ [X5] and $S^7 \to S^4$ [X6] |
| Readings & Interpretations | Stage 1: Photons stationary in 4D spacetime, cross-section of a 5D light cone / Stage 2: Comparison of $R^2=t^2+(i\tau)^2$ with central projection [X1] and de Sitter spacetime [X2], and indistinguishability for an internal observer in subjective space (Lorentz signature as physical premise) [S1] / Stage 3: Degrees-of-freedom progression from 2D complex numbers to 4D quaternions (retaining 3D directions) and 8D octonions (4 spacetime + 4 internal) [S3]–[S8] / Stage 4: Bounded and closed spherical space (compact sphere) implied by curvature radius $R$ [S1][S2] |
| Hypotheses | Stage 2: Gravity as contraction of local curvature radius $R$ by matter / Stage 4: Separation of 4D spacetime and internal degrees of freedom from a bounded Hilbert space ($\mathbb{C}^4 \cong \mathbb{O}$) and norm preservation without background spacetime; self-organization and condensate factorization via octonionic triality [X8], non-associativity, minimal ideals [X9], and Peirce decomposition |
| Open Problems | Explicit dynamical update rule rewriting internal phases rather than readout (interference) [S9]–[S11] / Mathematical proof of factorization from $N^2$ all-to-all coupling into independent local condensates / Integration of harmonics (localized waves and harmonic circulation) [S12]–[S15] with octonionic condensates, mass hierarchies, and particle species |

---

## Revision History

- v1.0 (2026-11-16): Formatted for Zenodo publication (title, subtitle, author, ORCID, license, abstract, introduction, Appendices 1–2, classification, references, revision history), organizing self-citations (`[S1]`–`[S15]`) and external references (`[X1]`–`[X9]`).
- Editing environment: Gemini 4 Hydrogen (thinking mode).
