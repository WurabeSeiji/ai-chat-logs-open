# The sine-Gordon Equation — Foundations of Topological Soliton Theory

## — Classification of Topological Solutions and Collision Properties in Nonlinear Wave Equations —

**Author:** Noriaki Kihara (WF System Co., Ltd.)

**Date:** April 2026

**DOI:** [10.5281/zenodo.19650966](https://doi.org/10.5281/zenodo.19650966)

---

## 0. Purpose

This paper surveys the mathematical structure of the **sine-Gordon equation**, one of the fundamental nonlinear wave equations. The sine-Gordon equation admits topological solitons as solutions and is widely applied in particle physics, condensed-matter physics, and field theory.

We review the definition of the sine-Gordon equation, the classification of its solutions (kink, antikink, breather), phase shifts upon collision, and physical applications, tracing back to the original literature. Furthermore, we prove that soliton solutions are exactly characterized by a finite number of parameters, and that all structural properties of solitons are preserved in the rectangular-wave limit $\sigma \to 0$. In addition, we explicitly identify the limitations of the sine-Gordon theory (§9.8) and present an alternative model that overcomes them: a formulation based on wave equations for perfectly elastic collisions at arbitrary phases in a perfect fluid with different specific gravities (§11).

---

## 1. Definition of the sine-Gordon Equation

### 1.1 The Equation

The 1+1-dimensional sine-Gordon equation is defined as follows [1, 2].

$$\frac{\partial^2 \phi}{\partial t^2} - \frac{\partial^2 \phi}{\partial x^2} + \sin(\phi) = 0$$

Here $\phi = \phi(x, t)$ is a scalar field (phase field), $x$ is the spatial coordinate, and $t$ is the time coordinate.

Natural units ($c = 1$, $m = 1$) are employed. Restoring the general mass parameter $m$ and the speed of light $c$,

$$\frac{1}{c^2}\frac{\partial^2 \phi}{\partial t^2} - \frac{\partial^2 \phi}{\partial x^2} + \frac{m^2 c^2}{\hbar^2}\sin(\phi) = 0$$

### 1.2 Origin of the Name

The name "sine-Gordon equation" derives from its analogy with the Klein-Gordon equation

$$\frac{\partial^2 \phi}{\partial t^2} - \frac{\partial^2 \phi}{\partial x^2} + m^2 \phi = 0$$

in which the linear term $m^2 \phi$ is replaced by the nonlinear term $\sin(\phi)$. The name is a portmanteau of "sine" and "Gordon," introduced by Rubinstein [3].

### 1.3 Lagrangian Density

The sine-Gordon equation is derived from the following Lagrangian density [4].

$$\mathcal{L} = \frac{1}{2}\left(\frac{\partial \phi}{\partial t}\right)^2 - \frac{1}{2}\left(\frac{\partial \phi}{\partial x}\right)^2 - (1 - \cos\phi)$$

The potential $V(\phi) = 1 - \cos\phi$ attains its minimum value of 0 at $\phi = 2n\pi$ ($n \in \mathbb{Z}$). The existence of these infinitely many degenerate vacua is the origin of topological solitons.

### 1.4 Historical Background

The sine-Gordon equation was independently introduced in the following contexts.

- **Crystal dislocation theory** (1939): Frenkel and Kontorova introduced it as a model describing the motion of dislocations in crystal lattices [1].
- **Differential geometry** (19th century): Bour [5] and Enneper [6] derived an isomorphic equation in the theory of surfaces with constant negative curvature (pseudospheres). Bianchi [7] studied the Bäcklund transformation of this equation.
- **Particle physics** (1960s): Skyrme used a sine-Gordon-type equation as a nonlinear field theory model for baryons [8]. Perring and Skyrme carried out the first numerical computation of kink-antikink collisions [9].

---

## 2. Classification of Solutions

### 2.1 Kink Solution (Winding Number $+1$)

The kink solution is a localized wave in which $\phi$ smoothly transitions from $0$ to $2\pi$ [2, 9].

$$\phi_K(x, t) = 4\arctan\left(\exp\left(\frac{x - vt}{\sqrt{1 - v^2}}\right)\right)$$

Here $v$ ($|v| < 1$) is the velocity of the kink.

**Properties:**

- $\phi \to 0$ as $x \to -\infty$, $\phi \to 2\pi$ as $x \to +\infty$
- The phase increases by $2\pi$ — **winding number (topological charge)** $N = +1$
- The denominator $\sqrt{1 - v^2}$ is the Lorentz factor, ensuring special-relativistic covariance
- Rest energy $E_0 = 8m$ (in natural units $E_0 = 8$)

### 2.2 Antikink Solution (Winding Number $-1$)

The antikink solution is a localized wave in which $\phi$ transitions from $2\pi$ to $0$ [2].

$$\phi_{\bar{K}}(x, t) = 4\arctan\left(\exp\left(-\frac{x - vt}{\sqrt{1 - v^2}}\right)\right)$$

**Properties:**

- $\phi \to 2\pi$ as $x \to -\infty$, $\phi \to 0$ as $x \to +\infty$
- The phase decreases by $2\pi$ — **winding number** $N = -1$
- Same mass and velocity characteristics as the kink

### 2.3 Definition and Conservation of Winding Number

**Definition 2.1.** The **winding number** (topological charge) of a field configuration $\phi(x, t)$ is defined as [4, 10]:

$$N = \frac{1}{2\pi}\int_{-\infty}^{+\infty} \frac{\partial \phi}{\partial x}\,dx = \frac{\phi(+\infty) - \phi(-\infty)}{2\pi}$$

**Proposition 2.1 (Conservation of winding number).** Under the time evolution of the sine-Gordon equation, the winding number $N$ is conserved.

*Proof.* The boundary values $\phi(\pm\infty)$ are fixed at vacuum values ($2n\pi$) and cannot change by an integer multiple of $2\pi$ under continuous time evolution. Therefore $N$ is conserved as a time-independent integer. $\square$

This conservation law guarantees the "particle nature" of solitons — creation and annihilation are possible only in pairs (kink + antikink).

### 2.4 Breather Solution (Bound Oscillatory State)

A breather is an oscillatory solution in which a kink and an antikink are bound together, carrying winding number $N = 0$ [2, 11].

$$\phi_B(x, t) = 4\arctan\left(\frac{\sqrt{1 - \omega^2}\;\cos(\omega t)}{\omega\;\cosh\left(\sqrt{1 - \omega^2}\;x\right)}\right)$$

Here $\omega$ ($0 < \omega < 1$) is the internal oscillation frequency of the breather.

**Properties:**

- A localized oscillatory solution that decays exponentially in space
- A smaller $\omega$ yields a broader spatial extent of the breather
- Rest mass of the breather: $M_B = 16\sqrt{1 - \omega^2}$
- In the limit $\omega \to 0$, the breather mass approaches that of a kink-antikink pair ($2 \times 8 = 16$)
- In the limit $\omega \to 1$, the breather mass approaches $0$, transitioning to linear small oscillations (mesons)

### 2.5 Summary of Solution Classification

| Solution type | Winding number $N$ | Localization | Oscillation | Mass |
|:--------:|:----------:|:------:|:----:|:----:|
| Kink | $+1$ | Localized | None | $8$ |
| Antikink | $-1$ | Localized | None | $8$ |
| Breather | $0$ | Localized | Yes | $16\sqrt{1-\omega^2}$ |
| Meson (linear wave) | $0$ | Non-localized | Yes | $\geq 2$ (threshold) |

---

## 3. Soliton Collisions and Phase Shifts

### 3.1 Kink-Kink Collision

The collision of two kinks with the same sign ($N = +1$ and $N = +1$) is repulsive, and the two-body solution is given by [2, 12]:

$$\phi_{KK}(x, t) = 4\arctan\left(\frac{v\,\sinh\left(\frac{x}{\sqrt{1-v^2}}\right)}{\cosh\left(\frac{vt}{\sqrt{1-v^2}}\right)}\right)$$

After the collision, the two kinks separate again, perfectly preserving their shape. However, compared to the case without collision, each kink's position is displaced by a **phase shift** $\Delta x$.

$$\Delta x = \frac{\sqrt{1 - v^2}}{v}\ln(v^2)$$

### 3.2 Kink-Antikink Collision

The collision of a kink and an antikink with opposite signs ($N = +1$ and $N = -1$) is attractive, and the two-body solution is given by [2, 9, 12]:

$$\phi_{K\bar{K}}(x, t) = 4\arctan\left(\frac{v^{-1}\,\cosh\left(\frac{vt}{\sqrt{1-v^2}}\right)}{\sinh\left(\frac{x}{\sqrt{1-v^2}}\right)}\right)$$

**Properties:**

- With sufficient energy, the kink and antikink pass through each other and separate again (perfectly elastic collision: velocities and shapes are identical to those before the collision, with only a phase shift)
- With insufficient energy, they form a breather (bound state) that oscillates perpetually without decay

**Note:** In the classical fully integrable sine-Gordon equation, annihilation of a kink-antikink pair into radiation waves (mesons) does not occur. This phenomenon arises only in quantum effects or in systems where integrability is broken (perturbations, discretization, higher dimensions).

### 3.3 Physical Meaning of the Phase Shift

**Proposition 3.1 (Shape preservation and phase shift).** In soliton collisions of the sine-Gordon equation,

1. The **shape** of each soliton after the collision is identical to that before the collision.
2. The **position** of each soliton after the collision is displaced by a phase shift $\Delta x$ compared to the case without collision.
3. The total winding number is conserved before and after the collision.

*Proof.* The sine-Gordon equation can be solved exactly by the inverse scattering method [12], and multi-soliton solutions asymptotically separate into a superposition of individual solitons. The asymptotic shape of each soliton is identical to that before the collision, but a phase shift is added. The winding number is conserved as a topological invariant. $\square$

This property is critically important for the theory. The **type** (winding number = orientation structure) of a soliton (= particle) is invariant under interactions, and only the **position** (phase) changes.

---

## 4. Complete Integrability and the Inverse Scattering Method

### 4.1 Lax Pair

The sine-Gordon equation is a **completely integrable system** and possesses a Lax pair [12, 13].

In light-cone coordinates $u = \frac{1}{2}(x - t)$, $v = \frac{1}{2}(x + t)$, the equation becomes

$$\frac{\partial^2 \phi}{\partial u\,\partial v} = \sin\phi$$

This equation is expressed as the compatibility condition of the following Lax pair.

$$\Psi_u = U\Psi, \quad \Psi_v = V\Psi$$

Here $U$ and $V$ are $2 \times 2$ matrices depending on the spectral parameter $\lambda$ [12].

### 4.2 Infinitely Many Conserved Quantities

As a consequence of complete integrability, the sine-Gordon equation possesses infinitely many independent conserved quantities [12, 14]. The first three are:

- **Energy:**
$$E = \int_{-\infty}^{+\infty}\left[\frac{1}{2}\left(\frac{\partial\phi}{\partial t}\right)^2 + \frac{1}{2}\left(\frac{\partial\phi}{\partial x}\right)^2 + (1-\cos\phi)\right]dx$$

- **Momentum:**
$$P = -\int_{-\infty}^{+\infty}\frac{\partial\phi}{\partial t}\frac{\partial\phi}{\partial x}\,dx$$

- **Winding number:**
$$N = \frac{1}{2\pi}\int_{-\infty}^{+\infty}\frac{\partial\phi}{\partial x}\,dx$$

### 4.3 Bäcklund Transformation

The sine-Gordon equation admits a Bäcklund transformation, which generates new solutions from known ones [7, 15].

Given a known solution $\phi_0$, the solution $\phi_1$ of the following system of equations is also a solution of the sine-Gordon equation.

$$\frac{\partial}{\partial u}\left(\frac{\phi_1 + \phi_0}{2}\right) = a\,\sin\left(\frac{\phi_1 - \phi_0}{2}\right)$$

$$\frac{\partial}{\partial v}\left(\frac{\phi_1 - \phi_0}{2}\right) = \frac{1}{a}\sin\left(\frac{\phi_1 + \phi_0}{2}\right)$$

Here $a$ is the parameter of the Bäcklund transformation. Starting from the vacuum solution $\phi_0 = 0$, a single Bäcklund transformation yields the kink solution, and two transformations yield the kink-kink solution or the breather solution [15].

---

## 5. Extension to Higher Dimensions

### 5.1 2+1 and 3+1 Dimensions

The direct extension of the sine-Gordon equation to higher dimensions is

$$\frac{\partial^2\phi}{\partial t^2} - \nabla^2\phi + \sin(\phi) = 0$$

However, in 2+1 dimensions and above, this equation is no longer completely integrable [16]. Furthermore, by Derrick's theorem [17], no stable localized solutions (solitons) composed solely of scalar fields exist in 3+1 dimensions.

### 5.2 Methods of Stabilization

The following methods are known for obtaining stable solitons in higher dimensions [10, 18].

1. **Skyrme model**: Adding higher-order derivative terms (Skyrme term) [8, 18].
2. **Coupling with gauge fields**: Coupling scalar fields to gauge fields.
3. **Topological stabilization**: When the field carries a nontrivial topological charge, energy minimization prohibits the collapse of the soliton [10].

### 5.3 Stabilization on Closed Spaces

When treating solitons on higher-dimensional spaces, a mechanism for circumventing Derrick's theorem is required. On a closed curve $S^1$ with circumference $L$, the premises of Derrick's theorem (asymptotic conditions on infinite space) are violated, and stable soliton solutions can exist [10].

---

## 6. Physical Applications

### 6.1 Condensed-Matter Physics

- **Josephson junctions**: The motion of magnetic flux quanta across tunnel junctions between superconductors is described by the sine-Gordon equation [19]. Kinks correspond to magnetic flux quanta (fluxons).
- **Crystal dislocations**: Dislocations in crystal lattices in the Frenkel-Kontorova model correspond to kink solutions [1].

### 6.2 Particle Physics

- **Skyrme model**: Baryons are described as topological solitons (Skyrmions) [8, 18]. The baryon number corresponds to the topological charge.
- **Equivalence of the quantum sine-Gordon model and the massive Thirring model**: Coleman [20] and Mandelstam [21] showed that the quantum theory of the sine-Gordon model is equivalent to the massive Thirring model (a fermionic field theory). The **boson-fermion correspondence**, in which bosonic solitons correspond to fermions, is one of the most remarkable results of 1+1-dimensional field theory.

### 6.3 Differential Geometry

- **Pseudospheres**: The compatibility conditions of the first and second fundamental forms of surfaces with constant negative curvature ($K = -1$) reduce to the sine-Gordon equation [5, 6]. Kink solutions correspond to geodesics on pseudospheres.

---

## 7. Implications for Physical Interpretation

### 7.1 Summary of Particle-like Properties of Solitons

Soliton solutions of the sine-Gordon equation naturally possess the following particle-like properties.

| Soliton property | Correspondence to particles |
|:-----------:|:----------------:|
| Winding number $N = \pm 1$ | Particle / antiparticle |
| Conservation of winding number | Conservation of particle number |
| Shape preservation in collisions | Identity of particles |
| Phase shift (position change after collision) | Scattering through interaction |
| Lorentz covariance | Special-relativistic structure |
| Breather ($N = 0$ bound state) | Meson-like bound state |
| Periodicity of potential ($2\pi$) | Discreteness of phase |
| Boson-fermion correspondence [20, 21] | Spin-statistics |

### 7.2 Significance of the Rectangular-Wave Description

As shown in §8–§9, all of these particle-like properties are preserved in the rectangular-wave limit $\sigma \to 0$. This means that the essential information of solitons is not contained in the smooth shape of the continuous field but is concentrated in the topological and scattering discrete structures.

---

## 8. Finite-Parameter Representation of Solitons

### 8.1 Motivation: From Continuous Fields to Finite Parameters

Solutions of the sine-Gordon equation are defined as continuous fields $\phi(x, t)$, but soliton solutions, by virtue of their localization, are **completely characterized by a finite number of parameters**. This is no coincidence. Solitons are **finite-mode solutions** in the mode decomposition of nonlinear fields, restricted to a finite-dimensional invariant manifold (moduli space) within an infinite-dimensional field theory [10, 22].

In this section, we give a formulation describing soliton solutions of the sine-Gordon equation with a finite number of parameters, making explicit how solutions of the continuous wave equation are reduced to a small number of discrete pieces of information.

### 8.2 Five-Parameter Representation of the Kink

We restate the kink solution from §2.1.

$$\phi_K(x, t) = 4\arctan\left(\exp\left(\frac{x - x_0 - vt}{\ell}\right)\right)$$

Here $\ell = \sqrt{1 - v^2}$ is the width parameter of the kink. This solution is **completely** determined by the following five parameters.

**Definition 8.1 (Soliton parameters).** The five parameters characterizing a 1+1-dimensional sine-Gordon kink are defined as follows.

| Symbol | Name | Definition | Range |
|:----:|:----:|:----:|:----:|
| $N$ | Winding number | Topological charge $\pm 1$ | $\{-1, +1\}$ |
| $A$ | Maximum wave height | Total change in $\phi$ $= 2\pi |N|$ | $2\pi$ (fixed) |
| $\theta_0$ | Central phase | Spatial phase of soliton center $= 2\pi x_0 / L$ | $[0, 2\pi)$ |
| $\alpha$ | Translational phase | Phase-angle representation of velocity $v = \sin\alpha$ | $(-\pi/2, \pi/2)$ |
| $\sigma$ | Spread phase | Phase-angle representation of width $\ell = \cos\alpha$ | $(0, 1]$ |

**Proposition 8.1 (Completeness of parameters).** A single kink (antikink) solution of the sine-Gordon equation is uniquely determined by the above five parameters $(N, A, \theta_0, \alpha, \sigma)$. Conversely, the five parameters are uniquely read off from any kink solution.

*Proof.* The general form of the kink solution is $\phi = 4\arctan(\exp(\pm(x - x_0 - vt)/\sqrt{1-v^2}))$, where the sign $\pm$ determines $N$, the displacement $x_0$ determines $\theta_0$, and the velocity $v$ determines $\alpha$. $A = 2\pi$ is fixed by $|N| = 1$, and $\sigma = \cos\alpha$ is determined subordinately from $\alpha$. The free parameters are three: $(N, \theta_0, \alpha)$, with $A$ and $\sigma$ derived from these. $\square$

**Note.** The free parameters are effectively three ($N, \theta_0, \alpha$), but we adopt the five-parameter representation for clarity of physical interpretation. $A$ is determined from the winding number, and $\sigma$ from the translational phase — both are dependent parameters.

### 8.3 Relation Between Translational Phase and Spread Phase

The following identity holds between the translational phase $\alpha$ and the spread phase $\sigma$.

$$\sin^2\alpha + \sigma^2 = v^2 + \ell^2 = 1$$

That is, $(\sin\alpha, \sigma) = (v, \ell)$ is a point on the unit circle.

**Physical interpretation:**

- $\alpha = 0$ (at rest): $\sigma = 1$ — the soliton has maximum width
- $\alpha \to \pm\pi/2$ (speed of light): $\sigma \to 0$ — the soliton width contracts to zero (Lorentz contraction)
- Faster motion yields narrower spatial extent — **motion and shape are unified by a single phase angle $\alpha$**

### 8.4 Six-Parameter Representation of the Breather

The breather solution of §2.4 requires the internal oscillation frequency $\omega$ in addition to the five kink parameters.

**Definition 8.2 (Breather parameters).**

| Symbol | Name | Definition | Range |
|:----:|:----:|:----:|:----:|
| $N$ | Winding number | $0$ (kink-antikink bound) | $\{0\}$ (fixed) |
| $A$ | Maximum wave height | $4\arctan(\sqrt{1-\omega^2}/\omega)$ | $(0, \pi)$ |
| $\theta_0$ | Central phase | Spatial phase of breather center | $[0, 2\pi)$ |
| $\alpha$ | Translational phase | Bulk motion $v = \sin\alpha$ | $(-\pi/2, \pi/2)$ |
| $\sigma$ | Spread phase | $1/\sqrt{1-\omega^2}$ | $[1, \infty)$ |
| $\omega$ | Internal phase velocity | Angular frequency of internal oscillation | $(0, 1)$ |

The breather has $N = 0$ but possesses an additional degree of freedom in the internal oscillation $\omega$.

### 8.5 Parameter Representation of Multi-Soliton States

A multi-body state consisting of $n$ solitons is represented as a collection of parameters for each soliton [10, 12].

**Definition 8.3 ($n$-soliton state).** A state consisting of $n$ solitons is defined by the following parameter set.

$$\mathcal{S}_n = \{(N_i, A_i, \theta_{0,i}, \alpha_i, \sigma_i)\}_{i=1}^{n}$$

**Proposition 8.2 (Asymptotic completeness).** A field configuration $\phi(x, t)$ consisting of $n$ well-separated solitons is asymptotically completely determined by the parameter set $\mathcal{S}_n$. In regions where the overlap between solitons is negligible,

$$\phi(x, t) \approx \sum_{i=1}^{n} \phi_i(x, t; N_i, \theta_{0,i}, \alpha_i)$$

holds, where $\phi_i$ is the single-soliton solution for each soliton.

*Proof.* By the inverse scattering method [12], the $n$-soliton solution is completely determined by $n$ scattering data (eigenvalues and norming constants). In the asymptotic region, the scattering data correspond one-to-one with the above parameters. $\square$

### 8.6 Description of Collisions in Parameter Space

The soliton collisions of §3 can be described as mappings in parameter space.

**Proposition 8.3 (Collision map).** The collision of two solitons $\mathcal{S}_1 = (N_1, \theta_{0,1}, \alpha_1)$ and $\mathcal{S}_2 = (N_2, \theta_{0,2}, \alpha_2)$ is described as the following parameter transformation.

$$\mathcal{S}_1' = (N_1,\; \theta_{0,1} + \Delta\theta_1,\; \alpha_1)$$
$$\mathcal{S}_2' = (N_2,\; \theta_{0,2} + \Delta\theta_2,\; \alpha_2)$$

That is:

1. **Winding number $N_i$ is invariant** — the particle type does not change
2. **Translational phase $\alpha_i$ is invariant** — velocity (momentum) is conserved
3. **Only the central phase $\theta_{0,i}$ changes** — displaced by the phase shift $\Delta\theta_i$

**Corollary 8.1.** All information about a soliton collision is concentrated in the phase shifts $\Delta\theta_i$. The phase shifts are computable from the incident parameters [12]:

$$\Delta\theta_1 = -\Delta\theta_2 = f(N_1, N_2, \alpha_1 - \alpha_2)$$

That is, the collision depends only on the relative translational phase $\alpha_1 - \alpha_2$.

### 8.7 Reduction from Continuous Fields to Parameter Representation

Summarizing the above, the solution space of the sine-Gordon equation is hierarchically structured as follows.

| Description level | Degrees of freedom | Information content |
|:----------:|:------:|:------:|
| Continuous field $\phi(x, t)$ | Infinite-dimensional | Real value at each spacetime point |
| $n$-soliton solution | $3n$-dimensional | Each soliton $(N_i, \theta_{0,i}, \alpha_i)$ |
| Post-collision change | $n$-dimensional | Phase shifts $\Delta\theta_i$ only |

**Proposition 8.4 (Finite-dimensionality of the soliton sector).** The sector of the solution space of the sine-Gordon equation containing $n$ solitons is described as a $3n$-dimensional finite-dimensional manifold (moduli space) $\mathcal{M}_n$ [10, 22]. The infinite degrees of freedom of the continuous field are **exactly** reduced to a finite number of parameters in the soliton sector.

*Proof.* In the inverse scattering method, the $n$-soliton sector is characterized by $n$ discrete eigenvalues $\zeta_i$ and associated norming constants $c_i$ [12, 13]. Each eigenvalue determines the winding number $N_i$ and velocity $v_i$ (= $\sin\alpha_i$), while the norming constant determines the initial position $x_{0,i}$ (= $\theta_{0,i}$). Excluding the radiation component (continuous spectrum), the degrees of freedom of the soliton sector are exactly $3n$. $\square$

This reduction is not an approximation. In the soliton sector, the infinite-dimensional partial differential equation is **exactly** reduced to a finite-dimensional system of ordinary differential equations.

### 8.8 Relation to Rectangular Waves

We re-examine the shape of the kink solution.

$$\phi_K(x) = 4\arctan(e^{x/\ell}) = 2\pi \cdot \frac{1}{1 + e^{-x/\ell}}$$

The right-hand side is the logistic function (sigmoid function), which converges to a **step function** (half-period of a rectangular wave) in the limit $\ell \to 0$.

$$\lim_{\ell \to 0} \phi_K(x) = \begin{cases} 0 & (x < 0) \\ 2\pi & (x > 0) \end{cases}$$

**Proposition 8.5 (Rectangular-wave limit).** In the limit $\sigma = \ell \to 0$ of the spread parameter (i.e., $\alpha \to \pm\pi/2$, the speed-of-light limit), the kink solution converges to a rectangular wave (step function). Conversely, a kink with finite width $\sigma > 0$ can be understood as a **smoothed** rectangular wave.

*Proof.* $\phi_K(x) = 2\pi/(1 + e^{-x/\ell})$ is a logistic function, and it converges pointwise to the Heaviside step function $2\pi \cdot H(x)$ as $\ell \to 0$. $\square$

**Corollary 8.2.** A soliton can be interpreted as a finite-mode approximation of a rectangular wave. The continuous $\tanh$-type shape is merely a smoothed version of a discrete step function. The spread parameter $\sigma$ controls the degree of smoothing.

From this perspective, the parameter representation of solitons is reinterpreted as follows.

| Parameter | Meaning in continuous field | Meaning in rectangular wave (discrete) |
|:----------:|:------------:|:--------------------:|
| $N = \pm 1$ | Winding number | Direction of step (ascending/descending) |
| $A = 2\pi$ | Total phase change | Height of step |
| $\theta_0$ | Position of soliton center | Position of step edge |
| $\alpha$ | Phase angle of velocity | Direction of edge motion |
| $\sigma$ | Width of smoothing | **Coincides with rectangular wave at $\sigma \to 0$** |

---

## 9. Preservation of Properties in the Rectangular-Wave Limit

### 9.1 Purpose

The properties of the sine-Gordon equation described in §1–§7 have been proved for the continuous field $\phi(x, t)$. As shown in §8.8, soliton solutions converge to rectangular waves (step functions) in the limit $\sigma \to 0$ of the spread parameter.

In this section, we systematically verify for each property in §1–§7 that (a) it holds for soliton solutions (this is trivial since they are exact solutions), and (b) it is preserved in the rectangular-wave limit $\sigma \to 0$.

If all properties are preserved in the rectangular-wave limit, then the rectangular-wave description is not a "convenient approximation" but an **exact limit that completely inherits** the structure of solitons.

### 9.2 Framework of Verification

Each property is verified in three stages.

1. **Validity for soliton solutions**: Does it satisfy the sine-Gordon equation as an exact solution?
2. **Analysis of $\sigma$-dependence**: Does the property depend on the spread parameter $\sigma$?
3. **Existence of the rectangular-wave limit**: Is the property preserved or does it diverge as $\sigma \to 0$?

### 9.3 Definition and Conservation of Winding Number (§2.3 → Rectangular Wave)

**Property from §2.3:** The winding number $N = [\phi(+\infty) - \phi(-\infty)] / 2\pi$ is conserved under time evolution.

**Validity for soliton solutions:** The kink solution gives $\phi(-\infty) = 0$, $\phi(+\infty) = 2\pi$, yielding $N = +1$. Since this is an exact solution of the sine-Gordon equation, conservation of $N$ is trivial. $\checkmark$

**$\sigma$-dependence:** The definition $N = [\phi(+\infty) - \phi(-\infty)] / 2\pi$ depends **only on boundary values** and is entirely independent of the internal shape of $\phi$ (and hence of $\sigma$).

**Proposition 9.1 ($\sigma$-independence of winding number).** The winding number $N$ does not depend on the spread parameter $\sigma$. It retains the same integer value in the rectangular-wave limit $\sigma \to 0$.

*Proof.* For any $\sigma > 0$, $\phi_K^\sigma(\pm\infty)$ is fixed at the vacuum values $0, 2\pi$. In the limit $\sigma \to 0$, $\phi_K^\sigma(x) \to 2\pi \cdot H(x)$ pointwise, and since $H(-\infty) = 0$, $H(+\infty) = 1$, we have $N = (2\pi \cdot 1 - 2\pi \cdot 0)/2\pi = 1$, unchanged. $\square$

**Winding number for rectangular waves:** For the rectangular wave $\phi(x) = 2\pi H(x - x_0)$, $N = +1$. If the step is ascending, $N = +1$; if descending, $N = -1$. The winding number is precisely the sign of the step. $\checkmark$

### 9.4 Conservation Law of Winding Number (Proposition 2.1 → Rectangular Wave)

**Property from §2.3:** $N$ is conserved under time evolution. Creation and annihilation are possible only as kink-antikink pairs.

**$\sigma$-dependence:** The proof of the conservation law (§2.3) relies solely on the fact that "boundary values are fixed at vacuum values and cannot change by an integer multiple of $2\pi$ under continuous time evolution." This argument is independent of $\sigma$.

**Conservation for rectangular waves:** The boundary values of a rectangular wave are $0$ and $2\pi$, taking only discrete values. Discrete values cannot change under continuous deformation. Therefore, conservation of the winding number for rectangular waves is **even more trivially** satisfied than for continuous fields.

**Proposition 9.2.** The conservation law of the winding number holds exactly in the rectangular-wave limit. Moreover, since $\phi$ itself takes only discrete values ($0$ or $2\pi$) for rectangular waves, conservation is guaranteed not only topologically but also by the **discreteness of values**. $\checkmark$

### 9.5 Existence of Kink and Antikink Solutions (§2.1, §2.2 → Rectangular Wave)

**Property from §2.1–2.2:** Localized solutions with winding number $\pm 1$ exist and propagate with velocity $v$.

**Validity for soliton solutions:** The kink $\phi_K = 4\arctan(\exp((x-vt)/\sigma))$ is an exact solution. $\checkmark$

**Rectangular-wave limit:** As $\sigma \to 0$, $\phi_K \to 2\pi \cdot H(x - vt)$. This is a step function moving with velocity $v$, preserving localization and propagation.

**Proposition 9.3 (Existence of rectangular-wave kink).** The rectangular-wave kink $\phi(x, t) = 2\pi \cdot H(x - vt)$ has the following properties.
1. Winding number $N = +1$
2. Moves with velocity $v$
3. Localized at a single spatial point ($x = vt$)
4. Holds as a **distributional solution** of the sine-Gordon equation

*Proof.* (1)–(3) are immediate from the definition. For (4): the distributional derivatives of $\phi = 2\pi H(x - vt)$ are $\phi_x = 2\pi\delta(x-vt)$, $\phi_{xx} = 2\pi\delta'(x-vt)$, $\phi_{tt} = 2\pi v^2\delta'(x-vt)$, and $\sin(\phi) = \sin(2\pi H(x-vt)) = 0$ (since $\phi = 0$ or $2\pi$ for $x \neq vt$, and $\sin = 0$ in both cases). The wave-equation part $\phi_{tt} - \phi_{xx} = 2\pi(v^2-1)\delta'(x-vt)$ vanishes only for $v = \pm 1$ (speed of light). For $v \neq 1$, the rectangular wave is not a classical solution, and smoothing with finite width $\sigma > 0$ is required. $\square$

**Note.** The fact that the rectangular wave is not an exact solution for $v \neq 1$ is consistent with the relation $\sigma = \cos\alpha = \sqrt{1-v^2}$ from §8.3. If $v < 1$ then $\sigma > 0$, and the soliton has finite width. A rectangular wave ($\sigma = 0$) corresponds only to $v = 1$ (speed of light). The fact that subluminal particles have finite "extent" is qualitatively consistent with the spread of wave packets in quantum mechanics.

### 9.6 Breather Solution (§2.4 → Rectangular Wave)

**Property from §2.4:** A breather exists as a bound state of a kink-antikink pair.

**Rectangular-wave limit:** The breather is an oscillatory solution with winding number $N = 0$. In the rectangular-wave limit, it corresponds to a **rectangular pulse** oscillating between $0$ and $2\pi$.

**Proposition 9.4.** The rectangular-wave limit of the breather is a spatially localized rectangular pulse oscillating between the two states of presence ($\phi = 2\pi$) and absence ($\phi = 0$). The mass $M_B = 16\sqrt{1-\omega^2}$ depends on $\omega$ and is not directly dependent on $\sigma$. $\checkmark$

### 9.7 Soliton Collisions and Phase Shift (§3 → Rectangular Wave)

**Property from §3:** After collision, (i) shape preservation, (ii) phase shift, (iii) conservation of winding number.

This is the most important property, and we verify each item individually.

**(i) Shape preservation:**

**Soliton solutions:** Each soliton after the collision has the same shape (same $\sigma$) as before [12]. $\checkmark$

**Rectangular-wave limit:** A rectangular wave has $\sigma = 0$, and the only shape parameter is "step height $2\pi$." If $\sigma = 0$ is preserved before and after the collision, shape preservation holds automatically. Since $\sigma$ is a conserved quantity in soliton collisions (Proposition 8.3: $\alpha_i$ invariant $\Rightarrow$ $\sigma_i = \cos\alpha_i$ invariant), $\sigma = 0$ is also preserved. $\checkmark$

**(ii) Phase shift:**

**Soliton solutions:** The phase shift after collision is $\Delta\theta = f(N_1, N_2, \alpha_1 - \alpha_2)$ [12]. $\checkmark$

**Rectangular-wave limit:** The phase-shift formula depends only on $\alpha_i$ (translational phase) and not on $\sigma_i$ (spread) (see §8.6 Corollary 8.1). Therefore, the phase-shift value does not change even as $\sigma \to 0$.

**Proposition 9.5 ($\sigma$-independence of phase shift).** The phase shift $\Delta\theta$ in a soliton collision does not depend on the spread parameters $\sigma_i$ of the individual solitons. In the rectangular-wave limit $\sigma \to 0$, the phase shift takes the same value as for solitons.

*Proof.* In the inverse scattering method, the phase shift is computed solely from the discrete eigenvalues $\zeta_i$ of the scattering data [12, 13]. $\zeta_i$ is determined by the velocity $v_i = \sin\alpha_i$ and winding number $N_i$ of the soliton, independent of the norming constants (position information) and the radiation component. The soliton width $\sigma_i$ is related to the norming constants but does not contribute to the phase shift. $\square$

**(iii) Conservation of winding number:** Already shown in §9.4. $\checkmark$

### 9.8 Complete Integrability and Conserved Quantities (§4 → Rectangular Wave)

**Property from §4:** Conservation of energy, momentum, and winding number.

**Energy:**

The energy of a soliton is $E = 8/\sigma = 8/\cos\alpha = 8/\sqrt{1-v^2}$. As $\sigma \to 0$, the energy diverges.

**Proposition 9.6 ($\sigma$-dependence of energy).** The energy of a soliton depends on $\sigma$, diverging as $E \to \infty$ when $\sigma \to 0$. This means that the rectangular-wave limit has infinite energy.

**Proposition 9.6a (Limitations of the theory).** The sine-Gordon soliton theory is based on the following implicit assumptions.

1. **Fixed starting point**: The initial value of the field is fixed at the vacuum value $\phi(-\infty) = 0$ (§2.1).
2. **Periodicity of the potential**: $V(\phi) = 1-\cos\phi$ has period $2\pi$, and the target space is effectively $S^1$. The size of one "period" of the potential ($2\pi$) is fixed by the theory and has no freedom to be changed.
3. **Forcing the terminus to a vacuum value**: The requirement of finite energy on infinite space $(-\infty, +\infty)$ forces $\phi(+\infty) = 2n\pi$ ($n \in \mathbb{Z}$). This is because $V(2n\pi) = 0$ (potential minimum), while $\phi = (2n+1)\pi$ gives $V = 2$ (potential maximum), and the energy diverges on infinite space, excluding such configurations.
4. **Exclusion of fractional winding numbers**: As a consequence of the above, the winding number $N$ is always an integer. Half-integer winding numbers ($\phi(+\infty) = (2n+1)\pi$) are automatically excluded by the finite-energy requirement.
5. **Absence of boundary reflection**: Since the space is infinite or closed ($S^1$), no physical boundaries exist, and boundary reflection of solitons does not occur.
6. **Perfectly elastic collisions and absence of interaction**: The theory is a completely integrable system (§4) with infinitely many conserved quantities. As a consequence of this excessive constraint, soliton collisions are described as perfectly elastic collisions, in which the velocity, shape, and energy of each soliton are **individually** conserved before and after the collision (§3.3). To avoid violating momentum conservation, only the position (phase shift $\Delta x$) of each soliton changes, and no exchange of velocities, deformation of shapes, or dissipation of energy occurs. This effectively means that solitons behave as **free particles that do not interact**. In the case of identical particles (e.g., kink-kink collision), it is in principle indistinguishable whether they "passed through each other" or "bounced off."

*Consequence.* This theory can be consistently applied only in the following two cases.

- **Infinite space** $(-\infty, +\infty)$: The terminus is at infinity and is never reached. The energy per kink is $E = 8$ (at rest), which is finite.
- **Closed curve** $S^1$ (circumference $L$): Periodic boundary conditions $\phi(x+L) = \phi(x) + 2\pi N$. No boundary exists.

**Open strings** (Dirichlet/Neumann boundary conditions on a finite interval $[0, L]$) fall outside the scope of this theory. The reason is that a finite interval $[0, L]$ is a contractible space ($\pi_1([0,L]) = 0$), so the winding number is not defined as a topological invariant, and the topological stability of solitons is not guaranteed. Furthermore, on a finite interval, $\phi(L) = (2n+1)\pi$ (half-integer winding number, potential maximum) is energetically admissible, so the restriction to integer winding numbers breaks down.

**Self-cancellation structure of energy:** The energy of each $2\pi$ rotation (one kink) is $E = 8$ (at rest) and is self-contained. This is because the potential $V(\phi) = 1-\cos\phi$ returns from $0$ to $0$ over one period. A kink-antikink pair ($4\pi$ rotation, winding number $N = 0$) can annihilate to return the energy to zero. The only observable persistent energy is the contribution $8|N|$ from the $|N|$ surplus kinks that cannot pair-annihilate.

**Momentum and winding number:** Momentum $P = -Ev$ has the same structure. The winding number is $\sigma$-independent as shown in §9.3. $\checkmark$

### 9.9 Extension to Higher Dimensions (§5 → Rectangular Wave)

**Property from §5:** By Derrick's theorem, no stable localized solutions exist in 3+1-dimensional infinite space.

**Impact on the rectangular-wave limit:** The premise of Derrick's theorem is the existence of **finite-energy solutions on infinite space** [17]. On a closed curve $S^1$ (circumference $L$), the premises of the theorem are violated for the following two reasons.

1. **Space is closed**: The range of $\lambda$ in the scale transformation $\phi(x) \to \phi(\lambda x)$ is restricted
2. **Topological protection**: The winding number $N \in \pi_1(S^1) = \mathbb{Z}$ is protected as a topological invariant

**Proposition 9.7.** On a closed curve $S^1$ (circumference $L$), localized solutions of any width, including rectangular waves ($\sigma = 0$), are topologically stable. Derrick's theorem does not apply to closed spaces. $\checkmark$

### 9.10 Summary of Verification Results

| Property | Section | Soliton | $\sigma$-dependence | Rectangular-wave limit | Remarks |
|:----:|:--:|:-------:|:--------------:|:----------:|:----:|
| Winding number definition | §2.3 | $\checkmark$ | **Independent** | $\checkmark$ | Boundary values only |
| Winding number conservation | §2.3 | $\checkmark$ | **Independent** | $\checkmark$ | Topological |
| Kink/antikink existence | §2.1–2 | $\checkmark$ | Shape only | $\checkmark$ | Distributional solution |
| Breather existence | §2.4 | $\checkmark$ | Shape only | $\checkmark$ | Rectangular pulse |
| Shape preservation after collision | §3.1–2 | $\checkmark$ | **Independent** | $\checkmark$ | $\sigma$ conserved |
| Phase shift after collision | §3.3 | $\checkmark$ | **Independent** | $\checkmark$ | $\alpha$ only |
| Winding number conservation (collision) | §3.3 | $\checkmark$ | **Independent** | $\checkmark$ | Topological |
| Energy conservation | §4.2 | $\checkmark$ | Dependent | $\checkmark$* | *Self-contained per $2\pi$ rotation |
| Momentum conservation | §4.2 | $\checkmark$ | Dependent | $\checkmark$* | *Same as above |
| Derrick avoidance | §5.1–2 | — | — | $\checkmark$ | Closed space $S^1$ |
| Parameter representation $(N, \theta_0, \alpha)$ | §8 | $\checkmark$ | **Independent** | $\checkmark$ | $\sigma$ is subordinate |

**Theorem 9.1 (Structural preservation in the rectangular-wave limit).** The following structural properties of soliton solutions of the sine-Gordon equation are all preserved in the rectangular-wave limit $\sigma \to 0$.

1. **Definition and conservation of winding number** (topological property: $\sigma$-independent)
2. **Shape preservation and phase shift after collision** (property of scattering data: $\sigma$-independent)
3. **Finite-parameter representation** $(N, \theta_0, \alpha)$ (structure of moduli space: $\sigma$-independent)
4. **Finiteness of energy and momentum** (by the self-cancellation structure of each $2\pi$ rotation. Proposition 9.6a)

*Proof.* (1) follows from Propositions 9.1–9.2, (2) from Proposition 9.5, (3) from Proposition 8.4 and the subordination $\sigma = \cos\alpha$, (4) from the energy self-cancellation structure in Proposition 9.6a. $\square$

**Corollary 9.1.** The rectangular-wave description is not a "crude approximation" of the soliton. It is an **exact limit that completely inherits** the topological and scattering structures of the soliton. The only information lost is the spread parameter $\sigma$ (details of the internal shape).

**Corollary 9.2 (Scope of applicability).** The validity of this theorem depends on the implicit assumptions of the theory made explicit in Proposition 9.6a — the fixed starting point $\phi = 0$, the forcing of the terminus to a vacuum value ($2n\pi$), and integer winding numbers. These assumptions hold only on infinite space $(-\infty,+\infty)$ or on a closed curve $S^1$. Extension to open strings (boundary conditions on a finite interval) requires admitting half-integer winding numbers and a theory of boundary reflection (Ghoshal-Zamolodchikov [23]), which is beyond the scope of this paper.

---

## 10. Conclusion

In this paper, we surveyed the mathematical structure of the sine-Gordon equation and demonstrated the following.

1. Soliton solutions (kink, antikink, breather) are **exactly** characterized by a finite number of parameters $(N, \theta_0, \alpha)$ (§8).

2. All structural properties of solitons — conservation of winding number, shape preservation and phase shift in collisions, finite-parameter representation — are preserved in the rectangular-wave limit $\sigma \to 0$ (§9).

3. The finiteness of energy is based on the structure in which the energy of each $2\pi$ rotation is self-contained (the potential $V(\phi) = 1-\cos\phi$ returns from $0$ to $0$ over one period). Kink-antikink pairs ($4\pi$ rotation) can pair-annihilate, and persistent energy consists only of the contribution from surplus kinks that cannot pair-annihilate.

4. Therefore, the transition from solitons to rectangular waves entails no loss of information. Among the infinite degrees of freedom of the continuous field, the physically meaningful information is concentrated in the finite number of parameters $(N, \theta_0, \alpha)$, and this concentration holds regardless of the soliton shape.

**Limitations of applicability:** The above conclusions hold under the implicit assumptions of the theory (Proposition 9.6a).

- **Spatial constraints**: The phase at the starting point is fixed at the vacuum value $\phi = 0$, and the phase at the terminus is forced to a vacuum value $\phi = 2n\pi$ (winding number $N$ is integer only). No boundary reflection exists. These conditions are satisfied only on infinite space $(-\infty, +\infty)$ or a closed curve $S^1$. Extension to open strings (boundary conditions on a finite interval) requires admitting half-integer winding numbers ($\phi = (2n+1)\pi$, potential maximum) and handling soliton reflection at boundaries, which is beyond the scope of this paper.
- **Collision constraints**: This theory is a completely integrable system, and soliton collisions are described as perfectly elastic collisions. To avoid violating momentum conservation, the velocity, shape, and energy of each soliton are individually conserved before and after the collision, and only the position (phase shift $\Delta x$) changes. This effectively means that solitons behave as free particles that do not interact, and dissipation of energy, exchange of velocities, and inelastic channels (pair-annihilation → radiation) do not occur within the classical theory.

The rectangular-wave description is not a "crude approximation" of the soliton but an exact limit that completely inherits the topological and scattering structures. The fact that solutions of nonlinear wave equations are reduced to a finite number of discrete parameters reveals a deep correspondence between the theory of continuous fields and discrete structures.

An alternative model that overcomes the above limitations — in particular, the inability to handle particles of different masses and the absence of velocity redistribution after collisions — is presented in §11.

---

## 11. Perfectly Elastic Collision Model at Arbitrary Phases in a Perfect Fluid

In §9–10, we identified the limitations of the sine-Gordon soliton theory. In particular, (i) all solitons have the same mass, and (ii) no velocity redistribution occurs after collisions (transmissive collisions), which are insufficient for describing physical elastic collisions.

In this section, we formulate an alternative model without these constraints: **perfectly elastic collisions at arbitrary phases in a perfect fluid with different specific gravities (masses)**, within the framework of wave equations.

This model is based on the wave equation for longitudinal waves (sound waves) in a perfect fluid. Waves in a perfect fluid are inherently longitudinal, and the propagation velocity depends on the density (specific gravity) of the medium ($c = \sqrt{K/\rho}$). The collision of two waves with different specific gravities is completely described by momentum conservation and energy conservation. Below, we formulate this physics in terms of a spring-mass system, which provides the most transparent representation.

### 11.1 Model Construction

Two perfect springs (zero dissipation) with spring constant $k$ and natural length $l$ are connected in series, with both ends fixed to walls separated by a distance $2l$. A point mass $M$ is placed at the junction of the two springs.

$$\text{Wall}—[\,k,\; l\,]—M—[\,k,\; l\,]—\text{Wall}$$

Letting $x$ denote the displacement from the equilibrium position (center), the restoring forces from the left spring $-kx$ and the right spring $-kx$ combine to give the equation of motion

$$M\ddot{x} + 2kx = 0 \tag{W1}$$

This is simple harmonic motion (longitudinal wave) with angular frequency $\omega = \sqrt{2k/M}$. The general solution is

$$x(t) = A\cos(\omega t + \phi_0) \tag{W2}$$

where $A$ is the amplitude and $\phi_0$ is the initial phase.

### 11.2 Two-Body Configuration and Definition of Collision

We consider two independent simple harmonic oscillators with different masses $M_1, M_2$.

$$\ddot{x}_1 + \omega_1^2 x_1 = 0, \qquad \omega_1 = \sqrt{\frac{2k}{M_1}} \tag{W3}$$

$$\ddot{x}_2 + \omega_2^2 x_2 = 0, \qquad \omega_2 = \sqrt{\frac{2k}{M_2}} \tag{W4}$$

When $M_1 \neq M_2$, $\omega_1 \neq \omega_2$, and the two systems oscillate independently. The difference in specific gravity of the perfect fluid corresponds to this mass difference $M_1 \neq M_2$.

**Definition of collision.** In this paper, a "collision" is defined as a mathematical process in which the states $(x_1, \dot{x}_1)$ and $(x_2, \dot{x}_2)$ of the two systems interact via momentum conservation and energy conservation at a given time $t = t_c$, instantaneously transitioning to new states $(x_1, \dot{x}_1')$ and $(x_2, \dot{x}_2')$. Positions $x_i$ do not change due to the collision; only velocities $v_i = \dot{x}_i$ change according to (W7)–(W8).

**Remark.** In this model, a collision is not a geometric description of physical contact but an abstraction: an "instantaneous state transition in which the states of two oscillatory systems are coupled via conservation laws." It replaces the role of the physical contact surface in longitudinal wave collisions in a perfect fluid with a state transformation via conservation laws. The collision time $t_c$ is an input parameter given externally and is not determined from within the model. This allows collisions at arbitrary phases $\phi_1, \phi_2$ to be treated in a unified manner.

### 11.3 Perfectly Elastic Collision at Arbitrary Phases

**Pre-collision state:** The state of each system at time $t_c$ is expressed as follows.

$$x_{1,c} = A_1\cos\phi_1, \qquad v_1 = -A_1\omega_1\sin\phi_1 \tag{W5}$$

$$x_{2,c} = A_2\cos\phi_2, \qquad v_2 = -A_2\omega_2\sin\phi_2 \tag{W6}$$

Here $\phi_i = \omega_i t_c + \phi_{0,i}$ is the phase of each system at the collision instant, and $A_i$ is the amplitude. The sign of the velocity follows from the natural differentiation of (W3)(W4).

**Conditions for perfectly elastic collision (momentum conservation + energy conservation):**

$$M_1 v_1 + M_2 v_2 = M_1 v_1' + M_2 v_2' \tag{W7}$$

$$\frac{1}{2}M_1 v_1^2 + \frac{1}{2}M_2 v_2^2 = \frac{1}{2}M_1 v_1'^2 + \frac{1}{2}M_2 v_2'^2 \tag{W8}$$

These two conditions uniquely determine the post-collision velocities.

**Theorem 11.1 (Elastic collision formula at arbitrary phases).** The post-collision velocities of a perfectly elastic collision between masses $M_1, M_2$ at phases $\phi_1, \phi_2$ are given by:

$$v_1' = \frac{(M_1 - M_2)\,v_1 + 2M_2\,v_2}{M_1 + M_2} \tag{W9}$$

$$v_2' = \frac{(M_2 - M_1)\,v_2 + 2M_1\,v_1}{M_1 + M_2} \tag{W10}$$

*Proof.* Obtained directly by solving the system of equations (W7) and (W8). $\square$

### 11.4 Post-Collision Oscillation State

The collision is instantaneous, and positions do not change. After the collision, each mass resumes oscillation in its original spring system with the collision-time position $x_{i,c}$ and the new velocity $v_i'$ as initial conditions.

$$x_1(t) = x_{1,c}\cos(\omega_1\tau) + \frac{v_1'}{\omega_1}\sin(\omega_1\tau) \qquad (\tau = t - t_c) \tag{W11}$$

$$x_2(t) = x_{2,c}\cos(\omega_2\tau) + \frac{v_2'}{\omega_2}\sin(\omega_2\tau) \tag{W12}$$

**Post-collision amplitude and phase:**

$$A_i' = \sqrt{x_{i,c}^2 + \left(\frac{v_i'}{\omega_i}\right)^2} \tag{W13}$$

$$\phi_i' = \arctan\!\left(\frac{-v_i'}{\omega_i\, x_{i,c}}\right) \tag{W14}$$

The angular frequency $\omega_i$ is determined solely by the spring constant and mass, and is not changed by the collision. What changes are only the amplitude $A_i$ and the phase $\phi_i$.

### 11.5 Proof of Energy Conservation

**Proposition 11.1.** In a perfectly elastic collision at arbitrary phases, the total energy of the system is conserved.

*Proof.* The total energy of each system is the sum of the spring's potential energy and kinetic energy.

$$E_i = \frac{1}{2}M_i\omega_i^2\,x_{i,c}^2 + \frac{1}{2}M_i\,v_i^2 = k\,x_{i,c}^2 + \frac{1}{2}M_i\,v_i^2$$

Since the collision is instantaneous, the potential energy $k\,x_{i,c}^2$ is unchanged before and after the collision. The kinetic energy is conserved by the elastic collision condition (W8). Therefore

$$E_1 + E_2 = E_1' + E_2' \qquad \square \tag{W15}$$

### 11.6 Characteristics of Amplitude and Phase Changes by Collision Phase

From the elastic collision formulas (W9)–(W10) and the post-collision amplitude/phase determination formulas (W19)(W14), the following characteristics emerge depending on the collision-instant phase $\phi_i$.

| Collision phase | Position $x_{i,c}$ | Velocity $|v_i|$ | Energy partition | Collision effect |
|:--|:--|:--|:--|:--|
| $\phi_i = \pi/2$ | $0$ | $A_i\omega_i$ (maximum) | All kinetic energy | Velocity is fully redistributed; only amplitude $A_i$ changes |
| $\phi_i = 0$ | $A_i$ (maximum) | $0$ | All potential energy | No velocity exchange occurs; state is unchanged before and after |
| $\phi_i = \pi/4$ | $A_i/\sqrt{2}$ | $A_i\omega_i/\sqrt{2}$ | Equally partitioned | Partial velocity redistribution and phase shift occur |

**Remark.** The "unchanged state" at $\phi_i = 0$ arises because the right-hand sides of (W9)(W10) are determined solely by $v_i$, and $v_i = 0$ yields $v_i' = 0$ as the solution. This does not mean that no collision occurs, but rather that the states with zero velocity, upon being coupled via conservation laws, remain invariant. In the general case $\phi_1 \neq \phi_2$, a nontrivial state transformation necessarily occurs as long as both velocities are not simultaneously zero.

### 11.7 Comparison with the sine-Gordon Theory

| Property | sine-Gordon (§1–10) | This model (§11) |
|:--|:--|:--|
| Different masses | Impossible (all kinks have identical mass) | **Possible** ($M_1 \neq M_2$) |
| Post-collision velocity redistribution | None (transmissive) | **Present** (elastic collision formula) |
| Collision at arbitrary phases | Undefined | **Fully defined** |
| Closed-string assumption | Required (implicit) | **Not required** |
| Integer winding number constraint | Required | **Not required** |
| Half-integer phase | Excluded ($\phi = \pi$ inadmissible) | **Admissible** |
| Wave equation | $\phi_{tt} - \phi_{xx} + \sin\phi = 0$ | $\ddot{x}_i + \omega_i^2 x_i = 0$ |
| Hidden assumptions | 6 (Proposition 9.6a) | **None** (conservation laws only) |
| Energy dissipation | Zero | Zero |
| Prerequisites of formulation | Momentum conservation + energy conservation + topological constraints + integrability + closed space + identical masses | **Momentum conservation + energy conservation only** |

### 11.8 Physical Interpretation

The physical correspondences of the parameters in this model are shown below.

| Model parameter | Physical correspondence |
|:--|:--|
| Mass $M_i$ | Specific gravity (density) of the perfect fluid |
| Amplitude $A_i$ | Wave amplitude (magnitude of energy) |
| Phase $\phi_i$ | Wave phase |
| Angular frequency $\omega_i = \sqrt{2k/M_i}$ | Natural frequency (depends on specific gravity) |
| Collision (state transition) | Instantaneous interaction of waves via conservation laws |

The process in which longitudinal waves in a perfect fluid undergo perfectly elastic collision is completely described by (W1)–(W15). The only assumptions are the two conservation laws of momentum and energy; no additional assumptions of topological constraints, integrability, compactness of space, or integer winding numbers are required.

### 11.9 Equivalence with the Rectangular-Wave Limit

In §9, we showed that the solitons of the sine-Gordon equation preserve their structure in the rectangular-wave limit. In this section, we prove that for the model of §11 as well, the collision results remain identical when the waveform changes from a sinusoidal (cosine) wave to a rectangular wave.

**Introduction of the waveform parameter.** We define a family $x^\sigma(t)$ of oscillation waveforms controlled by the transition-width parameter $\sigma$ ($0 < \sigma \leq 1$).

- $\sigma = 1$: Smooth sinusoidal (cosine) wave $x^1(t) = A\cos(\omega t)$
- $\sigma \to 0$: Rectangular wave $x^0(t) = A\,\mathrm{sgn}(\cos(\omega t))$

Intermediate $\sigma$ controls the steepness of the transition region. For all $\sigma$, the amplitude is $A$ and the period is $T = 2\pi/\omega$.

**Theorem 11.2 (Waveform independence of collision results).** Given the collision-time positions $x_{i,c}$ and total energies $E_i$, the post-collision velocities $v_i'$, amplitudes $A_i'$, and phases $\phi_i'$ do not depend on the waveform parameter $\sigma$.

*Proof.* We show this in three stages.

**Stage 1: The collision-time velocity is uniquely determined by energy and position.**

The total energy of the spring-mass system is conserved regardless of the waveform.

$$E_i = k\,x_i^2(t) + \frac{1}{2}M_i\,\dot{x}_i^2(t) = k\,A_i^2 \tag{W16}$$

This is independent of the waveform. The spring's potential energy is a function of position only, and the kinetic energy is the remainder. Therefore, the velocity at the collision position $x_{i,c}$ is

$$v_i^2 = \frac{2(E_i - k\,x_{i,c}^2)}{M_i} = \omega_i^2(A_i^2 - x_{i,c}^2) \tag{W17}$$

and is **uniquely** determined. This is independent of $\sigma$. Whether the wave is sinusoidal or rectangular, a mass passing through the same position with the same energy has the same velocity. $\checkmark$

**Stage 2: The elastic collision formula depends only on velocities.**

The collision formulas (W9)–(W10) of Theorem 11.1 are algebraic relations involving only $M_1, M_2, v_1, v_2$, containing no waveform information. Since Stage 1 showed that $v_i$ is $\sigma$-independent, $v_i'$ is also $\sigma$-independent. $\checkmark$

**Stage 3: Post-collision amplitudes and phases are also $\sigma$-independent.**

The post-collision total energy is

$$E_i' = k\,x_{i,c}^2 + \frac{1}{2}M_i\,v_i'^2 \tag{W18}$$

and the new amplitude is

$$A_i' = \sqrt{\frac{E_i'}{k}} = \sqrt{x_{i,c}^2 + \frac{v_i'^2}{\omega_i^2}} \tag{W19}$$

Since $x_{i,c}$, $v_i'$, and $\omega_i$ are all $\sigma$-independent, $A_i'$ is also $\sigma$-independent.

The initial phase $\phi_i' = \arctan(-v_i'/(\omega_i\,x_{i,c}))$ is similarly $\sigma$-independent. $\square$

**Corollary 11.1 (Exactness of the rectangular-wave limit).** The results of perfectly elastic collision in this model are **exactly identical** whether the oscillation is a smooth sinusoidal wave or a rectangular wave. The only information lost is the detail of the waveform (the shape of the transition region), while the physics of the collision — velocity redistribution, energy conservation, changes in amplitude and phase — does not depend on the waveform at all.

**Corollary 11.2.** This result reproduces, in a more general framework, the structural preservation in the rectangular-wave limit that was demonstrated for the sine-Gordon equation in §9. While the proof of structural preservation in sine-Gordon required assumptions of integrability and closed space, in this model the equivalence holds from **the single condition of energy conservation** alone (equation (W16)).

### 11.10 Verification Summary of This Model

| Property | Sinusoidal wave | Rectangular wave | $\sigma$-dependence | Remarks |
|:----:|:------:|:------:|:--------------:|:----:|
| Collision-time velocity $v_i$ | $\checkmark$ | $\checkmark$ | **Independent** | By (W17) |
| Post-collision velocity $v_i'$ | $\checkmark$ | $\checkmark$ | **Independent** | (W9)–(W10) |
| Post-collision amplitude $A_i'$ | $\checkmark$ | $\checkmark$ | **Independent** | (W19) |
| Post-collision phase $\phi_i'$ | $\checkmark$ | $\checkmark$ | **Independent** | (W14) |
| Energy conservation | $\checkmark$ | $\checkmark$ | **Independent** | (W15) |
| Momentum conservation | $\checkmark$ | $\checkmark$ | **Independent** | (W7) |
| Different masses | $\checkmark$ | $\checkmark$ | **Independent** | $M_1 \neq M_2$ |
| Arbitrary phases | $\checkmark$ | $\checkmark$ | **Independent** | $\phi_1, \phi_2$ arbitrary |

**Theorem 11.3 (Structural preservation of the perfectly elastic collision model).** The process of perfectly elastic collision at arbitrary phases $\phi_1, \phi_2$ of longitudinal waves in a perfect fluid with different specific gravities (masses $M_1 \neq M_2$) is completely described by wave equations (W3)–(W4) and elastic collision formulas (W9)–(W10), and all of the following properties are preserved in the rectangular-wave limit $\sigma \to 0$.

1. **Velocity redistribution** (elastic collision formula: $\sigma$-independent)
2. **Conservation of energy and momentum** ((W15) and (W7): $\sigma$-independent)
3. **Determination of amplitude and phase** ((W19) and (W14): $\sigma$-independent)

*Proof.* A direct consequence of Theorem 11.2. $\square$

The only assumptions of this entire formulation are the following two:

- **Momentum conservation**: $M_1 v_1 + M_2 v_2 = M_1 v_1' + M_2 v_2'$
- **Energy conservation**: $\frac{1}{2}M_1 v_1^2 + \frac{1}{2}M_2 v_2^2 = \frac{1}{2}M_1 v_1'^2 + \frac{1}{2}M_2 v_2'^2$

No additional assumptions of topological constraints, integrability, compactness of space, integer winding numbers, periodicity of potential, or identical masses are required, and the equivalence with the rectangular-wave limit follows automatically from energy conservation.

---

## References

[1] Frenkel, J. & Kontorova, T. "On the theory of plastic deformation and twinning." *Izv. Akad. Nauk SSSR, Ser. Fiz.* 1, 137–149 (1939).

[2] Scott, A.C., Chu, F.Y.F. & McLaughlin, D.W. "The Soliton: A New Concept in Applied Science." *Proc. IEEE* 61, 1443–1483 (1973).

[3] Rubinstein, J. "Sine-Gordon equation." *J. Math. Phys.* 11, 258–266 (1970).

[4] Rajaraman, R. *Solitons and Instantons: An Introduction to Solitons and Instantons in Quantum Field Theory.* North-Holland (1982).

[5] Bour, E. "Théorie de la déformation des surfaces." *J. Éc. Imp. Polytech.* 22, Cahier 39, 1–148 (1862).

[6] Enneper, A. "Über asymptotische Linien." *Nachr. Königl. Gesellsch. d. Wiss. Göttingen* 1870, 493–510 (1870).

[7] Bianchi, L. "Ricerche sulle superficie a curvatura costante e sulle elicoidi." *Ann. Sc. Norm. Super. Pisa* 2, 285–341 (1879).

[8] Skyrme, T.H.R. "A non-linear field theory." *Proc. R. Soc. Lond. A* 260, 127–138 (1961).

[9] Perring, J.K. & Skyrme, T.H.R. "A model unified field equation." *Nucl. Phys.* 31, 550–555 (1962).

[10] Manton, N.S. & Sutcliffe, P.M. *Topological Solitons.* Cambridge University Press (2004).

[11] Ablowitz, M.J., Kaup, D.J., Newell, A.C. & Segur, H. "Method for solving the sine-Gordon equation." *Phys. Rev. Lett.* 30, 1262–1264 (1973).

[12] Ablowitz, M.J., Kaup, D.J., Newell, A.C. & Segur, H. "The inverse scattering transform — Fourier analysis for nonlinear problems." *Stud. Appl. Math.* 53, 249–315 (1974).

[13] Faddeev, L.D. & Takhtajan, L.A. *Hamiltonian Methods in the Theory of Solitons.* Springer (1987).

[14] Zakharov, V.E. & Shabat, A.B. "Exact theory of two-dimensional self-focusing and one-dimensional self-modulation of waves in nonlinear media." *Sov. Phys. JETP* 34, 62–69 (1972).

[15] Rogers, C. & Schief, W.K. *Bäcklund and Darboux Transformations: Geometry and Modern Applications in Soliton Theory.* Cambridge University Press (2002).

[16] Gibbon, J.D., James, I.N. & Moroz, I.M. "The sine-Gordon equation as a model for a rapidly rotating baroclinic fluid." *Proc. R. Soc. Lond. A* 367, 219–237 (1979).

[17] Derrick, G.H. "Comments on nonlinear wave equations as models for elementary particles." *J. Math. Phys.* 5, 1252–1254 (1964).

[18] Manton, N.S. "Geometry of Skyrmions." *Commun. Math. Phys.* 111, 469–478 (1987).

[19] Barone, A. & Paternò, G. *Physics and Applications of the Josephson Effect.* Wiley (1982).

[20] Coleman, S. "Quantum sine-Gordon equation as the massive Thirring model." *Phys. Rev. D* 11, 2088–2097 (1975).

[21] Mandelstam, S. "Soliton operators for the quantized sine-Gordon equation." *Phys. Rev. D* 11, 3026–3030 (1975).

[22] Manton, N.S. "A remark on the scattering of BPS monopoles." *Phys. Lett. B* 110, 54–56 (1982).

[23] Ghoshal, S. & Zamolodchikov, A.B. "Boundary S matrix and boundary state in two-dimensional integrable quantum field theory." *Int. J. Mod. Phys. A* 9, 3841–3885 (1994).
