# The sine-Gordon Equation: Foundations of Topological Soliton Theory

## -- Classification of Topological Solutions and State Transformation Properties in Nonlinear Wave Equations --

**Author:** Noriaki Kihara (WF System Co., Ltd.)

**Date:** April 2026

**DOI:** 10.5281/zenodo.19650966

---

## 0. Purpose of This Paper

This paper provides an overview of the mathematical structure of the **sine-Gordon equation**, one of the nonlinear wave equations. The sine-Gordon equation is the most fundamental nonlinear equation that admits topological solitons as solutions, and it is widely applied in particle physics, condensed matter physics, and field theory.

In this paper, we organize the definition of the sine-Gordon equation, the classification of solutions (kinks, anti-kinks, breathers), phase shifts during collisions, and physical applications, tracing back to the original papers. Furthermore, we prove that soliton solutions are rigorously characterized by a finite number of parameters, and that all structural properties of solitons are preserved in the limit of the spread parameter $\sigma \to 0$ (rectangular wave limit). In addition, we define the transformation invariants of the sinusoidal-rectangular wave transformation (§11), explicitly state the limits of applicability of sine-Gordon theory (§9.8), and present an elastic state transformation model on phase space that preserves the transformation invariants as an alternative model that overcomes these limitations (§12).

---

## 1. Definition of the sine-Gordon Equation

### 1.1 The Equation

The 1+1-dimensional sine-Gordon equation is defined as follows [1, 2]:

$$\frac{\partial^2 \phi}{\partial t^2} - \frac{\partial^2 \phi}{\partial x^2} + \sin(\phi) = 0$$

Here $\phi = \phi(x, t)$ is a scalar field (phase field), and $x, t$ are dimensionless independent variables. In this paper, we treat the above equation as a dimensionless phase equation, and do not restore physical unit systems.

### 1.2 Origin of the Name

The name "sine-Gordon equation" derives from its similarity to the Klein-Gordon equation

$$\frac{\partial^2 \phi}{\partial t^2} - \frac{\partial^2 \phi}{\partial x^2} + m^2 \phi = 0$$

It is obtained by replacing the linear term $m^2 \phi$ with the nonlinear term $\sin(\phi)$, and the name is a portmanteau combining "sine" and "Gordon." This name was introduced by Rubinstein [3].

### 1.3 Lagrangian Density

The sine-Gordon equation is derived from the following Lagrangian density [4]:

$$\mathcal{L} = \frac{1}{2}\left(\frac{\partial \phi}{\partial t}\right)^2 - \frac{1}{2}\left(\frac{\partial \phi}{\partial x}\right)^2 - (1 - \cos\phi)$$

The potential $V(\phi) = 1 - \cos\phi$ attains its minimum value of 0 at $\phi = 2n\pi$ ($n \in \mathbb{Z}$). The existence of these infinitely many degenerate vacua is the origin of topological solitons.

### 1.4 Historical Background

The sine-Gordon equation was independently introduced in the following contexts:

- **Crystal dislocation theory** (1939): Frenkel and Kontorova introduced it as a model describing the motion of dislocations in crystal lattices [1].
- **Differential geometry** (19th century): In the theory of surfaces with constant negative curvature (pseudospheres), Bour [5] and Enneper [6] had derived an isomorphic equation. Bianchi [7] studied the Bäcklund transformation of this equation.
- **Particle physics** (1960s): Skyrme used a sine-Gordon type equation as a nonlinear field theory model for baryons [8]. Perring and Skyrme were the first to perform numerical calculations of kink-anti-kink collisions [9].

---

## 2. Classification of Solutions

### 2.1 Kink Solution (Winding Number $+1$)

The kink solution is a localized wave in which $\phi$ smoothly transitions from $0$ to $2\pi$ [2, 9]:

$$\phi_K(x, t) = 4\arctan\left(\exp\left(\frac{x - vt}{\sqrt{1 - v^2}}\right)\right)$$

Here $v$ ($|v| < 1$) is the velocity of the kink.

**Properties:**

- $x \to -\infty$: $\phi \to 0$; $x \to +\infty$: $\phi \to 2\pi$
- The phase increases by $2\pi$ -- **winding number (topological charge)** $N = +1$
- The denominator $\sqrt{1 - v^2}$ is the Lorentz factor, covariant under special relativity
- Rest energy $E_0 = 8m$ (in natural units $E_0 = 8$)

### 2.2 Anti-kink Solution (Winding Number $-1$)

The anti-kink solution is a localized wave in which $\phi$ transitions from $2\pi$ to $0$ [2]:

$$\phi_{\bar{K}}(x, t) = 4\arctan\left(\exp\left(-\frac{x - vt}{\sqrt{1 - v^2}}\right)\right)$$

**Properties:**

- $x \to -\infty$: $\phi \to 2\pi$; $x \to +\infty$: $\phi \to 0$
- The phase decreases by $2\pi$ -- **winding number** $N = -1$
- Has the same mass and velocity characteristics as the kink

### 2.3 Definition and Conservation of Winding Number

**Definition 2.1.** The **winding number** (topological charge) of a field configuration $\phi(x, t)$ is defined as follows [4, 10]:

$$N = \frac{1}{2\pi}\int_{-\infty}^{+\infty} \frac{\partial \phi}{\partial x}\,dx = \frac{\phi(+\infty) - \phi(-\infty)}{2\pi}$$

**Proposition 2.1 (Conservation of winding number).** Under the time evolution of the sine-Gordon equation, the winding number $N$ is conserved.

*Proof.* The boundary values $\phi(\pm\infty)$ are fixed at vacuum values ($2n\pi$), and cannot change by integer multiples of $2\pi$ through continuous time evolution. Therefore $N$ is conserved as a time-independent integer value. $\square$

This conservation law guarantees the "particle-like nature" of solitons -- creation and annihilation are possible only in pairs (kink + anti-kink).

### 2.4 Breather Solution (Bound Oscillating State)

A breather is an oscillating solution in which a kink and an anti-kink are bound, with winding number $N = 0$ [2, 11]:

$$\phi_B(x, t) = 4\arctan\left(\frac{\sqrt{1 - \omega^2}\;\cos(\omega t)}{\omega\;\cosh\left(\sqrt{1 - \omega^2}\;x\right)}\right)$$

Here $\omega$ ($0 < \omega < 1$) is the internal oscillation frequency of the breather.

**Properties:**

- A localized oscillating solution that decays exponentially in space
- The smaller the internal frequency $\omega$, the larger the spatial extent of the breather
- The rest mass of the breather is $M_B = 16\sqrt{1 - \omega^2}$
- In the limit $\omega \to 0$, the breather mass approaches that of a kink-anti-kink pair ($2 \times 8 = 16$)
- In the limit $\omega \to 1$, the breather mass approaches $0$, transitioning to linear small oscillations (mesons)

### 2.5 Summary of Solution Classification

| Solution Type | Winding Number $N$ | Localization | Oscillation | Mass |
|:--------:|:----------:|:------:|:----:|:----:|
| Kink | $+1$ | Localized | None | $8$ |
| Anti-kink | $-1$ | Localized | None | $8$ |
| Breather | $0$ | Localized | Yes | $16\sqrt{1-\omega^2}$ |
| Meson (linear wave) | $0$ | Non-localized | Yes | $\geq 2$ (threshold) |

---

## 3. Soliton Collisions and Phase Shifts

### 3.1 Kink-Kink Collision

The collision between kinks of the same sign ($N = +1$ and $N = +1$) is repulsive, and the two-body solution is given by [2, 12]:

$$\phi_{KK}(x, t) = 4\arctan\left(\frac{v\,\sinh\left(\frac{x}{\sqrt{1-v^2}}\right)}{\cosh\left(\frac{vt}{\sqrt{1-v^2}}\right)}\right)$$

After the collision, the two kinks separate again, completely preserving their shapes. However, compared to the case without collision, the position of each kink is shifted by a **phase shift** $\Delta x$:

$$\Delta x = \frac{\sqrt{1 - v^2}}{v}\ln(v^2)$$

### 3.2 Kink-Anti-kink Collision

The collision between a kink and an anti-kink of opposite signs ($N = +1$ and $N = -1$) is attractive, and the two-body solution is given by [2, 9, 12]:

$$\phi_{K\bar{K}}(x, t) = 4\arctan\left(\frac{v^{-1}\,\cosh\left(\frac{vt}{\sqrt{1-v^2}}\right)}{\sinh\left(\frac{x}{\sqrt{1-v^2}}\right)}\right)$$

**Properties:**

- With sufficient energy, the kink and anti-kink pass through each other and separate again (perfectly elastic collision: velocity and shape are identical before and after, only a phase shift occurs)
- With insufficient energy, they form a breather (bound state) (oscillating permanently without decay)

**Remark:** In the sine-Gordon equation as a classical completely integrable system, the annihilation of a kink-anti-kink pair followed by conversion to radiation waves (mesons) does not occur. This is a phenomenon that arises only through quantum effects or in systems where integrability is broken (perturbation, discretization, higher-dimensionalization).

### 3.3 Physical Meaning of the Phase Shift

**Proposition 3.1 (Shape preservation and phase shift).** In soliton collisions of the sine-Gordon equation:

1. The **shape** of each soliton after the collision is identical to that before the collision.
2. The **position** of each soliton after the collision is displaced by a phase shift $\Delta x$ compared to the case without collision.
3. The total winding number is conserved before and after the collision.

*Proof.* The sine-Gordon equation can be solved exactly by the inverse scattering method [12], and multi-soliton solutions asymptotically separate into superpositions of individual solitons. The asymptotic shape of each soliton is identical to that before the collision, but with an added phase shift. Since the winding number is a topological invariant, it is conserved. $\square$

This property is critically important for the present theory. The **type** of soliton (= particle) (winding number = orientational structure) is invariant under interactions, and only the **position** (phase) changes.

---

## 4. Complete Integrability and the Inverse Scattering Method

### 4.1 Lax Pair

The sine-Gordon equation is a **completely integrable system** and possesses a Lax pair [12, 13].

Using light-cone coordinates $u = \frac{1}{2}(x - t)$, $v = \frac{1}{2}(x + t)$, the equation becomes

$$\frac{\partial^2 \phi}{\partial u\,\partial v} = \sin\phi$$

This equation is expressed as the compatibility condition of the following Lax pair:

$$\Psi_u = U\Psi, \quad \Psi_v = V\Psi$$

Here $U$ and $V$ are $2 \times 2$ matrices depending on a spectral parameter $\lambda$ [12].

### 4.2 Infinitely Many Conserved Quantities

As a consequence of complete integrability, the sine-Gordon equation possesses infinitely many independent conserved quantities [12, 14]. The first three are as follows:

- **Energy:**
$$E = \int_{-\infty}^{+\infty}\left[\frac{1}{2}\left(\frac{\partial\phi}{\partial t}\right)^2 + \frac{1}{2}\left(\frac{\partial\phi}{\partial x}\right)^2 + (1-\cos\phi)\right]dx$$

- **Momentum:**
$$P = -\int_{-\infty}^{+\infty}\frac{\partial\phi}{\partial t}\frac{\partial\phi}{\partial x}\,dx$$

- **Winding number:**
$$N = \frac{1}{2\pi}\int_{-\infty}^{+\infty}\frac{\partial\phi}{\partial x}\,dx$$

### 4.3 Bäcklund Transformation

The sine-Gordon equation possesses a Bäcklund transformation, which can generate new solutions from known ones [7, 15].

Given a known solution $\phi_0$, the solution $\phi_1$ of the following system of equations is also a solution of the sine-Gordon equation:

$$\frac{\partial}{\partial u}\left(\frac{\phi_1 + \phi_0}{2}\right) = a\,\sin\left(\frac{\phi_1 - \phi_0}{2}\right)$$

$$\frac{\partial}{\partial v}\left(\frac{\phi_1 - \phi_0}{2}\right) = \frac{1}{a}\sin\left(\frac{\phi_1 + \phi_0}{2}\right)$$

Here $a$ is the parameter of the Bäcklund transformation. Starting from the vacuum solution $\phi_0 = 0$, one Bäcklund transformation yields the kink solution, and two transformations yield the kink-kink solution or breather solution [15].

---

## 5. Extension to Higher Dimensions

### 5.1 2+1 and 3+1 Dimensions

The direct extension of the sine-Gordon equation to higher dimensions is

$$\frac{\partial^2\phi}{\partial t^2} - \nabla^2\phi + \sin(\phi) = 0$$

However, in 2+1 dimensions and above, this equation is no longer completely integrable [16]. Furthermore, by Derrick's theorem [17], stable localized solutions (solitons) consisting only of scalar fields do not exist in 3+1 dimensions.

### 5.2 Methods of Stabilization

The following methods are known for obtaining stable solitons in higher dimensions [10, 18]:

1. **Skyrme model**: Adding higher-order derivative terms (Skyrme term) [8, 18].
2. **Coupling with gauge fields**: Coupling a gauge field to the scalar field.
3. **Topological stabilization**: When the field carries a non-trivial topological charge, energy minimization prohibits the collapse of the soliton [10].

### 5.3 Stabilization on Closed Spaces

When dealing with solitons on higher-dimensional spaces, a mechanism to circumvent Derrick's theorem is needed. On a closed curve $S^1$ of circumference $L$, the premises of Derrick's theorem (asymptotic conditions on infinite space) break down, so stable soliton solutions can exist [10].

---

## 6. Physical Applications

### 6.1 Condensed Matter Physics

- **Josephson junctions**: The motion of magnetic flux quanta in tunnel junctions between superconductors is described by the sine-Gordon equation [19]. Kinks correspond to magnetic flux quanta (fluxons).
- **Crystal dislocations**: Dislocations in crystal lattices in the Frenkel-Kontorova model correspond to kink solutions [1].

### 6.2 Particle Physics

- **Skyrme model**: Baryons are described as topological solitons (skyrmions) [8, 18]. The baryon number corresponds to the topological charge.
- **Equivalence between the quantum sine-Gordon model and the massive Thirring model**: Coleman [20] and Mandelstam [21] showed that the quantum theory of the sine-Gordon model is equivalent to the massive Thirring model (a fermionic field theory). The **boson-fermion correspondence**, in which bosonic solitons correspond to fermions, is one of the most remarkable results of 1+1-dimensional field theory.

### 6.3 Differential Geometry

- **Pseudospheres**: The compatibility conditions between the first fundamental form and the second fundamental form of surfaces with constant negative curvature ($K = -1$) reduce to the sine-Gordon equation [5, 6]. Kink solutions correspond to geodesics on pseudospheres.

---

## 7. Implications for Physical Interpretation

### 7.1 Summary of Particle-like Properties of Solitons

The soliton solutions of the sine-Gordon equation naturally possess the following particle-like properties:

| Soliton Property | Correspondence to Particles |
|:-----------:|:----------------:|
| Winding number $N = \pm 1$ | Particle / anti-particle |
| Conservation of winding number | Conservation of particle number |
| Shape preservation during collisions | Identity of particles |
| Phase shift (position change after collision) | Scattering through interactions |
| Lorentz covariance | Special relativistic structure |
| Breather ($N = 0$ bound state) | Meson-like bound state |
| Periodicity of the potential ($2\pi$) | Discreteness of phase |
| Boson-fermion correspondence [20, 21] | Spin-statistics |

### 7.2 Significance of the Rectangular Wave Description

As shown in §8--§9, all of these particle-like properties are preserved in the rectangular wave limit $\sigma \to 0$. This means that the essential information of solitons is not contained in the smooth shape of the continuous field, but is concentrated in the topological and scattering discrete structure.

---

## 8. Finite Parameter Representation of Solitons

### 8.1 Motivation: From Continuous Fields to Finite Parameters

Solutions of the sine-Gordon equation are defined as continuous fields $\phi(x, t)$, but due to their localization, soliton solutions are **completely characterized by a finite number of parameters**. This is not coincidental. Solitons are **finite-mode solutions** in the mode decomposition of nonlinear fields, solutions restricted to finite-dimensional invariant manifolds (moduli spaces) within the infinite-dimensional field theory [10, 22].

In this section, we provide a formulation that describes the soliton solutions of the sine-Gordon equation with a finite number of parameters, and explicitly show how solutions of the continuous wave equation are reduced to a small number of discrete pieces of information.

### 8.2 Five-Parameter Representation of the Kink

We reproduce the kink solution from §2.1:

$$\phi_K(x, t) = 4\arctan\left(\exp\left(\frac{x - x_0 - vt}{\ell}\right)\right)$$

Here $\ell = \sqrt{1 - v^2}$ is the width parameter of the kink. This solution is **completely** determined by the following five parameters.

**Definition 8.1 (Soliton parameters).** The five parameters characterizing a 1+1-dimensional sine-Gordon kink are defined as follows:

| Symbol | Name | Definition | Range |
|:----:|:----:|:----:|:----:|
| $N$ | Winding number | Topological charge $\pm 1$ | $\{-1, +1\}$ |
| $A$ | Maximum wave height | Total change in $\phi$ $= 2\pi |N|$ | $2\pi$ (fixed) |
| $\theta_0$ | Center phase | Spatial phase of soliton center $= 2\pi x_0 / L$ | $[0, 2\pi)$ |
| $\alpha$ | Velocity phase | Phase angle representation of velocity $v = \sin\alpha$ | $(-\pi/2, \pi/2)$ |
| $\sigma$ | Spread parameter | Phase angle representation of width $\ell = \cos\alpha$ | $(0, 1]$ |

**Proposition 8.1 (Completeness of parameters).** A single kink (anti-kink) solution of the sine-Gordon equation is uniquely determined by the above five parameters $(N, A, \theta_0, \alpha, \sigma)$. Conversely, the above five parameters can be uniquely read from any kink solution.

*Proof.* The general form of the kink solution is $\phi = 4\arctan(\exp(\pm(x - x_0 - vt)/\sqrt{1-v^2}))$, where the sign $\pm$ determines $N$, the displacement $x_0$ determines $\theta_0$, and the velocity $v$ determines $\alpha$. $A = 2\pi$ is fixed from the winding number $|N| = 1$, and $\sigma = \cos\alpha$ is determined subordinately from $\alpha$. The free parameters are the three quantities $(N, \theta_0, \alpha)$, and $A$ and $\sigma$ are derived from these. $\square$

**Remark.** The free parameters are essentially three ($N, \theta_0, \alpha$), but the five-parameter representation is adopted for clarity of physical interpretation. $A$ is a dependent parameter determined from the winding number, and $\sigma$ from the velocity phase.

### 8.3 Relationship Between Velocity Phase and Spread Parameter

The following identity holds between the velocity phase $\alpha$ and the spread parameter $\sigma$:

$$\sin^2\alpha + \sigma^2 = v^2 + \ell^2 = 1$$

That is, $(\sin\alpha, \sigma) = (v, \ell)$ is a point on the unit circle.

**Physical interpretation:**

- At $\alpha = 0$ (rest), $\sigma = 1$: the soliton has maximum width
- As $\alpha \to \pm\pi/2$ (speed of light), $\sigma \to 0$: the width of the soliton contracts to zero (Lorentz contraction)
- The faster the motion, the narrower the spatial extent -- **motion and shape are unified by a single phase angle $\alpha$**

### 8.4 Six-Parameter Representation of the Breather

The breather solution of §2.4 requires the internal oscillation frequency $\omega$ in addition to the five kink parameters.

**Definition 8.2 (Breather parameters).**

| Symbol | Name | Definition | Range |
|:----:|:----:|:----:|:----:|
| $N$ | Winding number | $0$ (kink-anti-kink bound state) | $\{0\}$ (fixed) |
| $A$ | Maximum wave height | $4\arctan(\sqrt{1-\omega^2}/\omega)$ | $(0, \pi)$ |
| $\theta_0$ | Center phase | Spatial phase of breather center | $[0, 2\pi)$ |
| $\alpha$ | Velocity phase | Overall motion $v = \sin\alpha$ | $(-\pi/2, \pi/2)$ |
| $\sigma$ | Spread parameter | $1/\sqrt{1-\omega^2}$ | $[1, \infty)$ |
| $\omega$ | Internal phase velocity | Angular frequency of internal oscillation | $(0, 1)$ |

The breather has $N = 0$ but possesses the additional degree of freedom of internal oscillation $\omega$.

### 8.5 Parameter Representation of Multi-Soliton States

A multi-body state consisting of $n$ solitons is represented as a set of parameters for each soliton [10, 12].

**Definition 8.3 ($n$-soliton state).** A state consisting of $n$ solitons is defined by the following parameter set:

$$\mathcal{S}_n = \{(N_i, A_i, \theta_{0,i}, \alpha_i, \sigma_i)\}_{i=1}^{n}$$

**Proposition 8.2 (Asymptotic completeness).** A field configuration $\phi(x, t)$ consisting of $n$ sufficiently separated solitons is asymptotically completely determined by the parameter set $\mathcal{S}_n$. In the region where overlap between solitons is negligible,

$$\phi(x, t) \approx \sum_{i=1}^{n} \phi_i(x, t; N_i, \theta_{0,i}, \alpha_i)$$

holds, where $\phi_i$ is the single-soliton solution for each soliton.

*Proof.* By the inverse scattering method [12], the $n$-soliton solution is completely determined by $n$ scattering data (eigenvalues and norming constants). In the asymptotic region, the scattering data and the above parameters are in one-to-one correspondence. $\square$

### 8.6 Description of Collisions in Parameter Space

The soliton collisions of §3 can be described as mappings on parameter space.

**Proposition 8.3 (Collision mapping).** The collision of two solitons $\mathcal{S}_1 = (N_1, \theta_{0,1}, \alpha_1)$ and $\mathcal{S}_2 = (N_2, \theta_{0,2}, \alpha_2)$ is described as the following parameter transformation:

$$\mathcal{S}_1' = (N_1,\; \theta_{0,1} + \Delta\theta_1,\; \alpha_1)$$
$$\mathcal{S}_2' = (N_2,\; \theta_{0,2} + \Delta\theta_2,\; \alpha_2)$$

That is:

1. **Winding numbers $N_i$ are invariant** -- the type of particle does not change
2. **Velocity phases $\alpha_i$ are invariant** -- velocity (momentum) is conserved
3. **Only the center phases $\theta_{0,i}$ change** -- they are displaced by the phase shift $\Delta\theta_i$

**Corollary 8.1.** All information about a soliton collision is concentrated in the phase shifts $\Delta\theta_i$. The phase shift is computable from the incident parameters [12] and takes the form

$$\Delta\theta_1 = -\Delta\theta_2 = f(N_1, N_2, \alpha_1 - \alpha_2)$$

That is, the collision depends only on the relative velocity phase $\alpha_1 - \alpha_2$.

### 8.7 Reduction from Continuous Fields to Parameter Representation

Summarizing the above, the solution space of the sine-Gordon equation is stratified as follows:

| Description Level | Degrees of Freedom | Information Content |
|:----------:|:------:|:------:|
| Continuous field $\phi(x, t)$ | Infinite-dimensional | Real value at each spacetime point |
| $n$-soliton solution | $3n$-dimensional | Each soliton $(N_i, \theta_{0,i}, \alpha_i)$ |
| Change after collision | $n$-dimensional | Only phase shifts $\Delta\theta_i$ |

**Proposition 8.4 (Finite-dimensionality of the soliton sector).** Within the solution space of the sine-Gordon equation, the sector containing $n$ solitons is described as a $3n$-dimensional finite-dimensional manifold (moduli space) $\mathcal{M}_n$ [10, 22]. The infinite degrees of freedom of the continuous field are **rigorously** reduced to a finite number of parameters in the soliton sector.

*Proof.* In the inverse scattering method, the $n$-soliton sector is characterized by $n$ discrete eigenvalues $\zeta_i$ and their associated norming constants $c_i$ [12, 13]. Each eigenvalue determines the winding number $N_i$ and velocity $v_i$ (= $\sin\alpha_i$), and the norming constants determine the initial positions $x_{0,i}$ (= $\theta_{0,i}$). Excluding the radiation component (continuous spectrum), the degrees of freedom of the soliton sector are exactly $3n$. $\square$

This reduction is not an approximation. In the soliton sector, the infinite-dimensional partial differential equation **rigorously** reduces to a finite-dimensional system of ordinary differential equations.

### 8.8 Relationship with Rectangular Waves

Let us re-examine the shape of the kink solution:

$$\phi_K(x) = 4\arctan(e^{x/\ell}) = 2\pi \cdot \frac{1}{1 + e^{-x/\ell}}$$

The right-hand side is a logistic function (sigmoid function), which converges to a **step function** (half-period of a rectangular wave) in the limit $\ell \to 0$:

$$\lim_{\ell \to 0} \phi_K(x) = \begin{cases} 0 & (x < 0) \\ 2\pi & (x > 0) \end{cases}$$

**Proposition 8.5 (Rectangular wave limit).** In the limit of the spread parameter $\sigma = \ell \to 0$ (i.e., $\alpha \to \pm\pi/2$, the speed-of-light limit), the kink solution converges to a rectangular wave (step function). Conversely, a kink with finite width $\sigma > 0$ can be understood as a **smoothing** of the rectangular wave.

*Proof.* $\phi_K(x) = 2\pi/(1 + e^{-x/\ell})$ is a logistic function, and converges pointwise to the Heaviside step function $2\pi \cdot H(x)$ as $\ell \to 0$. $\square$

**Corollary 8.2.** A soliton can be interpreted as a finite-mode approximation of a rectangular wave. The continuous $\tanh$-type shape is merely a smoothed version of a discrete step function. The spread parameter $\sigma$ controls the degree of smoothing.

From this viewpoint, the parameter representation of solitons is reinterpreted as follows:

| Parameter | Meaning in Continuous Field | Meaning in Rectangular Wave (Discrete) |
|:----------:|:------------:|:--------------------:|
| $N = \pm 1$ | Winding number | Direction of step (ascending/descending) |
| $A = 2\pi$ | Total phase change | Height of step |
| $\theta_0$ | Position of soliton center | Edge position of step |
| $\alpha$ | Phase angle of velocity | Direction of edge movement |
| $\sigma$ | Width of smoothing | **Equals rectangular wave at $\sigma \to 0$** |

---

## 9. Preservation of Properties in the Rectangular Wave Limit

### 9.1 Purpose

The properties of the sine-Gordon equation described in §1--§7 have been proved for the continuous field $\phi(x, t)$. As shown in §8.8, the soliton solution converges to a rectangular wave (step function) in the limit of the spread parameter $\sigma \to 0$.

In this section, for each property from §1--§7, we systematically verify that (a) it holds for soliton solutions (which is trivial since they are exact solutions), and (b) it is preserved in the rectangular wave limit $\sigma \to 0$.

If all properties are preserved in the rectangular wave limit, then the rectangular wave description is not a "convenient approximation" but rather a **rigorous limit that completely inherits all structural properties** of the soliton.

### 9.2 Framework for Verification

Each property is verified in the following three stages:

1. **Validity for soliton solutions**: Does it satisfy the sine-Gordon equation as an exact solution?
2. **$\sigma$-dependence analysis**: Does the property depend on the spread parameter $\sigma$?
3. **Existence of the rectangular wave limit**: Is the property preserved or does it diverge as $\sigma \to 0$?

### 9.3 Definition and Conservation of Winding Number (§2.3 → Rectangular Wave)

**Property from §2.3:** The winding number $N = [\phi(+\infty) - \phi(-\infty)] / 2\pi$ is conserved under time evolution.

**Validity for soliton solutions:** The kink solution gives $\phi(-\infty) = 0$, $\phi(+\infty) = 2\pi$, so $N = +1$. Since this is an exact solution of the sine-Gordon equation, the conservation of $N$ is trivial. $\checkmark$

**$\sigma$-dependence:** The definition of winding number $N = [\phi(+\infty) - \phi(-\infty)] / 2\pi$ depends **only on the boundary values** and does not depend at all on the internal shape of $\phi$ (hence on $\sigma$).

**Proposition 9.1 ($\sigma$-independence of the winding number).** The winding number $N$ does not depend on the spread parameter $\sigma$. Even in the rectangular wave limit $\sigma \to 0$, $N$ retains the same integer value.

*Proof.* For any $\sigma > 0$, $\phi_K^\sigma(\pm\infty)$ are fixed at the vacuum values $0, 2\pi$. In the limit $\sigma \to 0$, $\phi_K^\sigma(x) \to 2\pi \cdot H(x)$ converges pointwise, but since $H(-\infty) = 0$, $H(+\infty) = 1$, we have $N = (2\pi \cdot 1 - 2\pi \cdot 0)/2\pi = 1$, which is invariant. $\square$

**Winding number for rectangular waves:** For the rectangular wave $\phi(x) = 2\pi H(x - x_0)$, $N = +1$. If the step is ascending, $N = +1$; if descending, $N = -1$. The winding number is precisely the sign of the step itself. $\checkmark$

### 9.4 Conservation Law for Winding Number (Proposition 2.1 → Rectangular Wave)

**Property from §2.3:** $N$ is conserved under time evolution. Creation and annihilation are possible only in kink-anti-kink pairs.

**$\sigma$-dependence:** The proof of the conservation law (§2.3) relies only on the fact that "the boundary values are fixed at vacuum values and cannot change by integer multiples of $2\pi$ through continuous time evolution." This argument does not depend on $\sigma$.

**Conservation for rectangular waves:** The boundary values of the rectangular wave are $0$ and $2\pi$, taking only discrete values. Discrete values cannot change under continuous deformation. Therefore, the conservation of winding number for rectangular waves is **rather more obvious** than for the continuous field case.

**Proposition 9.2.** The conservation law for the winding number holds rigorously in the rectangular wave limit. Furthermore, since $\phi$ itself takes only discrete values ($0$ or $2\pi$) in the rectangular wave case, conservation is guaranteed not only topologically but also by the **discreteness of values**. $\checkmark$

### 9.5 Existence of Kink and Anti-kink Solutions (§2.1, §2.2 → Rectangular Wave)

**Property from §2.1--2.2:** Localized solutions with winding number $\pm 1$ exist and propagate with velocity $v$.

**Validity for soliton solutions:** The kink $\phi_K = 4\arctan(\exp((x-vt)/\sigma))$ is an exact solution. $\checkmark$

**Rectangular wave limit:** In the limit $\sigma \to 0$, $\phi_K \to 2\pi \cdot H(x - vt)$. This is a step function moving with velocity $v$, preserving localization and propagation properties.

**Proposition 9.3 (Existence of rectangular wave kink).** The rectangular wave kink $\phi(x, t) = 2\pi \cdot H(x - vt)$ has the following properties:
1. Winding number $N = +1$
2. Moves with velocity $v$
3. Localized at a single point in space ($x = vt$)
4. Holds as a **distributional solution** of the sine-Gordon equation

*Proof.* (1)--(3) are obvious from the definition. For (4): the distributional derivatives of $\phi = 2\pi H(x - vt)$ are $\phi_x = 2\pi\delta(x-vt)$, $\phi_{xx} = 2\pi\delta'(x-vt)$, $\phi_{tt} = 2\pi v^2\delta'(x-vt)$, and $\sin(\phi) = \sin(2\pi H(x-vt)) = 0$ (for $x \neq vt$, $\phi = 0$ or $2\pi$, both giving $\sin = 0$). The wave equation part $\phi_{tt} - \phi_{xx} = 2\pi(v^2-1)\delta'(x-vt)$ vanishes only when $v = \pm 1$ (speed of light). For $v \neq 1$, the rectangular wave is not a classical solution, and smoothing with finite width $\sigma > 0$ is necessary. $\square$

**Remark.** The fact that the rectangular wave is not an exact solution for $v \neq 1$ is consistent with the relation $\sigma = \cos\alpha = \sqrt{1-v^2}$ from §8.3. If $v < 1$, then $\sigma > 0$, and the soliton has finite width. The rectangular wave ($\sigma = 0$) corresponds only to $v = 1$ (speed of light). The fact that particles with velocities below the speed of light have a finite "spread" is qualitatively consistent with the spread of wave packets in quantum mechanics.

### 9.6 Breather Solution (§2.4 → Rectangular Wave)

**Property from §2.4:** A breather exists as a kink-anti-kink bound state.

**Rectangular wave limit:** The breather is an oscillating solution with winding number $N = 0$. In the rectangular wave limit, this corresponds to a **rectangular pulse** oscillating between $0$ and $2\pi$.

**Proposition 9.4.** The rectangular wave limit of the breather is the oscillation of a spatially localized rectangular pulse. It periodically transitions between two states: presence ($\phi = 2\pi$) and absence ($\phi = 0$) of the pulse. The mass $M_B = 16\sqrt{1-\omega^2}$ depends on $\omega$ and does not directly depend on $\sigma$. $\checkmark$

### 9.7 Soliton Collisions and Phase Shift (§3 → Rectangular Wave)

**Property from §3:** After collision: (i) shape preservation, (ii) phase shift, (iii) winding number conservation.

This is the most important property, and each item is verified individually.

**(i) Shape preservation:**

**Soliton solution:** The soliton after collision has the same shape (same $\sigma$) as before the collision [12]. $\checkmark$

**Rectangular wave limit:** The rectangular wave has $\sigma = 0$, and the only shape parameter is "step height $2\pi$." If $\sigma = 0$ is preserved before and after the collision, shape preservation is automatically satisfied. Since $\sigma$ is a conserved quantity in soliton collisions (Proposition 8.3: $\alpha_i$ invariant $\Rightarrow$ $\sigma_i = \cos\alpha_i$ invariant), $\sigma = 0$ is also conserved. $\checkmark$

**(ii) Phase shift:**

**Soliton solution:** The phase shift after collision is $\Delta\theta = f(N_1, N_2, \alpha_1 - \alpha_2)$ [12]. $\checkmark$

**Rectangular wave limit:** The phase shift formula depends only on $\alpha_i$ (velocity phase) and does not depend on $\sigma_i$ (spread) (see §8.6, Corollary 8.1). Therefore, even at $\sigma \to 0$, the value of the phase shift does not change.

**Proposition 9.5 ($\sigma$-independence of phase shift).** The phase shift $\Delta\theta$ of a soliton collision does not depend on the spread parameter $\sigma_i$ of each soliton. Even in the rectangular wave limit $\sigma \to 0$, the phase shift takes the same value as in the soliton case.

*Proof.* In the inverse scattering method, the phase shift is computed only from the discrete eigenvalues $\zeta_i$ of the scattering data [12, 13]. The $\zeta_i$ are determined by the soliton velocity $v_i = \sin\alpha_i$ and winding number $N_i$, and are independent of the norming constants (position information) and radiation component. The soliton width $\sigma_i$ is related to the norming constants but does not contribute to the phase shift. $\square$

**(iii) Winding number conservation:** Already shown in §9.4. $\checkmark$

### 9.8 Complete Integrability and Conserved Quantities (§4 → Rectangular Wave)

**Property from §4:** Conservation of energy, momentum, and winding number.

**Energy:**

The energy of a soliton is $E = 8/\sigma = 8/\cos\alpha = 8/\sqrt{1-v^2}$. As $\sigma \to 0$, the energy diverges.

**Proposition 9.6 ($\sigma$-dependence of energy).** The energy of a soliton depends on $\sigma$, and diverges as $E \to \infty$ when $\sigma \to 0$. This means that the rectangular wave limit has infinite energy.

**Proposition 9.6a (Limits of applicability of this theory).** The sine-Gordon soliton theory is based on the following implicit assumptions:

1. **Fixed origin**: The initial value of the field is fixed at the vacuum value $\phi(-\infty) = 0$ (§2.1).
2. **Periodicity of the potential**: $V(\phi) = 1-\cos\phi$ has period $2\pi$, and the target space is effectively $S^1$. The size of one "circuit" of the potential ($2\pi$) is fixed in the theory and has no freedom for modification.
3. **Forcing the terminus to vacuum values**: The requirement of finite energy on the infinite space $(-\infty, +\infty)$ forces $\phi(+\infty) = 2n\pi$ ($n \in \mathbb{Z}$). This is because $V(2n\pi) = 0$ (potential minimum), while $\phi = (2n+1)\pi$ gives $V = 2$ (potential maximum), so energy diverges on infinite space, and such configurations are excluded.
4. **Exclusion of fractional winding numbers**: As a consequence of the above, the winding number $N$ is always an integer. Half-integer winding numbers ($\phi(+\infty) = (2n+1)\pi$) are automatically excluded by the finite energy requirement.
5. **Absence of boundary reflection**: Since the space is infinite or closed ($S^1$), there are no physical boundaries, and boundary reflection of solitons does not occur.
6. **Perfectly elastic collisions and absence of interactions**: This theory is a completely integrable system (§4) possessing infinitely many conserved quantities. As a consequence of this excessive constraint, soliton collisions are described as perfectly elastic collisions, where the velocity, shape, and energy of each soliton are **individually** conserved before and after the collision (§3.3). To avoid violating momentum conservation, only the position of each soliton (phase shift $\Delta x$) changes, and exchange of velocities, deformation of shapes, and dissipation of energy never occur. This effectively means that solitons behave as **free particles that do not interact** with each other. For identical particles (e.g., kink-kink collision), "whether they passed through or repelled" is in principle indistinguishable.

*Consequence.* This theory is consistently applicable only in the following two cases:

- **Infinite space** $(-\infty, +\infty)$: The endpoints are at infinity and are never reached. The energy is $E = 8$ per kink (at rest) and finite.
- **Closed curve** $S^1$ (circumference $L$): Periodic boundary conditions $\phi(x+L) = \phi(x) + 2\pi N$. No boundaries exist.

An **open string** (Dirichlet / Neumann boundary conditions on a finite interval $[0, L]$) lies outside the scope of this theory. The reason is that a finite interval $[0, L]$ is a contractible space ($\pi_1([0,L]) = 0$), so the winding number is not defined as a topological invariant, and the topological stability of solitons is not guaranteed. Furthermore, on a finite interval, $\phi(L) = (2n+1)\pi$ (half-integer winding number, potential maximum) is energetically allowed, so the restriction to integer winding numbers does not hold.

**Self-cancellation structure of energy:** The energy of each $2\pi$ rotation (one kink) is self-contained at $E = 8$ (at rest). This is because the potential $V(\phi) = 1-\cos\phi$ returns from $0 \to 0$ over one period. A kink-anti-kink pair ($4\pi$ rotation, winding number $N = 0$) can return its energy to zero through pair annihilation. The only observable permanent energy is the contribution $8|N|$ from the remaining kinks that cannot undergo pair annihilation.

**Momentum and winding number:** The momentum $P = -Ev$ has a similar structure. The winding number is $\sigma$-independent as shown in §9.3. $\checkmark$

### 9.9 Extension to Higher Dimensions (§5 → Rectangular Wave)

**Property from §5:** By Derrick's theorem, stable localized solutions do not exist in infinite space of 3+1 dimensions.

**Effect on rectangular wave limit:** The premise of Derrick's theorem is the **existence of finite-energy solutions in infinite space** [17]. On a closed curve $S^1$ (circumference $L$), the premises of the theorem break down for the following two reasons:

1. **Space is closed**: The range of $\lambda$ in the scale transformation $\phi(x) \to \phi(\lambda x)$ is restricted
2. **Topological protection**: The winding number $N \in \pi_1(S^1) = \mathbb{Z}$ is protected as a topological invariant

**Proposition 9.7.** On a closed curve $S^1$ (circumference $L$), localized solutions of any width, including rectangular waves ($\sigma = 0$), are topologically stable. Derrick's theorem does not apply to closed spaces. $\checkmark$

### 9.10 Summary of Verification Results

| Property | Section | Soliton | $\sigma$-dependence | Rectangular Wave Limit | Remarks |
|:----:|:--:|:-------:|:--------------:|:----------:|:----:|
| Definition of winding number | §2.3 | $\checkmark$ | **Independent** | $\checkmark$ | Boundary values only |
| Conservation of winding number | §2.3 | $\checkmark$ | **Independent** | $\checkmark$ | Topological |
| Existence of kink/anti-kink | §2.1--2 | $\checkmark$ | Shape only | $\checkmark$ | Distributional solution |
| Existence of breather | §2.4 | $\checkmark$ | Shape only | $\checkmark$ | Rectangular pulse |
| Shape preservation after collision | §3.1--2 | $\checkmark$ | **Independent** | $\checkmark$ | $\sigma$ conserved |
| Phase shift after collision | §3.3 | $\checkmark$ | **Independent** | $\checkmark$ | $\alpha$ only |
| Winding number conservation (collision) | §3.3 | $\checkmark$ | **Independent** | $\checkmark$ | Topological |
| Energy conservation | §4.2 | $\checkmark$ | Dependent | $\checkmark$* | *Self-contained per $2\pi$ rotation |
| Momentum conservation | §4.2 | $\checkmark$ | Dependent | $\checkmark$* | *Same as above |
| Derrick avoidance | §5.1--2 | -- | -- | $\checkmark$ | Closed space $S^1$ |
| Parameter representation $(N, \theta_0, \alpha)$ | §8 | $\checkmark$ | **Independent** | $\checkmark$ | $\sigma$ is subordinate |

**Theorem 9.1 (Structural preservation of the rectangular wave limit).** The following structural properties of soliton solutions of the sine-Gordon equation are all preserved in the rectangular wave limit $\sigma \to 0$:

1. **Definition and conservation of winding number** (topological property: $\sigma$-independent)
2. **Shape preservation after collision and phase shift** (property of scattering data: $\sigma$-independent)
3. **Finite parameter representation** $(N, \theta_0, \alpha)$ (structure of moduli space: $\sigma$-independent)
4. **Finiteness of energy and momentum** (due to the self-cancellation structure of each $2\pi$ rotation; Proposition 9.6a)

*Proof.* (1) follows from Propositions 9.1--9.2, (2) from Proposition 9.5, (3) from Proposition 8.4 and the subordination $\sigma = \cos\alpha$, and (4) from the energy self-cancellation structure of Proposition 9.6a. $\square$

**Corollary 9.1.** The rectangular wave description is not a "coarse approximation" of the soliton. It is a **rigorous limit that completely inherits** the topological and scattering structure of the soliton; what is lost is only the spread parameter $\sigma$ (details of the internal shape).

**Corollary 9.2 (Limits of applicability).** The conditions for the validity of this theorem depend on the implicit assumptions of the theory made explicit in Proposition 9.6a -- the fixed origin $\phi = 0$, the forcing of the terminus to vacuum values ($2n\pi$), and integer winding numbers. These assumptions hold only on infinite space $(-\infty,+\infty)$ or on a closed curve $S^1$. Extension to open strings (boundary conditions on a finite interval) requires the admission of half-integer winding numbers and a theory of soliton reflection at boundaries (Ghoshal-Zamolodchikov [23]), which lies outside the scope of this paper.

---

## 10. Conclusion

In this paper, we have provided an overview of the mathematical structure of the sine-Gordon equation and demonstrated the following:

1. Soliton solutions (kinks, anti-kinks, breathers) are **rigorously** characterized by a finite number of parameters $(N, \theta_0, \alpha)$ (§8).

2. All structural properties of solitons -- conservation of winding number, shape preservation during collisions and phase shift, finite parameter representation -- are preserved in the rectangular wave limit $\sigma \to 0$ (§9).

3. The finiteness of energy is based on the structure in which the energy of each $2\pi$ rotation is self-contained (the potential $V(\phi) = 1-\cos\phi$ returns from $0 \to 0$ over one period). A kink-anti-kink pair ($4\pi$ rotation) can undergo pair annihilation, and the permanent energy is only the contribution from the remaining kinks that cannot undergo pair annihilation.

4. Therefore, the transition from solitons to rectangular waves involves no loss of information. Among the infinite degrees of freedom of the continuous field, the physically meaningful information is concentrated in the finite number of parameters $(N, \theta_0, \alpha)$, and this concentration holds regardless of the soliton shape.

**Limits of applicability:** The above conclusions hold under the implicit assumptions of this theory (Proposition 9.6a).

- **Spatial constraints**: The phase at the origin is fixed at the vacuum value $\phi = 0$, and the phase at the terminus is forced to the vacuum value $\phi = 2n\pi$ (the winding number $N$ is integers only). Boundary reflections do not exist. These conditions are satisfied only on infinite space $(-\infty, +\infty)$ or on a closed curve $S^1$. Extension to open strings (boundary conditions on a finite interval) requires admitting half-integer winding numbers ($\phi = (2n+1)\pi$, potential maximum) and handling soliton reflection at boundaries, which lies outside the scope of this paper.
- **Collision constraints**: This theory is a completely integrable system, and soliton collisions are described as perfectly elastic collisions. To avoid violating momentum conservation, the velocity, shape, and energy of each soliton are individually conserved before and after the collision, and only the position (phase shift $\Delta x$) changes. This effectively means that solitons behave as free particles that do not interact with each other, and dissipation of energy, exchange of velocities, and inelastic channels (pair annihilation → radiation) do not occur within the scope of classical theory.

The rectangular wave description is not a "coarse approximation" of the soliton but rather a rigorous limit that completely inherits its topological and scattering structure. The fact that solutions of nonlinear wave equations are reduced to a finite number of discrete parameters reveals a deep correspondence between continuous field theories and discrete structures.

The above limits of applicability -- in particular, the inability to handle systems with different inertia parameters and the absence of velocity redistribution after collisions -- are overcome by the model presented in §12. The definition of transformation invariants for the sinusoidal-rectangular wave transformation is given in §11.

---

## 11. Sinusoidal-Rectangular Wave Transformation and Transformation Invariants

### 11.1 Purpose

In §9, we showed that soliton solutions preserve their structure even in the rectangular wave limit. In this section, we rigorously define the **forward transform** and **inverse transform** between sinusoidal and rectangular waveforms, and derive quantities invariant under the transformation (transformation invariants).

The invariants defined in this section are dimensionless quantities, and no identification with specific physical quantities is made.

### 11.2 Definition of Waveforms

We define periodic waveforms $\phi(\theta)$ on phase space as functions on the phase variable $\theta \in [0, 2\pi)$. Here $\theta$ is an index of processing steps, not time or distance.

**Sinusoidal waveform:**

$$\phi_{\sin}(\theta) = A_{\sin} \cdot \sin(\theta) \tag{W1}$$

**Rectangular waveform:**

$$\phi_{\text{rect}}(\theta) = A_{\text{rect}} \cdot \mathrm{sgn}(\sin(\theta)) \tag{W2}$$

Here $A_{\sin}$ and $A_{\text{rect}}$ are the maximum wave heights (amplitudes) of the respective waveforms.

### 11.3 Definition of Waveform Coefficient

**Definition 11.1 (Half-period area).** The half-period area $S[\phi]$ of a waveform $\phi(\theta)$ is defined as the integral of the absolute value of the waveform over a half-period (phase width $\pi$):

$$S[\phi] = \int_0^{\pi} |\phi(\theta)|\, d\theta \tag{W3}$$

We compute the half-period area for each waveform.

**Sinusoidal wave:**

$$S[\phi_{\sin}] = \int_0^{\pi} A_{\sin} \sin(\theta)\, d\theta = A_{\sin} \left[-\cos(\theta)\right]_0^{\pi} = 2A_{\sin} \tag{W4}$$

**Rectangular wave:**

$$S[\phi_{\text{rect}}] = \int_0^{\pi} A_{\text{rect}} \cdot 1\, d\theta = \pi \cdot A_{\text{rect}} \tag{W5}$$

**Definition 11.2 (Waveform coefficient).** The amplitude transformation coefficient $\kappa$ is defined as the factor that equalizes the half-period areas of the sinusoidal and rectangular waves:

$$\kappa = \frac{2}{\pi} \approx 0.6366 \tag{W6}$$

### 11.4 Forward Transform and Inverse Transform

**Definition 11.3 (Forward transform $T$).** The transform $T$ from a sinusoidal waveform to a rectangular waveform is defined as follows:

$$T: \phi_{\sin}(\theta) \mapsto \phi_{\text{rect}}(\theta), \qquad A_{\text{rect}} = \kappa \cdot A_{\sin} = \frac{2}{\pi} A_{\sin} \tag{W7}$$

That is, the waveform is rectangularized while simultaneously reducing the maximum wave height by a factor of $\kappa$.

**Definition 11.4 (Inverse transform $T^{-1}$).** The transform $T^{-1}$ from a rectangular waveform to a sinusoidal waveform is defined as follows:

$$T^{-1}: \phi_{\text{rect}}(\theta) \mapsto \phi_{\sin}(\theta), \qquad A_{\sin} = \frac{1}{\kappa} \cdot A_{\text{rect}} = \frac{\pi}{2} A_{\text{rect}} \tag{W8}$$

That is, the waveform is sinusoidalized while simultaneously increasing the maximum wave height by a factor of $1/\kappa$.

**Proposition 11.1 (Existence of the inverse).** $T^{-1} \circ T = T \circ T^{-1} = \mathrm{id}$ (identity transformation) holds.

*Proof.* The amplitude of $T^{-1}(T(\phi_{\sin}))$ is $(\pi/2) \cdot (2/\pi) \cdot A_{\sin} = A_{\sin}$, and the waveform returns to a sinusoidal wave. $\square$

### 11.5 Transformation Invariant $I_1$ (Quantitative Invariant)

**Definition 11.5.** The quantitative invariant $I_1[\phi]$ of a waveform $\phi$ is defined as the half-period area:

$$I_1[\phi] = S[\phi] = \int_0^{\pi} |\phi(\theta)|\, d\theta \tag{W9}$$

**Theorem 11.1 (Transformation invariance of $I_1$).** The quantitative invariant $I_1$ is invariant under the forward transform $T$ and the inverse transform $T^{-1}$:

$$I_1[T(\phi_{\sin})] = I_1[\phi_{\sin}]$$

*Proof.*

$$I_1[\phi_{\sin}] = 2A_{\sin} \qquad \text{(from W4)}$$

$$I_1[T(\phi_{\sin})] = \pi \cdot A_{\text{rect}} = \pi \cdot \frac{2}{\pi} A_{\sin} = 2A_{\sin} \qquad \text{(from W5, W7)}$$

Both are equal. $\square$

### 11.6 Transformation Invariant $I_2$ (Translational Invariant)

$I_1$ is the half-period area of the waveform (Definition 11.1). As an invariant expressing how much this quantity "moves" when the phase is displaced by one step, we define $I_2$. Through the introduction of a normalization coefficient $\kappa_\Delta$, the value of $I_2$ is preserved under waveform transformation.

**Definition 11.6 (Phase displacement per processing step).** The phase displacement $\Delta\phi$ per processing step (phase increment $\delta$) of a waveform $\phi(\theta)$ is defined as:

$$\Delta\phi = \phi(\theta + \delta) - \phi(\theta) \tag{W10}$$

**Definition 11.7 (Translational invariant).** The translational invariant $I_2[\phi]$ of a waveform $\phi$ is defined as the product of the phase displacement and the quantitative invariant:

$$I_2[\phi] = \kappa_\Delta \cdot \Delta\phi \cdot I_1[\phi] \tag{W11}$$

Here $\kappa_\Delta$ is a normalization coefficient that depends on the waveform, defined using the same waveform coefficient $\kappa = 2/\pi$ as for $I_1$:

- Sinusoidal wave: $\kappa_\Delta = 1$
- Rectangular wave: $\kappa_\Delta = 1/\kappa = \pi/2$

**Theorem 11.2 (Transformation invariance of $I_2$).** The translational invariant $I_2$ is invariant under the forward transform $T$ and the inverse transform $T^{-1}$.

*Proof.* The forward transform $T$ multiplies the amplitude by $\kappa$ and multiplies the normalization coefficient by $1/\kappa$. $I_1$ is invariant by Theorem 11.1. Since the phase displacement is proportional to the amplitude, it is multiplied by $\kappa$, but $\kappa_\Delta$ is multiplied by $1/\kappa$, so the product $\kappa_\Delta \cdot \Delta\phi$ is invariant.

$$I_2[T(\phi)] = \frac{1}{\kappa} \cdot (\kappa \cdot \Delta\phi) \cdot I_1 = \Delta\phi \cdot I_1 = I_2[\phi] \qquad \square$$

### 11.7 Summary of Transformation Invariants

**Theorem 11.3 (Invariants of the sinusoidal-rectangular wave transformation).** Under the forward transform $T$ (Definition 11.3) and inverse transform $T^{-1}$ (Definition 11.4) using the waveform coefficient $\kappa = 2/\pi$, the following quantities are invariant:

| Invariant | Definition | Value | Property |
|:------:|:-----|:---|:-----|
| $I_0$ | Winding number $N = [\phi(2\pi) - \phi(0)] / 2\pi$ | Integer | Topological invariant |
| $I_1$ | Half-period area $\int_0^{\pi} |\phi(\theta)|\, d\theta$ | $2A_{\sin}$ | Positive definite, additive |
| $I_2$ | $\kappa_\Delta \cdot \Delta\phi \cdot I_1$ | $\Delta\phi \cdot 2A_{\sin}$ | Signed, additive |

*Proof.* $I_0$: depends only on boundary values and is independent of waveform (Proposition 9.1). $I_1$: Theorem 11.1. $I_2$: Theorem 11.2. All are transformation invariant due to the cancellation of the waveform coefficient $\kappa$. $\square$

**Corollary 11.1.** The rectangular wave representation is not a "coarse approximation" of the sinusoidal wave representation. The transform $T$ is a bijection that preserves the three transformation invariants $I_0, I_1, I_2$. The values of the invariants computed on the rectangular wave coincide exactly with the values on the corresponding sinusoidal wave.

**Corollary 11.2.** The invariants are determined solely from the waveform parameters (amplitude $A$, winding number $N$) and do not depend on the internal structure of the waveform (whether sinusoidal or rectangular).

---

## 12. Elastic State Transformation Model on Phase Space

In §9--10, we identified the limits of applicability of the sine-Gordon soliton theory. In particular, (i) the requirement that all solitons have the same inertia parameter, and (ii) the absence of velocity redistribution after collisions (transmissive collisions), are insufficient for describing general elastic state transformations.

In this section, as an alternative model without these constraints, we formulate **elastic state transformations at arbitrary phases for oscillating systems with different inertia parameters**, predicated on the preservation of the transformation invariants $I_1, I_2$ defined in §11. All parameters of this model are abstract quantities on phase space, and no identification with specific physical quantities is made.

### 12.1 Definition of the Oscillating System

We define a harmonic oscillating system on phase space. The equation of motion for a system with inertia parameter $M$ and restoring coefficient $k$ is

$$M\ddot{x} + 2kx = 0 \tag{W12}$$

Here $x$ is the displacement from the equilibrium position, and $\dot{x}$ is the rate of change of displacement per processing step. The system has angular frequency $\omega = \sqrt{2k/M}$, and the general solution is

$$x(\theta) = A\cos(\omega\theta + \phi_0) \tag{W13}$$

Here $A$ is the amplitude, $\phi_0$ is the initial phase, and $\theta$ is the index of processing steps.

### 12.2 Configuration of the Two-Body System and Definition of State Transformation

Consider two independent oscillating systems with different inertia parameters $M_1, M_2$:

$$\ddot{x}_1 + \omega_1^2 x_1 = 0, \qquad \omega_1 = \sqrt{\frac{2k}{M_1}} \tag{W14}$$

$$\ddot{x}_2 + \omega_2^2 x_2 = 0, \qquad \omega_2 = \sqrt{\frac{2k}{M_2}} \tag{W15}$$

When $M_1 \neq M_2$, we have $\omega_1 \neq \omega_2$, and the two systems oscillate independently.

**Definition of state transformation.** In this paper, a "state transformation" is defined as a mathematical process in which the states $(x_1, \dot{x}_1)$ and $(x_2, \dot{x}_2)$ of both systems are coupled through two conservation conditions at a certain step $\theta = \theta_c$, instantaneously transitioning to new states $(x_1, \dot{x}_1')$ and $(x_2, \dot{x}_2')$. The positions $x_i$ do not change through the state transformation; only the rates of change $v_i = \dot{x}_i$ change according to (W20)--(W21).

**Remark.** The state transformation is defined as an abstract operation: "an instantaneous state transition in which the states of two oscillating systems are coupled through conservation conditions." The step $\theta_c$ at which the state transformation occurs is an input parameter given externally, and is not determined from within the model. This allows state transformations at arbitrary phases $\phi_1, \phi_2$ to be treated in a unified manner.

### 12.3 Elastic State Transformation at Arbitrary Phases

**State immediately before the state transformation:** The state of each system at step $\theta_c$ is expressed as:

$$x_{1,c} = A_1\cos\phi_1, \qquad v_1 = -A_1\omega_1\sin\phi_1 \tag{W16}$$

$$x_{2,c} = A_2\cos\phi_2, \qquad v_2 = -A_2\omega_2\sin\phi_2 \tag{W17}$$

Here $\phi_i = \omega_i\theta_c + \phi_{0,i}$ is the phase of each system at the state transformation step, and $A_i$ is the amplitude.

**Conservation conditions for elastic state transformation:**

$$M_1 v_1 + M_2 v_2 = M_1 v_1' + M_2 v_2' \tag{W18}$$

$$\frac{1}{2}M_1 v_1^2 + \frac{1}{2}M_2 v_2^2 = \frac{1}{2}M_1 v_1'^2 + \frac{1}{2}M_2 v_2'^2 \tag{W19}$$

Equation (W18) represents the conservation of the first-order form of the rate of change $\sum M_i v_i$, which has the same structure as $I_2$ (translational invariant: area $\times$ displacement) of §11.6. Equation (W19) represents the conservation of the second-order form of the rate of change $\sum \frac{1}{2}M_i v_i^2$, which has the same structure as $I_1$ (quantitative invariant: positive-definite additive conserved quantity) of §11.5. Both are mathematically isomorphic conservation conditions, and from these two conditions, the rates of change after the state transformation are uniquely determined.

**Theorem 12.1 (Elastic state transformation formula at arbitrary phases).** The rates of change after transformation for systems with inertia parameters $M_1, M_2$ undergoing state transformation at phases $\phi_1, \phi_2$ are given by:

$$v_1' = \frac{(M_1 - M_2)\,v_1 + 2M_2\,v_2}{M_1 + M_2} \tag{W20}$$

$$v_2' = \frac{(M_2 - M_1)\,v_2 + 2M_1\,v_1}{M_1 + M_2} \tag{W21}$$

*Proof.* Obtained immediately by solving the simultaneous equations (W18) and (W19). $\square$

### 12.4 Oscillation State After State Transformation

The state transformation is instantaneous, and positions do not change. After the state transformation, each system resumes oscillation with the position $x_{i,c}$ at the time of transformation and the new rate of change $v_i'$ as initial conditions:

$$x_1(\theta) = x_{1,c}\cos(\omega_1\tau) + \frac{v_1'}{\omega_1}\sin(\omega_1\tau) \qquad (\tau = \theta - \theta_c) \tag{W22}$$

$$x_2(\theta) = x_{2,c}\cos(\omega_2\tau) + \frac{v_2'}{\omega_2}\sin(\omega_2\tau) \tag{W23}$$

**Amplitude and phase after state transformation:**

$$A_i' = \sqrt{x_{i,c}^2 + \left(\frac{v_i'}{\omega_i}\right)^2} \tag{W24}$$

$$\phi_i' = \arctan\!\left(\frac{-v_i'}{\omega_i\, x_{i,c}}\right) \tag{W25}$$

The angular frequency $\omega_i$ is determined solely by the restoring coefficient and inertia parameter, and therefore does not change through the state transformation. What changes are only the amplitude $A_i$ and the phase $\phi_i$.

### 12.5 Proof of $I_1$ Conservation

**Proposition 12.1.** In an elastic state transformation at arbitrary phases, the total quantitative invariant of the system is conserved.

*Proof.* The total conserved quantity for each system is the sum of a position-dependent term and a rate-of-change-dependent term:

$$E_i = \frac{1}{2}M_i\omega_i^2\,x_{i,c}^2 + \frac{1}{2}M_i\,v_i^2 = k\,x_{i,c}^2 + \frac{1}{2}M_i\,v_i^2$$

Since the state transformation is instantaneous, the position-dependent term $k\,x_{i,c}^2$ is invariant before and after the transformation. The rate-of-change-dependent term is conserved by the conservation condition (W19). Therefore

$$E_1 + E_2 = E_1' + E_2' \qquad \square \tag{W26}$$

### 12.6 Characteristics of Amplitude and Phase Changes by State Transformation Phase

From the state transformation formulas (W20)--(W21) and the post-transformation amplitude and phase determination formulas (W24)(W25), the following characteristics appear depending on the phase $\phi_i$ at the transformation step:

| Transformation Phase | Position $x_{i,c}$ | Rate of Change $|v_i|$ | Conserved Quantity Distribution | Effect of State Transformation |
|:--|:--|:--|:--|:--|
| $\phi_i = \pi/2$ | $0$ | $A_i\omega_i$ (maximum) | Entirely in rate-of-change-dependent terms | Rates of change are fully redistributed; only amplitudes $A_i$ change |
| $\phi_i = 0$ | $A_i$ (maximum) | $0$ | Entirely in position-dependent terms | No exchange of rates of change occurs; state is invariant before and after |
| $\phi_i = \pi/4$ | $A_i/\sqrt{2}$ | $A_i\omega_i/\sqrt{2}$ | Equal distribution between position and rate of change | Partial redistribution of rates of change and phase shift occur |

**Remark.** The "invariance of state" at $\phi_i = 0$ is because the right-hand sides of (W20)(W21) are determined solely by $v_i$, and when $v_i = 0$, $v_i' = 0$ is the solution. In the general case $\phi_1 \neq \phi_2$, as long as the rates of change of both systems are not simultaneously zero, a non-trivial state transformation necessarily occurs.

### 12.7 Comparison with sine-Gordon Theory

| Property | sine-Gordon (§1--10) | This Model (§12) |
|:--|:--|:--|
| Different inertia parameters | Not possible (all kinks identical) | **Possible** ($M_1 \neq M_2$) |
| Rate-of-change redistribution after transformation | None (transmissive) | **Present** (state transformation formula) |
| State transformation at arbitrary phases | Undefined | **Fully defined** |
| Closed string assumption | Required (implicit) | **Not required** |
| Integer constraint on winding number | Required | **Not required** |
| Half-integer phases | Excluded ($\phi = \pi$ not allowed) | **Allowed** |
| Wave equation | $\phi_{tt} - \phi_{xx} + \sin\phi = 0$ | $\ddot{x}_i + \omega_i^2 x_i = 0$ |
| Hidden assumptions | 6 (Proposition 9.6a) | **None** (conservation conditions only) |
| Dissipation | Zero | Zero |
| Premises of formulation | Conservation of the first-order form ($\simeq I_2$) + conservation of the second-order form ($\simeq I_1$) + topological constraint + integrability + closed space + identical inertia | **$I_2$ conservation + $I_1$ conservation only** |

### 12.8 Equivalence with the Rectangular Wave Limit

In §9, we showed that solitons of the sine-Gordon equation preserve their structure in the rectangular wave limit. In this section, we prove that in the model of §12 as well, the results of state transformation are identical regardless of whether the waveform changes from sinusoidal (cosine) to rectangular.

**Introduction of waveform parameters.** We define the oscillation waveform as a family $x^\sigma(\theta)$ controlled by a transition width parameter $\sigma$ ($0 < \sigma \leq 1$):

- $\sigma = 1$: smooth sinusoidal (cosine) wave $x^1(\theta) = A\cos(\omega\theta)$
- $\sigma \to 0$: rectangular wave $x^0(\theta) = A\,\mathrm{sgn}(\cos(\omega\theta))$

Intermediate $\sigma$ controls the steepness of the transition region. For all $\sigma$, the amplitude is $A$ and the period is $T = 2\pi/\omega$ in common.

**Theorem 12.2 (Waveform independence of state transformation results).** When the position $x_{i,c}$ and total conserved quantity $E_i$ at the time of state transformation are given, the post-transformation rate of change $v_i'$, amplitude $A_i'$, and phase $\phi_i'$ do not depend on the waveform parameter $\sigma$.

*Proof.* We show this in three stages.

**Stage 1: The rate of change at the time of transformation is uniquely determined from the conserved quantity and position.**

The total conserved quantity of the system is preserved regardless of the waveform:

$$E_i = k\,x_i^2(\theta) + \frac{1}{2}M_i\,\dot{x}_i^2(\theta) = k\,A_i^2 \tag{W27}$$

This does not depend on the waveform. The position-dependent term is a function of position only, and the rate-of-change-dependent term is the remainder. Therefore, the rate of change at the transformation position $x_{i,c}$ is

$$v_i^2 = \frac{2(E_i - k\,x_{i,c}^2)}{M_i} = \omega_i^2(A_i^2 - x_{i,c}^2) \tag{W28}$$

and is **uniquely** determined. This does not depend on $\sigma$. $\checkmark$

**Stage 2: The state transformation formula depends only on the rates of change.**

The state transformation formulas (W20)--(W21) of Theorem 12.1 are algebraic relations involving only $M_1, M_2, v_1, v_2$, and contain no waveform information. Since Stage 1 showed that $v_i$ is $\sigma$-independent, $v_i'$ is also $\sigma$-independent. $\checkmark$

**Stage 3: The post-transformation amplitude and phase are also $\sigma$-independent.**

The total conserved quantity after transformation is

$$E_i' = k\,x_{i,c}^2 + \frac{1}{2}M_i\,v_i'^2 \tag{W29}$$

and the new amplitude is

$$A_i' = \sqrt{\frac{E_i'}{k}} = \sqrt{x_{i,c}^2 + \frac{v_i'^2}{\omega_i^2}} \tag{W30}$$

Since all of $x_{i,c}$, $v_i'$, $\omega_i$ are $\sigma$-independent, $A_i'$ is also $\sigma$-independent.

The initial phase $\phi_i' = \arctan(-v_i'/(\omega_i\,x_{i,c}))$ is likewise $\sigma$-independent. $\square$

**Corollary 12.1 (Rigorous nature of the rectangular wave limit).** The results of elastic state transformation in this model are **exactly identical** whether the oscillation is a smooth sinusoidal wave or a rectangular wave. The only information that is lost is the detail of the waveform (the shape of the transition region), while the redistribution of rates of change, conservation of conserved quantities, and changes in amplitude and phase do not depend on the waveform at all.

**Corollary 12.2.** This result reproduces, in a more general framework, the structural preservation of the rectangular wave limit demonstrated for the sine-Gordon equation in §9. While the proof of structural preservation in the sine-Gordon case required the assumptions of integrability and closed space, in this model the equivalence holds with **only the single condition of $I_1$ conservation** (equation (W27)).

### 12.9 Verification Summary

| Property | Sinusoidal | Rectangular | $\sigma$-dependence | Remarks |
|:----:|:------:|:------:|:--------------:|:----:|
| Rate of change at transformation $v_i$ | $\checkmark$ | $\checkmark$ | **Independent** | By (W28) |
| Post-transformation rate of change $v_i'$ | $\checkmark$ | $\checkmark$ | **Independent** | (W20)--(W21) |
| Post-transformation amplitude $A_i'$ | $\checkmark$ | $\checkmark$ | **Independent** | (W30) |
| Post-transformation phase $\phi_i'$ | $\checkmark$ | $\checkmark$ | **Independent** | (W25) |
| $I_1$ conservation | $\checkmark$ | $\checkmark$ | **Independent** | (W26) |
| $I_2$ conservation | $\checkmark$ | $\checkmark$ | **Independent** | (W18) |
| Different inertia parameters | $\checkmark$ | $\checkmark$ | **Independent** | $M_1 \neq M_2$ |
| Arbitrary phases | $\checkmark$ | $\checkmark$ | **Independent** | $\phi_1, \phi_2$ arbitrary |

**Theorem 12.3 (Structural preservation of the elastic state transformation model).** The process in which two oscillating systems with different inertia parameters $M_1 \neq M_2$ undergo elastic state transformation at arbitrary phases $\phi_1, \phi_2$ is completely described by the wave equations (W14)--(W15) and the state transformation formulas (W20)--(W21), and all of the following properties are preserved in the rectangular wave limit $\sigma \to 0$:

1. **Redistribution of rates of change** (state transformation formula: $\sigma$-independent)
2. **Conservation of $I_1$ and $I_2$** ((W26) and (W18): $\sigma$-independent)
3. **Determination of amplitude and phase** ((W30) and (W25): $\sigma$-independent)

*Proof.* A direct consequence of Theorem 12.2. $\square$

The only assumptions of this formulation are the following two conservation conditions:

- **$I_2$ conservation**: $M_1 v_1 + M_2 v_2 = M_1 v_1' + M_2 v_2'$
- **$I_1$ conservation**: $\frac{1}{2}M_1 v_1^2 + \frac{1}{2}M_2 v_2^2 = \frac{1}{2}M_1 v_1'^2 + \frac{1}{2}M_2 v_2'^2$

Topological constraints, integrability, compactness of space, integrality of winding numbers, periodicity of the potential, and the assumption of identical inertia parameters are all unnecessary, and the equivalence with the rectangular wave limit follows automatically from $I_1$ conservation.

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
