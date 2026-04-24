# Considerations on the Relationship between the 6-Dimensional Hypercube and Central Projection

## —— The Structure of Three Curvature Radii R, R₀, R₁ ——

**Author**: Noriaki Kihara (WF System Co., Ltd.)
**Date**: April 2026
**DOI**: https://doi.org/10.5281/zenodo.19731614

---

## 0. Purpose of This Paper



---

## 1. Preliminaries: Assumptions and Premises of This Paper

### 1.0 Declaration of the Unit System

In this paper, in order to treat the cosmological scale (curvature radius $R_0$) and the particle physics scale (proton charge radius $R_1$) within a unified unit system, we adopt **geometrized units** ($c = 1$, $G = 1$).

$$c = 1, \quad G = 1$$

In this unit system, $[\text{length}] = [\text{time}] = [\text{mass}]$, and all scales can be compared uniformly in terms of length (meters).

**$\hbar = 1$ is not adopted.** In natural units ($c = \hbar = 1$), which are conventional in particle physics, length and $1/\text{energy}$ are identified, rendering the geometric meaning of the ratio $R_0 / R_1$ between the cosmological and particle scales opaque. Since the purpose of this paper is to formulate this ratio geometrically, $\hbar$ is retained explicitly as a physical constant.

### 1.1 Assumption on Spatial Curvature $R_0$

In modern cosmology, there are multiple definitions for quantities characterizing the curvature of the universe. Below, curvature radii based on current observational values (Planck 2018 [5]) are shown in geometrized units with $c = 1$.

| Definition | Curvature Radius | Value ($c=1$) |
|:---|:---|:---|
| $\Lambda$-derived | $R_\Lambda = \sqrt{3/\Lambda}$ | $5.5 \times 10^{17}$ (17.4 billion light-years) |
| Hubble radius | $R_H = c/H_0$ | $4.6 \times 10^{17}$ (14.5 billion light-years) |
| Spatial curvature | $R_K = 1/\sqrt{\lvert K \rvert}$ | $> 7.3 \times 10^{18}$ (> 230 billion light-years) |

In this paper, we adopt the spatial curvature radius $R_K$ as $R_0$. The derivation is given below.

#### Derivation of the Spatial Curvature Radius $R_K$

The Gaussian curvature $K$ of space in the FLRW metric is defined in terms of the density parameter $\Omega_k$ and the Hubble constant $H_0$ as

$$K = -\,\Omega_k \frac{H_0^2}{c^2}$$

Here, $\Omega_k = 1 - \Omega_{\mathrm{total}}$ is a dimensionless quantity representing the deviation of the total energy density from the critical density, and the Planck 2018 observations report

$$\Omega_k = 0.001 \pm 0.002 \quad (68\% \text{ CL})$$

The central value is consistent with $\Omega_k = 0$ (spatially flat), and curvature has not been detected at the current observational precision.

Using the $2\sigma$ upper bound $\lvert \Omega_k \rvert < 0.004$, the upper bound on the Gaussian curvature is

$$\lvert K \rvert < 0.004 \times \frac{H_0^2}{c^2} = 0.004 \times 5.31 \times 10^{-53} \text{ m}^{-2} = 2.1 \times 10^{-55} \text{ m}^{-2}$$

The corresponding lower bound on the curvature radius is

$$R_K = \frac{1}{\sqrt{\lvert K \rvert}} > \frac{1}{\sqrt{2.1 \times 10^{-55}}} = 2.2 \times 10^{27} \text{ m}$$

Converting to $c = 1$ units:

$$R_0 \equiv R_K > \frac{2.2 \times 10^{27}}{c} = \frac{2.2 \times 10^{27}}{3.0 \times 10^{8}} = 7.3 \times 10^{18}$$

This corresponds to 230 billion light-years, which is approximately 50 times the radius of the observable universe (about 46.5 billion light-years). While space is locally nearly flat, the possibility of a finite curvature radius has not been excluded.

In this paper, we adopt this lower bound as an assumption:

$$R_0 = 7.3 \times 10^{18} \quad (c = 1)$$

### 1.2 Assumption on Particle Curvature $R_1$

In modern particle physics, there are also multiple definitions for quantities characterizing the "size" of elementary particles. Below, experimental values for representative particles are shown.

| Particle | Definition of Radius | SI Value | $c=1$ |
|:---|:---|:---|:---|
| Proton | Charge radius | $8.41 \times 10^{-16}$ m (0.841 fm) | $2.81 \times 10^{-24}$ |
| Neutron | Charge radius | $\sqrt{\langle r^2 \rangle} \approx 8.4 \times 10^{-16}$ m | $2.80 \times 10^{-24}$ |
| Charged pion | Charge radius | $6.6 \times 10^{-16}$ m (0.66 fm) | $2.20 \times 10^{-24}$ |
| Electron | Classical electron radius | $2.82 \times 10^{-15}$ m | $9.40 \times 10^{-24}$ |
| Electron | Compton wavelength $\bar{\lambda}_e = \hbar / m_e c$ | $3.86 \times 10^{-13}$ m | $1.29 \times 10^{-21}$ |
| Electron | Experimental upper bound (point particle) | $< 10^{-18}$ m | $< 3.3 \times 10^{-27}$ |
| Photon | Experimental upper bound (point particle) | $< 10^{-18}$ m | $< 3.3 \times 10^{-27}$ |
| Neutrino | Experimental upper bound (point particle) | $< 10^{-18}$ m | $< 3.3 \times 10^{-27}$ |

**Notes:**
- The proton charge radius is based on the PRad (2019) [6] result and the CODATA 2018 value $r_p = 0.8414 \pm 0.0019$ fm.
- The electron, photon, and neutrino have no detected internal structure at current experimental precision and are treated as point particles.
- The classical electron radius $r_e = e^2 / (4\pi\varepsilon_0 m_e c^2)$ is a scale for scattering cross sections and differs from the "size" of the particle.

In this paper, we adopt the **proton charge radius** as $R_1$, as the proton is the most fundamental particle whose internal structure has been experimentally confirmed.

$$R_1 = 2.81 \times 10^{-24} \quad (c = 1)$$

---

## 2. Modeling of Elementary Particles

**On the geometric interpretation of set signs.** In paper [4], set sign $0$ was defined as "coordinate value $= 0$," but this does not mean that the particle is located at the origin on that axis. Geometrically, set sign $0$ means that **the particle has no component in that axis direction**, i.e., it is **orthogonal to that axis**. Therefore, for a particle to have set sign $0$ on a given axis means that the state of the particle is independent of the value on that axis.

**On the limitations of the physical interpretation of set signs.** The set signs listed for each particle in this chapter are combinatorial structures defined in paper [4]. What the nonzero components ($+$, $-$, $\pm$) on each axis physically mean — that is, what physical quantities they correspond to and what values they take — remains an open problem at present. Below, we describe only the structural features of the set signs, and the physical interpretation is examined in later chapters.

### 2.1 Modeling of the Photon as a 6-Dimensional Hypercube

In paper [4] §10.10, the photon $\gamma$ is defined by the following set signs:

$$\gamma: \quad (x, y, z, t, R, Q) = (0, 0, 0, 0, +, 0)$$

The meaning of each axis is made explicit:

| Axis | Set Sign | Meaning |
|:---:|:---:|:---|
| $x$ | $0$ | Orthogonal to spatial axis 1: independent of the value on the $x$ axis |
| $y$ | $0$ | Orthogonal to spatial axis 2: independent of the value on the $y$ axis |
| $z$ | $0$ | Orthogonal to spatial axis 3: independent of the value on the $z$ axis |
| $t$ | $0$ | Orthogonal to the time axis: independent of the value on the $t$ axis |
| $R$ | $+$ | Curvature radius axis: has a positive component ($R > 0$) |
| $Q$ | $0$ | Color charge label: $Q = 0$ (color-neutral) |

The photon is a spin-1 boson ($n = 1$; only the $R$ axis is nonzero) and is orthogonal to all of the $x, y, z, t$ axes (independent of the values on these axes). It has $Q = 0$ (color-neutral).

What is noteworthy here is that the sole nonzero component characterizing the photon is on the $R$ axis, and its set sign is "$+$", i.e., $R > 0$. In the definition of set signs in paper [4] (§1.3.3), $+$ requires only that the coordinate value is positive, and **no specific value is prescribed**.

That is, the value on the $R$ axis is any positive real number satisfying $R > 0$, and $R$ is the **only free parameter** in the 6-dimensional hypercube model of the photon.

The cosmological curvature radius $R_0 = 7.3 \times 10^{18}$ established in §1.1 and the proton charge radius $R_1 = 2.81 \times 10^{-24}$ established in §1.2 both lie within the range $R > 0$, but their values differ by 42 orders of magnitude. Which value to assume for the photon's $R$ is an open problem in this model.

### 2.2 Modeling of the Graviton as a 6-Dimensional Hypercube

In paper [4] §10.10, the graviton $G$ is defined by the following set signs:

$$G: \quad (x, y, z, t, R, Q) = (0, 0, 0, \pm, \pm, 0)$$

The meaning of each axis is made explicit:

| Axis | Set Sign | Meaning |
|:---:|:---:|:---|
| $x$ | $0$ | Orthogonal to spatial axis 1: independent of the value on the $x$ axis |
| $y$ | $0$ | Orthogonal to spatial axis 2: independent of the value on the $y$ axis |
| $z$ | $0$ | Orthogonal to spatial axis 3: independent of the value on the $z$ axis |
| $t$ | $\pm$ | Time axis: has a nonzero component ($t \neq 0$) |
| $R$ | $\pm$ | Curvature radius axis: has a nonzero component ($R \neq 0$) |
| $Q$ | $0$ | Color charge label: $Q = 0$ (color-neutral) |

The graviton is a spin-2 boson ($n = 2$; both the $t$ and $R$ axes are nonzero) and is orthogonal to all of the $x, y, z$ axes. It has $Q = 0$ (color-neutral).

A noteworthy point in comparison with the photon is that the graviton has nonzero components on **both** the $t$ and $R$ axes. Whereas the photon has only the $R$ axis as a nonzero component, the graviton possesses a structure in which the $t$ and $R$ axes are coupled. That is, the mediator of gravity is involved in both time and curvature radius.

Also, as with the photon, the set sign on the $R$ axis is "$\pm$" ($R \neq 0$), and the specific value of $R$ is not prescribed.

### 2.3 Modeling of the Higgs Boson as a 6-Dimensional Hypercube

In paper [4] §10.10, the Higgs boson $H^0$ is defined by the following set signs:

$$H^0: \quad (x, y, z, t, R, Q) = (\pm, 0, 0, 0, 0, 0)$$

The meaning of each axis is made explicit:

| Axis | Set Sign | Meaning |
|:---:|:---:|:---|
| $x$ | $\pm$ | Spatial axis 1: has a nonzero component ($x \neq 0$) |
| $y$ | $0$ | Orthogonal to spatial axis 2: independent of the value on the $y$ axis |
| $z$ | $0$ | Orthogonal to spatial axis 3: independent of the value on the $z$ axis |
| $t$ | $0$ | Orthogonal to the time axis: independent of the value on the $t$ axis |
| $R$ | $0$ | Orthogonal to the curvature radius axis: independent of the value on the $R$ axis |
| $Q$ | $0$ | Color charge label: $Q = 0$ (color-neutral) |

The Higgs boson is a spin-0 boson ($n = 0$; orthogonal to all of the $t, R, Q$ axes) and has a signed length component ($\pm$) on only one of the $x, y, z$ axes. Since the $x, y, z$ axes are symmetric (paper [4], Proposition 1.1), which axis is nonzero is merely a choice of labeling convention.

### 2.4 Modeling of the Electron as a 6-Dimensional Hypercube

In paper [4] §10.10, the electron $e^-$ is defined by the following set signs:

$$e^-: \quad (x, y, z, t, R, Q) = (+, +, 0, 0, 0, 0)$$

The meaning of each axis is made explicit:

| Axis | Set Sign | Meaning |
|:---:|:---:|:---|
| $x$ | $+$ | Spatial axis 1: has a positive component ($x > 0$) |
| $y$ | $+$ | Spatial axis 2: has a positive component ($y > 0$) |
| $z$ | $0$ | Orthogonal to spatial axis 3: independent of the value on the $z$ axis |
| $t$ | $0$ | Orthogonal to the time axis: independent of the value on the $t$ axis |
| $R$ | $0$ | Orthogonal to the curvature radius axis: independent of the value on the $R$ axis |
| $Q$ | $0$ | Color charge label: $Q = 0$ (color-neutral) |

The electron is a spin-1/2 fermion (two of the $xyz$ axes are nonzero, $n = 0$), orthogonal to the $t$ and $R$ axes, with $Q = 0$ (lepton, color-neutral), and belongs to the first generation (the $x$ and $y$ axes are nonzero). Since $x = +$ and $y = +$, the electron has an unsigned area (positive $\times$ positive = positive) on the $xy$ plane.

### 2.5 Modeling of the Proton as a 6-Dimensional Hypercube

The proton $p$ is a composite particle composed of three quarks $u_r, u_g, d_b$. Based on paper [4] §10.10, the set signs of each constituent quark are listed below.

**Up quark $u_r$ (color charge $r$):**

$$u_r: \quad (x, y, z, t, R, Q) = (+, +, 0, 0, 0, r)$$

| Axis | Set Sign | Meaning |
|:---:|:---:|:---|
| $x$ | $+$ | Spatial axis 1: has a positive component ($x > 0$) |
| $y$ | $+$ | Spatial axis 2: has a positive component ($y > 0$) |
| $z$ | $0$ | Orthogonal to spatial axis 3: independent of the value on the $z$ axis |
| $t$ | $0$ | Orthogonal to the time axis: independent of the value on the $t$ axis |
| $R$ | $0$ | Orthogonal to the curvature radius axis: independent of the value on the $R$ axis |
| $Q$ | $r$ | Color charge label: red |

**Up quark $u_g$ (color charge $g$):**

$$u_g: \quad (x, y, z, t, R, Q) = (+, +, 0, 0, 0, g)$$

| Axis | Set Sign | Meaning |
|:---:|:---:|:---|
| $x$ | $+$ | Spatial axis 1: has a positive component ($x > 0$) |
| $y$ | $+$ | Spatial axis 2: has a positive component ($y > 0$) |
| $z$ | $0$ | Orthogonal to spatial axis 3: independent of the value on the $z$ axis |
| $t$ | $0$ | Orthogonal to the time axis: independent of the value on the $t$ axis |
| $R$ | $0$ | Orthogonal to the curvature radius axis: independent of the value on the $R$ axis |
| $Q$ | $g$ | Color charge label: green |

**Down quark $d_b$ (color charge $b$):**

$$d_b: \quad (x, y, z, t, R, Q) = (+, -, 0, 0, 0, b)$$

| Axis | Set Sign | Meaning |
|:---:|:---:|:---|
| $x$ | $+$ | Spatial axis 1: has a positive component ($x > 0$) |
| $y$ | $-$ | Spatial axis 2: has a negative component ($y < 0$) |
| $z$ | $0$ | Orthogonal to spatial axis 3: independent of the value on the $z$ axis |
| $t$ | $0$ | Orthogonal to the time axis: independent of the value on the $t$ axis |
| $R$ | $0$ | Orthogonal to the curvature radius axis: independent of the value on the $R$ axis |
| $Q$ | $b$ | Color charge label: blue |

All three quarks are orthogonal to the $R$ axis (independent of the value on the $R$ axis).

Focusing on the sign of the area, the up quark $u$ has $x = +$ and $y = +$, so it has an unsigned area (positive $\times$ positive = positive) on the $xy$ plane, while the down quark $d$ has $x = +$ and $y = -$, so it has a signed area (positive $\times$ negative = negative). Among the three constituent quarks of the proton $p = u_r + u_g + d_b$, only the down quark has a negative area.

### 2.6 Comparison Focusing on the $R$ Axis

Comparing the particles discussed so far with a focus on the $R$ axis:

| Particle | Spin | $R$ Axis | $t$ Axis | Spatial Axes |
|:---|:---:|:---:|:---:|:---|
| Photon $\gamma$ | 1 | $+$ | $0$ | All $0$ |
| Graviton $G$ | 2 | $\pm$ | $\pm$ | All $0$ |
| Higgs $H^0$ | 0 | $0$ | $0$ | $x = \pm$ |
| Electron $e^-$ | 1/2 | $0$ | $0$ | $x = +, y = +$ |
| Quarks $u, d$ | 1/2 | $0$ | $0$ | $x = +, y = \pm$ |

Particles that have components on the $R$ axis (the photon and graviton) are all massless and orthogonal to all spatial axes. On the other hand, particles that are orthogonal to the $R$ axis (the Higgs boson, electron, and quarks) have mass and possess nonzero components on the spatial axes.

---

## 3. 

---

## References

[1] Noriaki Kihara, "Geometric Formulation of 4-Dimensional Spacetime via Gnomonic Orthoprojection," 2025.
[2] Noriaki Kihara, "Construction of Subjective Space via Central Projection and Formulation of Interactions in Phase Space," 2026.
[3] Noriaki Kihara, "Global Coverage of Central Projection in Discrete Spaces and the Number-Theoretic Necessity of 5-Dimensional Background Space," 2026.
[4] Noriaki Kihara, "Set-Theoretic Structure of 6-Dimensional Hypercubes and Their Combinatorial Properties," 2026.
