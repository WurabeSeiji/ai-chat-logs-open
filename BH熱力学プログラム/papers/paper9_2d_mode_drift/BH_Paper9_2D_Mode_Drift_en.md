# BH Thermodynamics Program Paper 9: Geometric Origin of the 0.036 Drift in the α Self-Consistency Equation
## ― Physical Interpretation via 2D Surface Vibration Modes in a 4D Hyperball

**Author**: Noriaki Kihara
**Affiliation**: WF System Co., Ltd.
**ORCID**: [0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020)
**Version**: v1
**Date**: May 20, 2026
**License**: CC BY 4.0
**Concept DOI**: [10.5281/zenodo.20319436](https://doi.org/10.5281/zenodo.20319436)
**v1 Version DOI**: [10.5281/zenodo.20319437](https://doi.org/10.5281/zenodo.20319437)
**Zenodo page**: https://zenodo.org/records/20319437

---

## Nature of this Paper

**This is an observational/interpretive paper, not a first-principles proof paper.**

We present one interpretation of the physical origin of the coefficient $\pi^2/2$ in the α identity $\alpha^{-1} = 137 + (\pi^2/2)\alpha$ of Paper 7 [BH7], based on the uncertainty principle of standard quantum mechanics and dimensional analysis of vibration modes.

This paper does **not** claim:

- That $\pi^2/2$ admits no interpretation other than the one given here
- That the high-precision residual of the observed value $\alpha^{-1} = 137.035999...$ is fully derived by the present interpretation
- To exclude other possible interpretations
- To "prove", "strengthen", or "extend" the claims of Paper 7

What this paper offers is one way of reading. Paper 7's claims do not depend on this paper. The numerical fact observed in Paper 7 (that $\alpha^{-1} = 137 + (\pi^2/2)\alpha$ holds at 8.7 ppb precision) is unchanged by whether one adopts this interpretation.

## Scope

The scope of this paper:

- **Subject**: $\alpha^{-1}(Q^2 \to 0) = 137.036$ (geometric origin of α in the Thomson limit)
- **Method**: Standard QM uncertainty principle + dimensional analysis of vibration modes in 4D space
- **Conclusion**: The 0.036 drift can be interpreted as the zero-point expectation value of 2D surface vibration modes

This paper does **not** address:

- The high-energy running of α ($\alpha^{-1}(M_Z) \approx 127.95$)
- Geometric rederivation of vacuum polarization from Standard Model lepton/quark loops
- Unification of coupling constants (α with $\alpha_s$, $\alpha_w$, etc.)
- Connection with gravity

These are already well-described by Standard Model vacuum polarization theory, and the present geometric interpretation takes the position of providing a **boundary condition (value at $Q^2 = 0$)** to that framework. The geometric reformulation of the high-energy side is left as future work (§8).

## Abstract

In Paper 7 [BH7], the self-consistent equation using unit-cube packing in a 4D ball,

$$\alpha^{-1} = 137 + \frac{\pi^2}{2}\alpha,$$

was reported to agree with the observed value $\alpha^{-1}(0) = 137.035999084$ at a relative precision of 8.7 ppb. This paper presents the following physical interpretation of the structure of this equation:

1. **Although α is dimensionless, its physical contributions appear through area-dimensional quantities** (§2): Since the scattering cross section $\sigma_T \propto \alpha^2$, the physical role of α is mediated by area.
2. **Classifying 4D space vibration modes by dimension, only 2D face modes can contribute to α under isotropic averaging** (§3): 0D/1D/3D/4D modes carry directionality and vanish under averaging, while 2D area remains as a scalar quantity.
3. **The position phase space measure of 2D faces of 137 hypercubes distributed in a 4D unit ball equals $\int_{B_4(1)} dV = \pi^2/2$** (§4): While this coincides numerically with $V_4(1)$, in our interpretation it is read as "the integral of position degrees of freedom of 2D faces".
4. **The W7 self-consistent equation $\alpha^{-1} = 137 + (\pi^2/2)\alpha$ can be interpreted as a self-consistent correction from zero-point vibrations of 2D face modes of the 137 hypercubes** (§5).
5. **This structure is identical to the plaquette action of Wilson lattice gauge theory** (§6): The isomorphism shown in Paper 8 [BH8] is realized at the geometric level.
6. **There exists a theoretical lower bound on the observed value of α, arising from the width of the 2D mode distribution** (§7).

Our interpretation is completed within the framework of standard quantum mechanics and does not require new assumptions (discrete spacetime, complex integer lattices, etc.).

---

## §1 Introduction

Paper 7 [BH7] presented the α self-consistent equation

$$\alpha^{-1} = N(1) + V_4(1)\alpha = 137 + \frac{\pi^2}{2}\alpha$$

using the 137 unit 4D cubes (tesseracts) inscribed in the 4D ball $B_4(R=3)$ and the volume of the 4D unit ball $V_4(1) = \pi^2/2$. This equation agrees with the observed value $\alpha^{-1}(0) = 137.035999084(21)$ at a precision of 8.7 ppb.

The structure of this self-consistent equation left the following unresolved questions:

**Question 1**: Why does $\pi^2/2$ appear as a coefficient of α?
- Numerically it coincides with $V_4(1)$, but it is physically unclear why the volume of a 4D unit ball should appear
- Is it a mere mathematical identity, or does it carry physical meaning?

**Question 2**: Geometric origin of the 0.036 drift
- The geometric structure giving the integer 137 (W6 paper) is clear
- But the origin of "0.036" was only explained as a closure of the self-consistency

**Question 3**: Physical content of the isomorphism with Wilson lattice gauge theory shown in W8 [BH8]
- Structural correspondence is established
- But the physical content of the correspondence (why it holds) remains unclear

This paper presents the following interpretation in response to these questions:

> **α is physically an area-dimensional quantity ($\sigma \propto \alpha^2$), and only 2D face modes can contribute to α among the 4D space vibration modes. The position phase space measure of 2D face vibration modes of the 137 hypercubes gives $\pi^2/2$, which physically explains the coefficient structure of the W7 self-consistent equation.**

This interpretation takes the standard QM uncertainty principle as its starting point and requires no new assumptions.

---

## §2 Dimensional Analysis of α and Area Dimensionality

### 2.1 α itself is dimensionless

The fine-structure constant α is defined in SI units as

$$\alpha = \frac{e^2}{4\pi\epsilon_0 \hbar c} \approx 7.297 \times 10^{-3}$$

and is dimensionless.

### 2.2 Physical contributions appear through α²

However, in many cases the quantities through which α is physically involved appear as $\alpha^2$:

**Thomson scattering cross section**:
$$\sigma_T = \frac{8\pi}{3} r_e^2 = \frac{8\pi}{3} \left(\alpha \cdot \frac{\hbar}{mc}\right)^2 \propto \alpha^2 \cdot [\text{length}]^2$$

**Rutherford scattering differential cross section**:
$$\frac{d\sigma}{d\Omega} \propto \alpha^2 \cdot [\text{length}]^2$$

**Bohr radius**:
$$a_0 = \frac{\hbar}{\alpha mc} \propto \alpha^{-1} \cdot [\text{length}]$$

**Classical electron radius**:
$$r_e = \alpha \cdot \frac{\hbar}{mc} \propto \alpha \cdot [\text{length}]$$

Thus, α is a quantity directly linked to **cross section (area dimension)**.

### 2.3 Correspondence with Geometric Structures in 4D Space

When considering geometric structures in 4D space to which α naturally couples:

- σ (cross section) ∝ α² ∝ area
- → α is a quantity that couples to **2-dimensional structures**

This is consistent with the fact that in Wilson lattice gauge theory, the gauge coupling appears as the action coefficient of plaquettes (2D faces, 2-forms).

### 2.4 Starting Point of This Paper

From the above:

> **When considering contributions to α from vibration modes in 4D space, the area dimensionality of α restricts the contributing modes to those with "2-dimensional area".**

This is the central observation of this paper.

---

## §3 Dimensional Classification of Vibration Modes in 4D Space

### 3.1 Vibration Modes of 137 Hypercubes and the Circumscribed 4-Sphere

The system treated in Paper 7:
- 4D ball $B_4(3)$ of radius $R = 3$
- 137 unit 4D cubes (tesseracts) packed inside

This system has the following vibration modes:

**Internal structure of one tesseract**:

| Structural element | Count | Independent vibration modes |
|---|---|---|
| Vertex (0-cell) | 16 | 16 × 4 = 64 directions |
| Edge (1-cell) | 32 | Stretching 32, transverse, etc. |
| Face (2-cell, square) | 24 | Bulging, twisting |
| Solid (3-cell, cube) | 8 | Volume vibration |
| 4-cell (whole tesseract) | 1 | Whole breathing mode |

**Across all 137 tesseracts**:
- Each structural element's vibration multiplied by 137
- Interaction modes between adjacent cubes
- Phonon-like spectrum as a whole

**Circumscribed 4-sphere $S^3$ (radius 3)**:

| Mode | Description |
|---|---|
| Radial breathing | Vibration of R = 3 |
| Spherical harmonics $Y_{\ell m n}$ on $S^3$ | Angular vibrations ($\ell = 0, 1, 2, \ldots$) |
| 4-axis distortion | Sphere → ellipsoid deformation |

### 3.2 Survival Conditions under Isotropic Averaging

Among these modes, what matters for contributions to α is the **survival condition under isotropic averaging**.

α is dimensionless and scalar, invariant under SO(4) symmetry. Therefore contributions to α are limited to **modes that do not vanish under isotropic averaging (averaging over 4-axis directions)**.

| Mode dimension | Mathematical property | Result of isotropic averaging |
|---|---|---|
| 0D (point) | Vector position fluctuation | **Zero** (translation symmetry) |
| 1D (line segment) | Has direction | **Zero** (line symmetry, sign flips on reversal) |
| **2D (face)** | **Area (scalar quantity, no direction)** | **Survives** |
| 3D (volume, oriented) | Volume element $dV = dx \wedge dy \wedge dz$ is pseudoscalar (has orientation) | **Zero** (sign flips on reversal) |
| 4D (4-volume) | 4-volume element is pseudoscalar | **Zero** |

**Key fact**: **Area is a scalar quantity**, and area is invariant under flipping the orientation of the face. In contrast, length and volume (oriented volume elements) change sign under reversal of orientation.

### 3.3 Selection Principle

From the above:

> **Among vibration modes in 4D space, only 2D face modes can give significant contributions to α under isotropic averaging.**

This is consistent with the area dimensionality of α (§2).

### 3.4 Consistency with Wilson Lattice Gauge Theory

This conclusion agrees with the structure of **2-form gauge fields (plaquette action)** in Wilson lattice gauge theory:

- 0-form (scalar field): gauge-invariant point-like structure
- 1-form (vector potential $A_\mu$): gauge-dependent
- **2-form (field strength $F_{\mu\nu}$): gauge-invariant, plaquette action $S = \beta \sum_p \text{Re}\,\text{tr}(U_p)$**
- 3-form, 4-form: higher-order structures

That α "couples to faces" is a standard fact in Wilson formalism.

---

## §4 Position Phase Space of 2D Face Modes of 137 Hypercubes

### 4.1 Position Degrees of Freedom of 2D Face Modes

To fully specify a 2D face vibration mode requires:

- **Position**: Where the 2D face is located in the 4D ball $B_4(R=3)$
- **Direction**: One of 6 possible 2-plane directions in 4D space ($\binom{4}{2} = 6$)
- **Amplitude**: Excitation level of the vibration

### 4.2 Integration of Position Degrees of Freedom

For a fixed direction and amplitude, integrating the position degrees of freedom of 2D faces within the 4D unit ball (scale-normalized) gives:

$$\text{2D face position phase space} = \int_{B_4(1)} dV = \frac{\pi^2}{2}$$

This coincides numerically with **$V_4(1) = \pi^2/2$**.

### 4.3 Interpretation: Reading of Physical Meaning

| Interpretation | Physical meaning of $\pi^2/2$ |
|---|---|
| Old interpretation (W7) | Volume $V_4(1)$ of 4D unit ball |
| **This paper** | **Measure of position degrees of freedom for placing 2D face modes in 4D unit ball** |

The two are **numerically identical, but physically different**. Under our interpretation, $\pi^2/2$ has a necessary reason to appear as the coefficient of α: since α couples to 2D face modes, the volume of their position phase space appears as the geometric coefficient of the coupling constant.

### 4.4 Direction and Amplitude Contributions

The direction contribution from 6 possible 2-planes may appear as a combinatorial factor separately. In this paper, we treat the position integral as the dominant factor, and consider the directional factor as an $O(1)$ correction even if numerically effective.

The amplitude appears as α itself (excitation probability amplitude):

$$\text{Collective contribution of 2D face modes} = \underbrace{\frac{\pi^2}{2}}_{\text{position phase space}} \times \underbrace{\alpha}_{\text{amplitude per face}} = \frac{\pi^2}{2} \alpha$$

This physically explains the α coefficient term in the W7 self-consistent equation.

---

## §5 Physical Interpretation of the W7 Self-Consistent Equation

### 5.1 Structure of the Equation

The self-consistent equation of Paper 7 [BH7]:

$$\alpha^{-1} = 137 + \frac{\pi^2}{2}\alpha$$

### 5.2 Each Term under Our Interpretation

| Term | Mathematical content | Physical interpretation (this paper) |
|---|---|---|
| **137** | Number of unit cubes packed in 4D ball $B_4(3)$ | **Contribution when all vibration modes are in the ground state (lowest energy)**: integer-theoretic fact of 4D ℤ⁴ lattice |
| **$\pi^2/2$** | $V_4(1)$, or in our interpretation, 2D face mode position phase space | **Position degrees of freedom for placing 2D face modes of 137 hypercubes in the 4D unit ball** |
| **α** (in coefficient term) | Coupling constant | **Excitation amplitude of each 2D face mode** |
| **$(\pi^2/2)\alpha$** | 0.036 (numerically) | **Collective contribution of zero-point vibrations of 2D face modes** (from standard QM uncertainty principle) |

### 5.3 Connection with Standard QM Uncertainty Principle

Each 2D face mode has the following zero-point energy under the harmonic oscillator approximation:

$$E_{\text{zero-point}}^{(i)} = \frac{1}{2}\hbar\omega_i$$

where $\omega_i$ is the eigen-frequency of the i-th 2D face mode.

The collective contribution of all 2D face modes is:

$$\sum_i \frac{1}{2}\hbar\omega_i = (\text{position phase space measure}) \times (\text{amplitude density}) = \frac{\pi^2}{2} \cdot \alpha$$

Here the α on the right side represents the self-consistency condition that the average excitation amplitude of each mode equals α itself.

### 5.4 Self-Consistency

From the structure that α is the amplitude of 2D face modes and that the zero-point contribution of 2D face modes provides a correction to $\alpha^{-1}$, the self-consistent equation

$$\alpha^{-1} = 137 + \frac{\pi^2}{2}\alpha \quad \Leftrightarrow \quad \frac{\pi^2}{2}\alpha^2 + 137\alpha - 1 = 0$$

is naturally derived.

### 5.5 Numerical Agreement

The positive root of this equation:

$$\alpha^{-1} = 137.036010988\ldots$$

The relative error with the CODATA 2018 value $\alpha^{-1}(0) = 137.035999084(21)$ is **8.7 ppb**. This is the fact already observed in Paper 7, and our interpretation provides a **physical framework** explaining this numerical agreement.

---

## §6 Strengthening of the Isomorphism with Wilson Lattice Gauge Theory

Paper 8 [BH8] showed that there is a structural correspondence between the W7 self-consistent equation and Wilson lattice gauge theory based on Schläfli duality and $B_4$ equivariance. Our 2D face mode interpretation clarifies the **physical content** of this correspondence.

### 6.1 Correspondence

| Wilson lattice gauge theory | W7 self-consistency (this paper's interpretation) |
|---|---|
| Coupling constant $\beta = 1/g^2$ | $\alpha^{-1}$ |
| Plaquette (2-form structure) | 2D face vibration mode |
| Number of plaquettes $\sum_p$ | $\pi^2/2$ (position phase space measure of 2D faces) |
| Plaquette action $\beta \sum_p \text{Re}\,\text{tr}(U_p)$ | $(\pi^2/2)\alpha$ |
| Cell (4D hard core) | 137 (number of unit cube packing) |
| Total action | $\alpha^{-1} = 137 + (\pi^2/2)\alpha$ |

### 6.2 Clarification of Physical Content

The isomorphism in Paper 8 was established as a **structural correspondence**, but with this paper it is **physically positioned** as:

> **The W7 self-consistent equation is the realization of the plaquette action structure of Wilson lattice gauge theory on the 4D ball geometry.**

### 6.3 Implications

With this:
- The correspondence between W7 and W8 is elevated from "mathematical isomorphism" to "physical isomorphism"
- α coupling to "2D faces" aligns with the natural interpretation of Wilson formalism
- A roadmap appears for **geometric reconstruction** of vacuum polarization in standard QFT

---

## §7 Theoretical Lower Bound on Observed Width of α

### 7.1 Distribution Width of Zero-Point Vibrations

Each 2D face mode's zero-point vibration has a finite distribution width from the standard QM uncertainty principle:

$$\Delta n_i^{\text{zero}} = \frac{1}{2}$$

(zero-point fluctuation of the photon number operator)

### 7.2 Propagation to α

The fluctuation in the collective contribution of all 2D face modes is:

$$(\Delta \alpha^{-1})^2 \approx \sum_i \left(\frac{\partial \alpha^{-1}}{\partial n_i}\right)^2 (\Delta n_i)^2$$

A concrete calculation is beyond the scope of this paper (requires determination of mode eigen-frequencies), but structurally:

> **There exists a theoretical lower bound on the observed value of α, arising from the width of the 2D face mode distribution.**

### 7.3 Correspondence with CODATA Measurement

$$\alpha^{-1}(0) = 137.035999084(21)$$

The CODATA uncertainty $(21) \approx 2.1 \times 10^{-8}$ is experimental precision, and whether our theoretical lower bound is even smaller requires further consideration of mode frequencies.

### 7.4 Implications

By our interpretation:
- α is not a "sharp point" but the "center value of a distribution"
- The observed 8.7 ppb deviation (W7) leaves room for higher-order mode contributions
- Future precision measurements may make our theoretical lower bound observable

---

## §8 Open Problems: Extension to High Energies

Our interpretation **addresses only the Thomson limit $\alpha^{-1}(Q^2 \to 0) = 137.036$**.

### 8.1 Running of α in Standard QED

In the Standard Model, α runs with the energy scale $Q^2$:

- $\alpha^{-1}(0) = 137.036$ (Thomson limit)
- $\alpha^{-1}(M_Z^2) \approx 127.95$ (Z boson mass scale)
- The difference $\approx 9.08$ units is explained by vacuum polarization from lepton/quark loops

This running is **well-described** by precision measurements of standard QED, and this paper does not compete with it.

### 8.2 Position of This Paper

Our geometric interpretation provides a physical origin for the **boundary condition (value at $Q^2 = 0$)** of α.
Running is left to standard QED.

The two are hierarchically complementary:

| Layer | Content | Responsibility |
|---|---|---|
| 1 | 137 (integer, geometric invariant) | Paper 6 [BH6]: integer theory of packing |
| 2 | $(\pi^2/2)\alpha \approx 0.036$ | **This paper: 2D face modes** |
| 3 | $\Delta\alpha_{SM}(Q^2)$ | Standard QED vacuum polarization |

### 8.3 Future Tasks

The relationship between layer 3 (standard QED vacuum polarization) and our geometric interpretation leaves the following open questions:

1. **Additional modes activated at high energies**: 3D and 4D vibration modes may start contributing to α at high energies
2. **Geometric counterparts of SM loops**: Whether each of electron loop, muon loop, quark loop corresponds to a specific geometric mode
3. **Geometric derivation of $\alpha^{-1}(M_Z) \approx 128$**: Whether our framework can be extended to high energies
4. **α behavior at GUT scale**: Geometric prediction of coupling constants at grand unification scale

The geometric reformulation of these is left as future research. **The claim of this paper (geometric origin of the Thomson limit) holds independently of the mechanism of running**.

---

## §9 Conclusion

This paper presented the following as a physical interpretation of the α self-consistent equation $\alpha^{-1} = 137 + (\pi^2/2)\alpha$ in Paper 7 [BH7]:

1. **α is dimensionless but physically area-dimensional** ($\sigma \propto \alpha^2$)
2. **Among 4D space vibration modes, only 2D face modes can contribute to α under isotropic averaging**
3. **Position phase space measure of 2D faces of 137 hypercubes distributed in the 4D unit ball = $\pi^2/2$**
4. **The W7 self-consistent equation can be interpreted as a self-consistent contribution of zero-point vibrations of 2D face modes**
5. **This structure is identical to Wilson lattice gauge theory plaquette action** (physical content of W8 isomorphism)
6. **The observed value of α has a theoretical lower bound arising from the 2D mode distribution width**

Our interpretation is completed within the framework of the standard QM uncertainty principle, and requires no new assumptions (discrete spacetime, extra dimensions, etc.).

**The scope of this paper** is only the Thomson limit of α. The high-energy running is left to standard QED, and its geometric reformulation is left as future tasks (§8).

### Position of This Paper

Relationship among Paper 7 [BH7], Paper 8 [BH8], and this paper:

| Paper | Content | Area |
|---|---|---|
| BH7 (α identity) | Discovery of $\alpha^{-1} = 137 + (\pi^2/2)\alpha$ | Algebraic observation |
| BH8 (Wilson isomorphism) | Proof of structural correspondence | Mathematical correspondence |
| **This paper (BH9)** | **Physical interpretation via 2D face modes** | **Physical content** |

With this trilogy, the geometric origin of α is positioned in three layers: **observation → structure → physical content**.

---

## References

[BH6] Kihara, N. (2026). *Synthesis and the 4+1 to 3+1 Reduction*. Zenodo. Concept DOI: 10.5281/zenodo.19837597.

[BH7] Kihara, N. (2026). *A Geometric Identity for the Fine-Structure Constant: From the 4D Unit Ball Volume and its Cube-Packing Deficit*. Zenodo. Concept DOI: 10.5281/zenodo.19869266.

[BH8] Kihara, N. (2026). *Chain Complex Structure on the 4D Hypercubic Lattice: Structural Correspondence between Kihara Cube-Packing and Wilson Lattice Gauge Theory*. Zenodo. Concept DOI: 10.5281/zenodo.19880467.

[BH7-Supp] Kihara, N. (2026). *Paper 7 Supplement: Geometric Observation on the Second-Order Correction Term of the α Identity*. Zenodo. Concept DOI: 10.5281/zenodo.19933729.

[Wilson1974] Wilson, K. G. (1974). *Confinement of quarks*. Phys. Rev. D **10**, 2445.

[CODATA2018] Tiesinga, E. et al. (2021). *CODATA Recommended Values of the Fundamental Physical Constants: 2018*. Rev. Mod. Phys. **93**, 025010.

---

## Revision History

- **v1 (May 20, 2026)**: Initial version
