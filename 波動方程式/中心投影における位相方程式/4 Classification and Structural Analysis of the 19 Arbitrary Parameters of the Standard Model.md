# Classification and Structural Analysis of the 19 Arbitrary Parameters of the Standard Model

## — Evaluation Based on Structural Requirements for the Theory of Everything —

**Author:** Noriaki Kihara (WF System Co., Ltd.)

**Date:** April 2026

**DOI:** [10.5281/zenodo.19604965](https://doi.org/10.5281/zenodo.19604965)

**Prerequisite Paper:** [Formulation of Structural Requirements That a Theory of Everything Must Satisfy](https://doi.org/10.5281/zenodo.19601592)

---

## 0. Purpose of This Paper

This paper classifies and analyzes the 19 arbitrary parameters (free parameters) contained in the Standard Model of particle physics from the perspective of structural requirements formulated in a prior work.

Specifically, we address the following questions:

1. Is each parameter structurally derivable, or is it truly arbitrary?
2. Which parameters should a Theory of Everything (TOE) internally generate (= TOE KPIs)?
3. What are the dependency relationships among the parameters?

---

## 1. Parameter Catalog

We present the full picture of parameters involved in the fundamental theories of modern physics (Standard Model + General Relativity). By convention, the Standard Model counts "19 arbitrary parameters," but for evaluating a Theory of Everything, we include all parameters as verification targets, including those treated as "fundamental constants."

The values of each parameter are those in the low-energy regime currently confirmed by experiment (roughly the eV to TeV scale, i.e., the range accessible to atomic and particle experiments).

### 1.0 Fundamental Physical Constants (Currently Treated as "Constants")

The following are treated as constants in current physics and are eliminated in natural units ($c = \hbar = k_B = 1$), so they are not counted among the "19 arbitrary parameters." However, among candidates for grand unified theories and quantum gravity, there exist hypotheses that these vary over cosmic history or are derivable from more fundamental quantities. We include them as verification targets for a Theory of Everything.

These constants are classified by the dimensional composition of their constituent quantities. This classification clarifies which layer of physics each constant belongs to.

#### Group I: Pure Spacetime ($L, T$)

| # | Parameter | Symbol | Value | Dimensions | Conversion |
|---|-----------|--------|-------|------------|------------|
| P1 | Speed of light | c | 299,792,458 m/s (definition) | $[L T^{-1}]$ | Space $L$ $\leftrightarrow$ Time $T$: $L = cT$ |

- Neither mass, charge, nor temperature is involved. A constant of **pure geometry**.
- The foundation of special relativity and the sole constant that unifies spacetime.
- TOE question: Can it be derived from the geometrical structure of spacetime? Variable speed of light (VSL) theories hypothesize it may vary.

#### Group II: Spacetime + Mass ($M, L, T$)

| # | Parameter | Symbol | Value | Dimensions | Conversion |
|---|-----------|--------|-------|------------|------------|
| P2 | Planck constant | $\hbar$ | 1.055 $\times$ 10$^{-34}$ J$\cdot$s (definition) | $[M L^2 T^{-1}]$ ($M^{+1}$) | Energy $E$ $\leftrightarrow$ Frequency $\nu$: $E = \hbar\omega$ |
| P3 | Gravitational constant | G | 6.674 $\times$ 10$^{-11}$ m$^3$/(kg$\cdot$s$^2$) | $[L^3 M^{-1} T^{-2}]$ ($M^{-1}$) | Mass $M$ $\leftrightarrow$ Spacetime curvature $L^{-2}$: $R \sim GM/(c^2 r)$ |

- Both are composed of the same three dimensions $M, L, T$, but **the sign of $M$ is opposite**.
  - $\hbar$: $M^{+1}$ — as mass increases, energy (frequency) increases
  - $G$: $M^{-1}$ — as mass increases, the gravitational constant's contribution inversely decreases (the curvature radius grows)
- This **opposite behavior with respect to mass** already suggests at the level of dimensional analysis the fundamental reason why quantum mechanics and general relativity are difficult to unify.
- $\hbar$ is the foundation of quantum mechanics (discretization), $G$ is the foundation of general relativity (geometrization).
- $G$ has a **duality**: it is both a dimensional conversion constant and the gravitational coupling constant (strength of the force).

#### Group III: Spacetime + Mass + Temperature ($M, L, T, \Theta$)

| # | Parameter | Symbol | Value | Dimensions | Conversion |
|---|-----------|--------|-------|------------|------------|
| P4 | Boltzmann constant | $k_B$ | 1.381 $\times$ 10$^{-23}$ J/K (definition) | $[M L^2 T^{-2} \Theta^{-1}]$ | Temperature $\Theta$ $\leftrightarrow$ Energy $E$: $E = k_B T$ |

- Adds the temperature dimension $\Theta$ to Group II ($M, L, T$).
- The constant that bridges statistical mechanics (many-body, macroscopic) and microscopic physics.
- Bridges informational entropy and thermodynamic entropy. In an information-theoretic TOE, $k_B = 1$ is natural, and it may not be an independent parameter.

#### Group IV: Charge Only ($Q$)

| # | Parameter | Symbol | Value | Dimensions | Conversion |
|---|-----------|--------|-------|------------|------------|
| P5 | Elementary charge | e | 1.602 $\times$ 10$^{-19}$ C (definition) | $[Q]$ | Minimum unit of charge (standard of quantization) |

- Depends on none of $L, T, M, \Theta$. **Intrinsically independent of spacetime geometry.**
- Constitutes the electromagnetic coupling constant $\alpha = e^2/(4\pi\epsilon_0\hbar c)$.
- TOE question: Why is charge quantized? By Dirac's quantization condition, if magnetic monopoles exist, charge quantization becomes inevitable.

#### Group V: Charge + Spacetime + Mass ($M, L, T, Q$)

| # | Parameter | Symbol | Value | Dimensions | Conversion |
|---|-----------|--------|-------|------------|------------|
| P6 | Vacuum permeability/permittivity | $\mu_0$, $\epsilon_0$ | $c = 1/\sqrt{\mu_0 \epsilon_0}$ | $\mu_0$: $[M L Q^{-2}]$, $\epsilon_0$: $[Q^2 T^2 M^{-1} L^{-3}]$ | Connects charge $Q$ to spacetime $(L,T)$ and mass $M$ |

- Connects the charge $Q$ of Group IV to the spacetime and mass of Groups I–II.
- Dependent on $c$ and thus not independent, but plays the structural role of **determining how isolated charge connects to spacetime**.
- If the vacuum is not "empty" (has quantum structure), the meaning of these constants may change.

#### Dimensional Hierarchy

From the above classification, the fundamental dimensions have the following hierarchy:

$$\underbrace{L, T}_{\text{Spacetime (geometry)}} \xrightarrow{+M} \underbrace{M, L, T}_{\text{Matter (QM / gravity)}} \xrightarrow{+\Theta} \underbrace{M, L, T, \Theta}_{\text{Statistics (information)}} $$

$$\underbrace{Q}_{\text{Charge (isolated)}} \xrightarrow{\mu_0, \epsilon_0} \underbrace{M, L, T, Q}_{\text{Electromagnetism}}$$

In natural units ($c = \hbar = k_B = G = 1$), the distinctions among $L, T, M, \Theta$ are eliminated, and all physical quantities reduce to dimensionless pure numbers. The charge dimension $Q$ is separately eliminated by setting $e = 1$.

**The "numbers" that remain are the true arbitrary parameters that a Theory of Everything must explain.**

**Note:** The 2019 SI unit reform made c, $\hbar$, $k_B$, e all **exact values** (definitions). This does not mean "their values were discovered" but rather "the units were aligned to these constants." The physical origin of these constants remains unresolved.

### 1.1 Independence Analysis of the "19 Arbitrary Parameters"

Below we enumerate the 19 parameters conventionally called "arbitrary" in the Standard Model, plus gravity, cosmology, and neutrino parameters, and analyze whether each is **a truly independent parameter or a dependent quantity derivable from the fundamental constants of §1.0 or other parameters**.

#### A. Coupling Constants of the Forces (Conventionally 4)

**#1 Electromagnetic coupling constant $\alpha \approx 1/137.036$**

$$\alpha = \frac{e^2}{4\pi\epsilon_0\hbar c}$$

- **Independent? $\times$ Dependent.** Merely a combination of the fundamental constants P1(c), P2($\hbar$), P4(e), P6($\epsilon_0$) from §1.0.
- $\alpha \approx 1/137$ is called a "mysterious number," but substituting the values of e, $\epsilon_0$, $\hbar$, c automatically yields this value. There is nothing to derive.
- The true question is "why does the elementary charge $e$ have this value" — $\alpha$ is not the substance of the question.
- **Should not be counted as an arbitrary parameter. Double-counted with $e$ (P4).**

**#2 Strong force coupling constant $\alpha_s \approx 0.12$ (at 91 GeV)**

$$\alpha_s(Q^2) \approx \frac{12\pi}{(33-2N_f)\ln(Q^2/\Lambda_{QCD}^2)} \quad \text{(1-loop approximation)}$$

- **Independent? $\triangle$ Conditional.** The value of $\alpha_s$ itself depends on the energy scale $Q$ and is not a constant. The sole independent parameter is $\Lambda_{QCD} \approx 200$ MeV (QCD scale).
- However, this derivation formula is a 1-loop perturbative expansion approximation that breaks down at low energies (large distances). No exact solution has been found (Millennium Prize Problem).
- The spatial dependence is logarithmic, not $1/r^2$ (power law). Qualitatively different from gravity and electromagnetism.
- **$\Lambda_{QCD}$ is recognized as an independent parameter, but $\alpha_s$ itself is a dependent quantity.**

**#3 Weak force coupling constant $\alpha_w \approx 1/30$**

$$\alpha_w = \frac{g^2}{4\pi} , \quad G_F = \frac{\sqrt{2}}{8} \frac{g^2}{M_W^2} = \frac{1}{\sqrt{2}v^2}$$

- **Independent? $\triangle$ Conditional.** The Fermi constant $G_F$ is derived from the Higgs vacuum expectation value $v$ ($G_F = 1/(\sqrt{2}v^2)$).
- In the electroweak unification theory, the weak and electromagnetic forces share the same origin (SU(2)$\times$U(1) symmetry) and are mixed by the Weinberg angle $\theta_W$.
- **$v$ (Higgs VEV) and $\theta_W$ are the independent parameters. $\alpha_w$ itself is a dependent quantity.**

**#4 Gravitational coupling constant $\alpha_G \approx 5.9 \times 10^{-39}$**

$$\alpha_G = \frac{Gm_p^2}{\hbar c}$$

- **Independent? $\times$ Dependent.** A combination of P1(c), P2($\hbar$), P5(G) from §1.0 and the particle mass $m_p$.
- Moreover, "which mass $m$ to use" is arbitrary (proton? electron? Planck mass?). The definition itself is not unique.
- **Should not be counted as an arbitrary parameter. Double-counted with G (P5) and particle masses.**

**Summary of coupling constants:** Among the "four coupling constants," only **$\Lambda_{QCD}$** and **the Weinberg angle $\theta_W$** can be counted as independent parameters. $\alpha$ and $\alpha_G$ are combinations of the fundamental constants from §1.0 and particle masses, and are not independent.

#### B. Cosmological Constant (1)

**#5 Cosmological constant $\Lambda \approx 1.1 \times 10^{-52}$ m$^{-2}$**

- **Independent? $\checkmark$ Independent.** No method of deriving it from other parameters is known.
- The theoretical prediction from quantum field vacuum zero-point energy deviates from the observed value by 120 orders of magnitude.
- Its dimension is $[L^{-2}]$ (curvature). It determines the geometrical properties of spacetime.
- **A truly independent arbitrary parameter.**

#### C. Particle Masses (9)

In the Standard Model, particle masses are determined through Yukawa coupling constants $y_i$ with the Higgs field via $m_i = y_i \cdot v / \sqrt{2}$.

| # | Parameter | Mass | Yukawa coupling $y_i = m_i\sqrt{2}/v$ | Independent? |
|---|-----------|------|---------------------------------------|-------------|
| 6 | Electron $m_e$ | 0.511 MeV | $y_e \approx 2.9 \times 10^{-6}$ | $\triangle$ $y_e$ is independent. $m_e$ is $y_e \times v$ |
| 7 | Muon $m_\mu$ | 105.7 MeV | $y_\mu \approx 6.1 \times 10^{-4}$ | Same |
| 8 | Tau $m_\tau$ | 1776.9 MeV | $y_\tau \approx 1.0 \times 10^{-2}$ | Same |
| 9 | Up $m_u$ | $\approx$ 2.2 MeV | $y_u \approx 1.3 \times 10^{-5}$ | Same |
| 10 | Down $m_d$ | $\approx$ 4.7 MeV | $y_d \approx 2.7 \times 10^{-5}$ | Same |
| 11 | Charm $m_c$ | $\approx$ 1275 MeV | $y_c \approx 7.3 \times 10^{-3}$ | Same |
| 12 | Strange $m_s$ | $\approx$ 95 MeV | $y_s \approx 5.5 \times 10^{-4}$ | Same |
| 13 | Top $m_t$ | $\approx$ 173200 MeV | $y_t \approx 1.0$ (nearly 1) | Same |
| 14 | Bottom $m_b$ | $\approx$ 4180 MeV | $y_b \approx 2.4 \times 10^{-2}$ | Same |

- **Independent? $\triangle$ Conditional.** The masses $m_i$ themselves are dependent quantities derived via $y_i \times v/\sqrt{2}$. The truly independent quantities are the Yukawa coupling constants $y_i$ (9) and the Higgs VEV $v$ (1).
- Notably, $y_t \approx 1$ (top quark) is the only "natural" value, while the remaining 8 are unnaturally small at $10^{-6}$ to $10^{-2}$. **Why the Yukawa couplings exhibit such hierarchy is unsolved.**
- **The 9 masses decompose into 9 Yukawa coupling constants $y_i$ + $v$. They should not be double-counted as masses.**

#### D. CKM Matrix Parameters (4)

| # | Parameter | Value | Independent? |
|---|-----------|-------|-------------|
| 15 | Mixing angle $\theta_{12}$ | $\approx$ 13.0$^\circ$ (Cabibbo angle) | $\checkmark$ Independent |
| 16 | Mixing angle $\theta_{23}$ | $\approx$ 2.4$^\circ$ | $\checkmark$ Independent |
| 17 | Mixing angle $\theta_{13}$ | $\approx$ 0.2$^\circ$ | $\checkmark$ Independent |
| 18 | CP phase $\delta_{CKM}$ | $\approx$ 1.2 rad | $\checkmark$ Independent |

- **Independent? $\checkmark$ Independent.** No method of derivation from other parameters is known.
- However, if the generation structure of quarks (why 3 generations) is elucidated, they might be determined geometrically or group-theoretically.
- **Parameters specific to the weak force. Unnecessary if the weak force did not exist.**

#### E. Higgs Field Parameters (2)

**#19 Higgs boson mass $m_H \approx 125.1$ GeV**

$$m_H = \sqrt{2\lambda} \cdot v$$

- **Independent? $\triangle$ Conditional.** Derived from the Higgs self-coupling constant $\lambda$ and VEV $v$. The truly independent quantities are $\lambda$ and $v$.
- Treating $m_H$ as an "arbitrary parameter" merely packages $\lambda \times v$ into a single number.

**#20 Vacuum expectation value $v \approx 246.2$ GeV**

$$v = (\sqrt{2}G_F)^{-1/2} = \frac{2M_W}{g}$$

- **Independent? $\checkmark$ Independent.** Determines the scale of electroweak symmetry breaking. The reference for all particle mass scales.
- The range of the weak force (W/Z boson masses) is also determined by this value.
- **One of the most important independent parameters.**

#### F. QCD Vacuum Parameter (1)

**#21 QCD vacuum angle $\theta_{QCD} < 10^{-10}$**

- **Independent? $\checkmark$ Independent (but the value is nearly zero).**
- The experimental upper bound is $< 10^{-10}$, and it may be physically zero.
- If zero, some symmetry (Peccei-Quinn symmetry) underlies it, making it not an arbitrary parameter but **structurally zero**.
- **Whether it is an "arbitrary parameter" or "structurally zero" is undetermined.**

#### G. Neutrino Parameters (7+)

| # | Parameter | Value | Independent? |
|---|-----------|-------|-------------|
| 22 | $\Delta m^2_{21}$ | $\approx 7.5 \times 10^{-5}$ eV$^2$ | $\checkmark$ Independent (squared mass difference; absolute mass unmeasured) |
| 23 | $\Delta m^2_{32}$ | $\approx 2.5 \times 10^{-3}$ eV$^2$ | $\checkmark$ Independent |
| 24 | PMNS $\theta_{12}$ | $\approx$ 33.4$^\circ$ | $\checkmark$ Independent |
| 25 | PMNS $\theta_{23}$ | $\approx$ 49$^\circ$ | $\checkmark$ Independent |
| 26 | PMNS $\theta_{13}$ | $\approx$ 8.5$^\circ$ | $\checkmark$ Independent |
| 27 | PMNS $\delta$ | Undetermined | $\checkmark$ Independent (unmeasured) |
| 28 | Majorana phases $\alpha_1$, $\alpha_2$ | Unmeasured | ? Whether they exist is also undetermined |

### 1.8 Summary of Independence Analysis

Through the above analysis, the reality of the "19 arbitrary parameters" is as follows.

#### Dependent Quantities (Derivable from Others; Should Not Be Counted as Independent Parameters)

| Parameter | Derivation Formula | Double-Counted With |
|-----------|-------------------|-------------------|
| $\alpha \approx 1/137$ | $e^2/(4\pi\epsilon_0\hbar c)$ | P1(c), P2($\hbar$), P4(e), P6($\epsilon_0$) |
| $\alpha_G \approx 10^{-39}$ | $Gm^2/(\hbar c)$ | P1(c), P2($\hbar$), P5(G) + particle mass. Moreover, $m$ is arbitrary, making the definition non-unique |
| $\alpha_s$ (at specific scale) | $12\pi/((33-2N_f)\ln(Q^2/\Lambda_{QCD}^2))$ | Reduces to $\Lambda_{QCD}$ |
| $\alpha_w$ | $g^2/(4\pi)$ | Reduces to $v$ and $\theta_W$ |
| 9 particle masses | $m_i = y_i \cdot v/\sqrt{2}$ | Reduces to Yukawa couplings $y_i$ + $v$ |
| Higgs mass $m_H$ | $m_H = \sqrt{2\lambda} \cdot v$ | Reduces to Higgs self-coupling $\lambda$ + $v$ |

#### Truly Independent Parameters

| Category | Parameters | Count |
|----------|-----------|-------|
| Fundamental physical constants | c, $\hbar$, $k_B$, e, G, ($\mu_0\epsilon_0$) | 6 ($\mu_0\epsilon_0$ dependent on c) |
| Cosmological constant | $\Lambda$ | 1 |
| Strong force scale | $\Lambda_{QCD}$ | 1 |
| Electroweak mixing | Weinberg angle $\theta_W$ | 1 |
| Yukawa couplings | $y_e, y_\mu, y_\tau, y_u, y_d, y_c, y_s, y_t, y_b$ | 9 |
| Higgs | Self-coupling $\lambda$, vacuum expectation value $v$ | 2 |
| CKM matrix | $\theta_{12}, \theta_{23}, \theta_{13}, \delta_{CKM}$ | 4 |
| QCD vacuum | $\theta_{QCD}$ | 1 (possibly structurally zero) |
| Neutrino | $\Delta m^2_{21}, \Delta m^2_{32}, \theta_{12}, \theta_{23}, \theta_{13}, \delta, (\alpha_1, \alpha_2)$ | 7+ |
| **Total** | | **$\approx$32** (5 fundamental constants + $\approx$27 arbitrary parameters) |

#### Conclusion of This Section (Provisional)

The conventional count of "19 arbitrary parameters" is inaccurate in the following respects:

1. **Double-counting:** $\alpha$, $\alpha_G$, $\alpha_w$, particle masses, and $m_H$ are combinations of other parameters and should not be counted as independent.
2. **Under-counting:** Gravity (G, $\Lambda$), neutrinos (7+), and hidden parameters ($\Lambda_{QCD}$, $\theta_W$, Yukawa couplings, $\lambda$) are not included.
3. **Exclusion of fundamental constants:** c, $\hbar$, G, $k_B$, e are excluded as "constants," but the origin of their values is unresolved.

As a KPI for the Theory of Everything, we evaluate the following for these approximately 32 parameters:

$$\text{KPI} = \frac{\text{Number of parameters whose theoretically derived values match experiment}}{\text{Total number of parameters under verification ($\approx$32)}}$$

The ideal is to reproduce all $\approx$32 from minimal input (0–2). The current Standard Model's KPI is effectively zero (it merely returns its inputs).

---

## 2. Structural Analysis of the Coupling Constants

As shown in §1.1, the "four coupling constants" are all dependent quantities, but the structure of their derivation formulas conceals important problems.

### 2.1 Electromagnetic Coupling Constant $\alpha$

$$\alpha = \frac{e^2}{4\pi\epsilon_0\hbar c} \approx \frac{1}{137.036}$$

**Dimensional analysis of constituents:**

| Element | Dimensions | §1.0 Classification |
|---------|------------|---------------------|
| $e^2$ | $[Q^2]$ | Group IV (Charge) |
| $\epsilon_0$ | $[Q^2 T^2 M^{-1} L^{-3}]$ | Group V (Charge $\leftrightarrow$ Spacetime) |
| $\hbar$ | $[M L^2 T^{-1}]$ | Group II (Spacetime + Mass) |
| $c$ | $[L T^{-1}]$ | Group I (Pure Spacetime) |

Dimensional verification: $[Q^2] / ([Q^2 T^2 M^{-1} L^{-3}] \cdot [M L^2 T^{-1}] \cdot [L T^{-1}]) = [1]$ (dimensionless) $\checkmark$

**Structural meaning:** $\alpha$ represents the ratio of how charge (Group IV) couples with spacetime + mass (Groups I, II, V). It uses 4 of the 5 fundamental constants (c, $\hbar$, e, $\epsilon_0$) to produce a single dimensionless number. **There is no mystery in the value of $\alpha$; if mystery exists, it lies in the value of $e$.**

### 2.2 Strong Force Coupling Constant $\alpha_s$

$$\alpha_s(Q^2) \approx \frac{12\pi}{(33-2N_f)\ln(Q^2/\Lambda_{QCD}^2)} \quad \text{(1-loop approximation)}$$

**Structural problems:**

1. **Not a constant.** A function of the energy scale $Q$ ($\approx$ inverse of spatial resolution $\hbar/r$). The term "coupling constant" is misleading.
2. **Approximate formula.** A result of 1-loop perturbative expansion, not an exact solution. Higher-order corrections increase complexity.
3. **Breaks down at low energy.** As $Q \to \Lambda_{QCD}$, the logarithm diverges and computation becomes impossible. The confinement regime relies on lattice QCD (numerical computation). The existence of an exact solution is unproven (Millennium Prize Problem).
4. **Spatial dependence is not geometrical.** Compared to the $1/r^2$ of gravity and electromagnetism (inverse of area in 3D space), $\alpha_s$ is logarithmic. This does not arise naturally from dimensional analysis but appears only as a result of quantum corrections (virtual particle pair creation and annihilation).

**The only independent parameter is $\Lambda_{QCD}$.** However, $\Lambda_{QCD}$ plays the role of a "hidden dimensional conversion constant that connects color charge (a discrete label) to the spatial dimension $L$" (unclassified in §1.0).

### 2.3 Weak Force Coupling Constant $\alpha_w$

$$\alpha_w = \frac{g^2}{4\pi}, \quad \sin^2\theta_W = 1 - \frac{M_W^2}{M_Z^2} \approx 0.231$$

$$G_F = \frac{1}{\sqrt{2}v^2} \approx 1.166 \times 10^{-5} \text{ GeV}^{-2}$$

**Structural features:**

1. **Already unified with electromagnetism (electroweak unification).** In SU(2)$\times$U(1) gauge theory, the weak and electromagnetic forces are mixed by the Weinberg angle $\theta_W$. At low energies they appear as separate forces, but above $\approx$ 100 GeV they are unified.
2. **Finite range.** Due to the masses of the W/Z bosons ($M_W \approx 80.4$ GeV, $M_Z \approx 91.2$ GeV), the range of the force is limited to $\sim \hbar/(M_W c) \approx 10^{-18}$ m. It follows a Yukawa-type potential $e^{-M_W r}/r$.
3. **Reduces to the Higgs VEV $v$.** Since $G_F = 1/(\sqrt{2}v^2)$, the scale of the weak force is entirely determined by $v$.

**The independent parameters are $v$ and $\theta_W$.** $\alpha_w$ itself is a dependent quantity.

### 2.4 Gravitational Coupling Constant $\alpha_G$

$$\alpha_G = \frac{Gm^2}{\hbar c}$$

**Structural problems:**

1. **The choice of mass $m$ is arbitrary.** The other three coupling constants are determined solely by properties of the force, but $\alpha_G$ changes depending on "which particle's mass is used."

   | Mass choice | Value of $\alpha_G$ |
   |-------------|-------------------|
   | Electron mass $m_e$ | $\approx 1.8 \times 10^{-45}$ |
   | Proton mass $m_p$ | $\approx 5.9 \times 10^{-39}$ |
   | Planck mass $m_{Pl}$ | $= 1$ (by definition) |

   A "constant" whose definition is not even unique should not be called a constant.

2. **A combination of §1.0 constants.** G (P5), $\hbar$ (P2), c (P1) are all fundamental constants. $\alpha_G$ contains no new information and is merely a double-counting of particle masses and fundamental constants.

3. **"Gravity is weak" is a question of curvature radius.** The smallness of $\alpha_G$ does not mean "gravity is weak" but reflects that the curvature radius of spacetime is overwhelmingly larger than the curvature radius of electromagnetic acceleration. Comparing them at the same scale is itself a category error.

### 2.5 Fundamental Problem in Comparing the Four Coupling Constants

Conventional physics lines up the four coupling constants as follows and poses the "hierarchy problem":

| Force | Coupling constant | Relative ratio |
|-------|------------------|---------------|
| Strong | $\alpha_s \approx 1$ | 1 |
| Electromagnetic | $\alpha \approx 1/137$ | $10^{-2}$ |
| Weak | $\alpha_w \approx 1/30$ | $10^{-2}$ |
| Gravity | $\alpha_G \approx 10^{-39}$ | $10^{-39}$ |

However, our analysis reveals that this comparison is invalid for the following reasons:

1. **The objects of comparison are not of the same kind.** $\alpha$ is a combination of fundamental constants, $\alpha_s$ is a scale-dependent function, $\alpha_w$ is a dependent quantity of $v$ and $\theta_W$, and $\alpha_G$ is an indeterminate quantity depending on mass selection. Placing these in the same table is itself meaningless.
2. **Characteristic scales differ.** Each force operates at its own spatial scale ($10^{-18}$ m to cosmic scale), and there is no basis for comparing them at a single energy scale.
3. **Spatial dependencies are qualitatively different.** Gravity and electromagnetism follow $1/r^2$ (geometrical), the strong force follows a logarithmic function (quantum corrections), and the weak force follows exponential decay (Yukawa type). Comparing the "strength" of quantities with different mathematical structures by a single number is physically meaningless.
4. **The "hierarchy problem" is a pseudo-problem.** From points 1–3 above, the question "why is gravity alone 39 orders of magnitude weaker?" is itself improperly formulated.

---

## 3. Structural Analysis of Particle Masses

### 3.1 Hierarchy of Yukawa Coupling Constants

As shown in §1.1 C, particle masses $m_i$ are products of the Yukawa coupling constants $y_i$ and the Higgs VEV $v$.

$$m_i = \frac{y_i \cdot v}{\sqrt{2}}$$

The truly independent quantities are $y_i$, whose hierarchy is as follows:

| Generation | Lepton | $y$ | Quark (up-type) | $y$ | Quark (down-type) | $y$ |
|------------|--------|-----|----------------|-----|-------------------|-----|
| 1st | $e$ | 2.9 $\times$ 10$^{-6}$ | $u$ | 1.3 $\times$ 10$^{-5}$ | $d$ | 2.7 $\times$ 10$^{-5}$ |
| 2nd | $\mu$ | 6.1 $\times$ 10$^{-4}$ | $c$ | 7.3 $\times$ 10$^{-3}$ | $s$ | 5.5 $\times$ 10$^{-4}$ |
| 3rd | $\tau$ | 1.0 $\times$ 10$^{-2}$ | $t$ | **$\approx$ 1** | $b$ | 2.4 $\times$ 10$^{-2}$ |

**Notable structures:**

1. **Only $y_t \approx 1$ (top quark) is "natural."** The remaining 8 are unnaturally small at $10^{-6}$ to $10^{-2}$. Why $y_t$ alone matches the Higgs VEV scale is unsolved.
2. **Regularity across generations.** Within the same column (leptons, up-type quarks, down-type quarks), the values roughly increase by 2–3 orders of magnitude per generation. Whether this is coincidence or reflects some structure is unknown.
3. **Are all 9 independent?** If the inter-generation ratios have regularity, the independent parameters might be reducible to "reference mass + ratio rule."

### 3.2 The Problem of Mass Origin

The Standard Model states that "the Higgs mechanism gives mass," but more precisely, the Higgs mechanism provides the **mechanism** of mass (spontaneous symmetry breaking) without explaining the **values** of masses (Yukawa coupling constants $y_i$). Why the 9 $y_i$ have these values is in principle unanswerable within the framework of the Standard Model.

---

## 4. Structural Analysis of Mixing Angles and CP Phases

### 4.1 CKM Matrix (Quarks)

The CKM matrix is the rotation matrix between the weak force gauge eigenstates and mass eigenstates. A 3-generation unitary matrix is completely described by 3 rotation angles and 1 phase.

| Parameter | Value | Physical meaning |
|-----------|-------|-----------------|
| $\theta_{12} \approx 13.0^\circ$ | Cabibbo angle | Mixing between 1st–2nd generations. Largest |
| $\theta_{23} \approx 2.4^\circ$ | | Between 2nd–3rd generations. Small |
| $\theta_{13} \approx 0.2^\circ$ | | Between 1st–3rd generations. Extremely small |
| $\delta_{CKM} \approx 1.2$ rad | CP phase | Matter–antimatter asymmetry |

**Structural problems:**

1. **Why 3 generations?** The CKM matrix is 3$\times$3 because quarks have 3 generations, but why there are 3 generations is unresolved. With 2 generations the CP phase would be unnecessary; with 4+ generations additional parameters would be required.
2. **Hierarchy of mixing angles.** There is a hierarchy $\theta_{12} > \theta_{23} > \theta_{13}$. This may correlate with the Yukawa coupling hierarchy (inter-generation mass ratios).
3. **Specific to the weak force.** The CKM matrix would be unnecessary without the weak force. If the existence of the weak force itself is derived from structural requirements, the mixing angles may also be structurally determined.

### 4.2 PMNS Matrix (Neutrinos)

The neutrino version of the CKM matrix. Same structure but vastly different values.

| Parameter | CKM (Quarks) | PMNS (Neutrinos) |
|-----------|-------------|-------------------|
| $\theta_{12}$ | 13.0$^\circ$ (small) | 33.4$^\circ$ (large) |
| $\theta_{23}$ | 2.4$^\circ$ (small) | 49$^\circ$ (nearly maximal mixing) |
| $\theta_{13}$ | 0.2$^\circ$ (tiny) | 8.5$^\circ$ (small but finite) |

**Striking contrast between CKM and PMNS:** Quark mixing is small (diagonal components dominate) while neutrino mixing is large (off-diagonal components are prominent). Why the values differ so greatly despite the same mathematical structure (3$\times$3 unitary matrix) is unsolved.

---

## 5. Structural Analysis of the Higgs Field and Cosmological Constant

### 5.1 Higgs Field

**Independent parameters:** Self-coupling constant $\lambda \approx 0.13$, vacuum expectation value $v \approx 246.2$ GeV

$$V(\phi) = -\mu^2|\phi|^2 + \lambda|\phi|^4, \quad v = \frac{\mu}{\sqrt{\lambda}}, \quad m_H = \sqrt{2\lambda}\cdot v$$

**Structural role:**

- $v$ determines the scale of electroweak symmetry breaking and is the most important parameter as the reference for all particle masses.
- $\lambda$ determines the mass of the Higgs boson itself, but there is no theoretical basis for the value $\approx 0.13$.
- The origins of $v$ and $\lambda$ cannot be explained within the framework of the Standard Model. "Why symmetry breaks" is merely assumed structurally.

### 5.2 Cosmological Constant $\Lambda$

$$G_{\mu\nu} + \Lambda g_{\mu\nu} = \frac{8\pi G}{c^4} T_{\mu\nu}$$

**Independent parameter:** $\Lambda \approx 1.1 \times 10^{-52}$ m$^{-2}$

#### Derivation of the Observed Value (Geometrical Method)

In the Einstein equation, using the energy density parameter $\Omega_\Lambda \approx 0.685$ and the Hubble constant $H_0 \approx 67.4$ km/s/Mpc obtained from observations of the accelerating expansion of the universe (Type Ia supernovae, CMB, BAO),

$$\Lambda = 3 \Omega_\Lambda \frac{H_0^2}{c^2} \approx 3 \times 0.685 \times \frac{(2.18 \times 10^{-18} \; \text{s}^{-1})^2}{(3.0 \times 10^8 \; \text{m/s})^2} \approx 1.1 \times 10^{-52} \; \text{m}^{-2}$$

However, this derivation is **circular**. $H_0$ and $\Omega_\Lambda$ are values obtained from the observed curvature of the universe, and using them to determine $\Lambda$ is merely substitution. We are not "deriving" the value of $\Lambda$ from geometry but **merely assigning the name $\Lambda$ to the observed curvature**.

#### Quantum Field Theory Theoretical Prediction

In quantum field theory, the vacuum possesses a finite energy density as the sum of zero-point energies of all field modes. Integrating the zero-point energy $\frac{1}{2}\hbar\omega$ of each mode up to an ultraviolet cutoff $k_{\max}$, the vacuum energy density is

$$\rho_{\text{vac}} = \frac{1}{(2\pi)^3} \int_0^{k_{\max}} \frac{1}{2} \hbar \omega_k \cdot 4\pi k^2 \, dk = \frac{\hbar}{4\pi^2 c} \int_0^{k_{\max}} \frac{1}{2} \sqrt{k^2 c^2 + m^2 c^4 / \hbar^2} \cdot k^2 \, dk$$

For a massless field ($m = 0$) with cutoff at the Planck scale $k_{\max} = k_P = \sqrt{\frac{c^3}{\hbar G}}$,

$$\rho_{\text{vac}} \sim \frac{\hbar c}{16\pi^2} k_P^4 = \frac{c^7}{16\pi^2 \hbar G^2}$$

Converting to a cosmological constant ($\Lambda_{\text{QFT}} = \frac{8\pi G}{c^4} \rho_{\text{vac}}$),

$$\Lambda_{\text{QFT}} \sim \frac{c^3}{2\pi \hbar G} = \frac{1}{2\pi \ell_P^2}$$

where $\ell_P = \sqrt{\hbar G / c^3} \approx 1.616 \times 10^{-35}$ m is the Planck length. Therefore,

$$\Lambda_{\text{QFT}} \sim \frac{1}{2\pi \times (1.616 \times 10^{-35})^2} \approx 6.1 \times 10^{68} \; \text{m}^{-2}$$

#### Discrepancy

| | Value (m$^{-2}$) | Derivation method |
|---|-----------|---------|
| Observed value | $1.1 \times 10^{-52}$ | Observed data of accelerating cosmic expansion (circular) |
| QFT theoretical value | $\sim 6 \times 10^{68}$ | Planck-scale integration of vacuum zero-point energy |
| **Discrepancy** | **$\approx$120 orders of magnitude** | |

**Structural problems:**

1. **The 120 orders of magnitude problem.** As shown above, the theoretical prediction from quantum field theory vacuum energy deviates from the observed value by approximately 120 orders of magnitude. The worst theoretical prediction in all of physics. The essence of this discrepancy lies in quantum field theory's assumption that "all modes possess zero-point energy up to the Planck scale." Whether this assumption itself is correct has not been verified.
2. **Dimension is curvature.** $[\Lambda] = [L^{-2}]$, representing intrinsic curvature of spacetime. Same geometrical category as G ($M$ $\leftrightarrow$ curvature conversion).
3. **Why not zero?** Many physicists considered $\Lambda = 0$ to be "natural," but the accelerating expansion of the universe was discovered in 1998, confirming $\Lambda > 0$.
4. **Cutoff dependence.** The QFT theoretical value is proportional to the 4th power of the cutoff $k_{\max}$. There is no physical basis for choosing the Planck scale, and changing $k_{\max}$ changes the theoretical value arbitrarily. That is, the 120-order discrepancy does not mean "the theoretical prediction was wrong" but indicates that **the formulation itself is incomplete**.

### 5.3 QCD Vacuum Angle $\theta_{QCD}$

**Independent parameter:** $\theta_{QCD} < 10^{-10}$

Experimentally nearly zero, but there is no theoretical necessity for it to be zero. Why it is near zero is unsolved as the "strong CP problem." If a Peccei-Quinn symmetry exists, it would be structurally zero, in which case the axion particle is predicted (not yet discovered).

---

## 6. Classification Results

### 6.1 Comparison of Formulation Levels of the Four Forces

Our analysis reveals that the four forces differ greatly in their theoretical completeness.

| Item | Gravity | Electromagnetism | Weak force | Strong force |
|------|---------|-----------------|------------|-------------|
| Spatial dependence | $1/r^2$ (exact) | $1/r^2$ (exact) | $e^{-Mr}/r$ (Yukawa type) | Logarithmic (1-loop approx.) |
| Computable at all distances? | $\checkmark$ | $\checkmark$ | $\checkmark$ | $\times$ (breaks down at low energy) |
| Exact solution exists? | $\checkmark$ | $\checkmark$ | $\checkmark$ | $\times$ (Millennium Prize Problem) |
| Dimension of charge | $M$ (continuous) | $Q$ (continuous, quantized) | Discrete ($\pm 1/2$) | Discrete (3 colors) |
| Discretization of charge explained? | Not addressed | Unsolved (Dirac condition) | Assumed (SU(2)) | Assumed (SU(3)) |
| Dimensional conversion constant | G (explicit) | e, $\mu_0\epsilon_0$ (explicit) | $v$ (treated separately) | $\Lambda_{QCD}$ (hidden) |
| Number of independent parameters | G, $\Lambda$ (2) | e (1) | $v$, $\theta_W$ (2) | $\Lambda_{QCD}$ (1) |

### 6.2 Parameter Dependency Relationships

Below we show the dependency relationships among the approximately 32 parameters. Arrows mean "B is derived from A."

```
§1.0 Fundamental Constants (Independent)
├── c (P1) ─────────────────┐
├── ℏ (P2) ──────────────┐  │
├── k_B (P3) ────────    │  │
├── e (P4) ──────────┐   │  │
├── G (P5) ──────┐   │   │  │
└── μ₀ε₀ (P6) ──┤   │   │  │
                 │   │   │  │
             ┌───┘   │   │  │
             ▼       ▼   ▼  ▼
          α_G ←── m_i   α = e²/(4πε₀ℏc)  ← Dependent (not independent)
                  ▲
                  │
          y_i × v/√2 = m_i
          ▲       ▲
          │       │
    Yukawa        Higgs VEV
    y_i (9)       v (1)  ← Independent parameters
                  │
                  ├── α_w = f(v, θ_W)  ← Dependent
                  ├── G_F = 1/(√2 v²)  ← Dependent
                  └── m_H = √(2λ)·v    ← Dependent
                       ▲
                       │
                  λ (1)  ← Independent parameter

Other Independent Parameters:
├── Λ_QCD (1) ── → α_s(Q²) = f(Λ_QCD, Q)  ← Dependent
├── θ_W (1)
├── Λ (1) ── Cosmological constant
├── θ_QCD (1)
├── CKM: θ₁₂, θ₂₃, θ₁₃, δ (4)
└── PMNS: Δm²₂₁, Δm²₃₂, θ₁₂, θ₂₃, θ₁₃, δ, (α₁, α₂) (7+)
```

### 6.3 TOE KPI (Acceptance Test Specification)

We define the evaluation criteria for a Theory of Everything as follows.

**Test cases:** Approximately 32 independent parameters

| ID | Category | Parameters | Count |
|----|----------|-----------|-------|
| T01–T05 | Fundamental constants | c, $\hbar$, $k_B$, e, G | 5 |
| T06 | Cosmological constant | $\Lambda$ | 1 |
| T07 | QCD scale | $\Lambda_{QCD}$ | 1 |
| T08 | Electroweak mixing | $\theta_W$ | 1 |
| T09–T17 | Yukawa couplings | $y_e, y_\mu, y_\tau, y_u, y_d, y_c, y_s, y_t, y_b$ | 9 |
| T18–T19 | Higgs | $\lambda$, $v$ | 2 |
| T20–T23 | CKM matrix | $\theta_{12}, \theta_{23}, \theta_{13}, \delta_{CKM}$ | 4 |
| T24 | QCD vacuum | $\theta_{QCD}$ | 1 |
| T25–T31 | Neutrino | $\Delta m^2_{21}, \Delta m^2_{32}, \theta_{12}, \theta_{23}, \theta_{13}, \delta_{PMNS}, (\alpha_1, \alpha_2)$ | 7+ |
| | **Total** | | **$\approx$32** |

**Evaluation method:**

$$\text{KPI} = \frac{\text{Number of values derived by the theory matching experiment within significant digits}}{32}$$

**Evaluation scale:**

| Score | Assessment |
|-------|-----------|
| 0/32 | Current Standard Model (merely returns inputs) |
| 1–5/32 | Partial success (e.g., explaining $\theta_{QCD} = 0$, unification of coupling constants) |
| 6–15/32 | Significant progress (e.g., explaining Yukawa coupling hierarchy) |
| 16–25/32 | Grand Unified Theory level |
| 26–31/32 | Theory of Everything candidate |
| **32/32** | **Completion of the Theory of Everything** |

**Bonus by number of input parameters:**

| Inputs | Assessment |
|--------|-----------|
| 32 (all inputs) | Nothing derived (current status) |
| 10–31 | Partial reduction |
| 3–9 | Significant reduction |
| 1–2 | Excellent TOE candidate |
| **0** | **Ultimate TOE (no initial conditions needed)** |

---

## 7. Conclusions

### 7.1 Reassessment of the "19 Arbitrary Parameters"

Our analysis reveals the following:

1. **The conventional "19" is inaccurate.** Both double-counting ($\alpha$, $\alpha_G$, $\alpha_w$, particle masses, $m_H$) and under-counting (G, $\Lambda$, $\Lambda_{QCD}$, $\theta_W$, Yukawa couplings, $\lambda$, neutrinos) exist. The accurate count is approximately 32 independent parameters.

2. **Exclusion of fundamental physical constants is unjustified.** c, $\hbar$, $k_B$, e, G are excluded as "constants," but the physical origin of their values is unresolved and they should be included as TOE verification targets.

3. **The four forces are not formulated at equal levels.** Gravity and electromagnetism possess exact analytical solutions ($1/r^2$) and explicit dimensional conversion constants, while the strong force has only approximate solutions (logarithmic) that break down at low energies, and the weak force hides its dimensional conversion constants within the Higgs field.

4. **There is no mystery in $\alpha \approx 1/137$.** It is a dependent quantity obtained automatically by substituting values into $e^2/(4\pi\epsilon_0\hbar c)$, and the assessment that "it would be amazing to derive it" reflects an incorrectly posed question. What should be asked is the origin of the elementary charge $e$.

5. **The "hierarchy problem" may be a pseudo-problem.** The four coupling constants differ in definition, spatial dependence, and charge dimensions. Comparing them at a single scale and lamenting that "gravity alone is weak" is the result of forcing quantities of different categories onto a common standard.

### 7.2 TOE KPI — Priority by Verifiability

Verifying all approximately 32 independent parameters at once is impossible at the current stage. The reason is that the theoretical maturity of the four forces is fundamentally different (see §6.1).

#### Verifiable Domain: Unification of Gravity and Electromagnetism

Gravity and electromagnetism satisfy the conditions necessary for verification in the following respects:

- Both possess exact analytical solutions ($1/r^2$)
- The dimensions of charge are explicitly defined (mass $M$, charge $Q$)
- Dimensional conversion constants are explicit ($G$, $e$, $\mu_0\epsilon_0$)
- The independent parameters involved are identified ($c$, $\hbar$, $G$, $e$, $\mu_0\epsilon_0$, $\Lambda$)

Therefore, the parameter relationships that a unified theory should derive are clear, and pass/fail judgment through comparison with experimental values is possible.

**KPIs for Gravity + Electromagnetism (Verifiable):**

| ID | Verification item | Pass criterion |
|----|------------------|---------------|
| T01 | Relationship between $G$ and $e$ | Derive one from the other and match experimental value |
| T02 | Value of $\alpha = e^2/(4\pi\epsilon_0\hbar c)$ | Derive $e$ from structure and match $\alpha \approx 1/137.036$ |
| T03 | Value of $\Lambda$ | Derive $\Lambda \approx 1.1 \times 10^{-52}$ m$^{-2}$ from theory |
| T04 | Relationship between mass $M$ and charge $Q$ | Derive the charge quantization condition |

#### Domain Difficult to Verify at Present: Weak and Strong Forces

For the weak and strong forces, the verification criterion for "having been unified" is itself unclear.

- **Strong force:** No exact solution exists (Millennium Prize Problem). Perturbation theory breaks down at low energy, and the spatial dependence of the coupling constant $\alpha_s$ is described only in 1-loop approximation. The "correct answer" against which to compare the output of a unified theory exists only as an approximation.
- **Weak force:** The charge (weak isospin) is a discrete label with no SI dimension. The dimensional conversion constant is hidden in the Higgs field vacuum expectation value $v$, and the separation of independent parameters is itself ambiguous.

Claiming to have "unified" incomplete theories merely means **bundling unformulated things while they remain unformulated**, and there is no way to verify it.

**Conclusion:** The first verification step for a Theory of Everything is the construction of a unified theory of gravity and electromagnetism, and the verification of T01–T04 above. Unification of the weak and strong forces becomes a verification target only after the rigorous formulation of each theory individually is completed.

### 7.3 Limitations of This Paper

This paper has classified existing parameters and organized the problems, but has not derived the parameters themselves. Structural derivation of each parameter becomes possible only at the stage of implementing the structural requirements formulated in the prior work [1] into concrete geometry.

---

## References

[1] N. Kihara, "Formulation of Structural Requirements That a Theory of Everything Must Satisfy," 2026. DOI: [10.5281/zenodo.19601592](https://doi.org/10.5281/zenodo.19601592)

[2] Particle Data Group, "Review of Particle Physics," *Phys. Rev. D* **110**, 030001 (2024).

[3] S. Weinberg, "A Model of Leptons," *Phys. Rev. Lett.* **19**, 1264 (1967).

[4] P. W. Higgs, "Broken Symmetries and the Masses of Gauge Bosons," *Phys. Rev. Lett.* **13**, 508 (1964).

[5] G. 't Hooft, "Naturalness, Chiral Symmetry, and Spontaneous Chiral Symmetry Breaking," *NATO Adv. Study Inst. Ser. B Phys.* **59**, 135 (1980).

[6] A. Einstein, "Die Feldgleichungen der Gravitation," *Sitzungsberichte der Preussischen Akademie der Wissenschaften zu Berlin*, 844–847 (1915).

[7] D. J. Gross and F. Wilczek, "Ultraviolet Behavior of Non-Abelian Gauge Theories," *Phys. Rev. Lett.* **30**, 1343 (1973).

[8] Planck Collaboration, "Planck 2018 results. VI. Cosmological parameters," *Astron. Astrophys.* **641**, A6 (2020).
